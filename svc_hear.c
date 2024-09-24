#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <float.h>

#include <gsl/gsl_histogram.h>
#include <gsl/gsl_statistics.h>

#define ALU               5.29177210903e-11 // SI: m
#define RAMTOAMU          1822.888485332
#define Hartree           4.3597447222071e-18 // SI: J 
#define Boltzmann         1.380649e-23 // SI: J * K^(-1)
#define Boltzmann_Hartree Boltzmann/Hartree // a.u. 
#define AVOGADRO          6.022140857 * 1e23 // mol^(-1)
#define HTOCM             2.1947463136320e5  // 1 Hartree in cm-1

#define MT_GENERATE_CODE_IN_HEADER 0
#include "mtwist.h"
/*
 * double mt_drand(void)
 *   Return a pseudorandom double in [0,1) with 32 bits of randomness
 *
 * uint32_t mt_lrand(void);
 *   Generate 32-bit random value 
 *
 */
double generate_normal(double sigma) 
/*
 * Generate normally distributed variable using Box-Muller method
 */
{
    double U = mt_drand();
    double V = mt_drand();
    return sigma * sqrt(-2 * log(U)) * cos(2.0 * M_PI * V);
}

#define ARENA_IMPLEMENTATION
#include "arena.h"
static Arena arena = {0};

#define COMMON_IMPLEMENTATION
#include "common.h"

#define HEAR_IMPLEMENTATION
#include "V_HeAr.h"

#include "raylib.h"

#define COLOR_BACKGROUND GetColor(0x181818FF)

#define FONT_SIZE_LOAD 160 
Font font = {0};

// ---------------------------------------------------------
#define m_He (4.00260325413 * RAMTOAMU)
#define m_Ar (39.9623831237 * RAMTOAMU)
#define mu m_He * m_Ar / (m_He + m_Ar)  // a.u.
                                         
#define lam 1.0/(2.0*mu) // hbar^2/2m

#define RMAX_COLLECT 30.0 // a.u.
#define RMAX_CUBE    35.0 // a.u.
//#define TT 300.0 // K
// ---------------------------------------------------------

#define XC(coords, i) coords[3*i + 0]
#define YC(coords, i) coords[3*i + 1]
#define ZC(coords, i) coords[3*i + 2]

Path path = {0};

typedef struct {
    double *items;
    size_t count;
    size_t capacity;
} SVC_Trace;

typedef struct {
    double *items;
    size_t count;
    size_t capacity;
} PositionTrace;

typedef struct {
    SVC_Trace svcs;
    PositionTrace positions;
} PIMC_Trace;

PIMC_Trace trace = {0};

void alloc_beads(Path *path)
{
    path->beads = (double**) malloc(path->numTimeSlices * sizeof(double*));
    for (size_t i = 0; i < path->numTimeSlices; ++i) {
        path->beads[i] = (double*) malloc(path->numParticles * 3*sizeof(double));
        memset(path->beads[i], 0.0, path->numParticles*3*sizeof(double));
    }
}

double PotentialAction(Path path, size_t tslice) 
{
    assert(tslice < path.numTimeSlices);

    double pot = 0.0;
    for (size_t ptcl = 0; ptcl < path.numParticles; ++ptcl) {
        double r = XC(path.beads[tslice], ptcl)*XC(path.beads[tslice], ptcl) + 
                   YC(path.beads[tslice], ptcl)*YC(path.beads[tslice], ptcl) + 
                   ZC(path.beads[tslice], ptcl)*ZC(path.beads[tslice], ptcl); 
        r = sqrt(r);
        pot = pot + V_HeAr(r);
    }

    return path.tau * pot;
}

double PotentialEnergy(Path path) {
    double result = 0.0;

    for (size_t tslice = 0; tslice < path.numTimeSlices; ++tslice) {
        for (size_t ptcl = 0; ptcl < path.numParticles; ++ptcl) {
            double r = XC(path.beads[tslice], ptcl)*XC(path.beads[tslice], ptcl) + 
                       YC(path.beads[tslice], ptcl)*YC(path.beads[tslice], ptcl) + 
                       ZC(path.beads[tslice], ptcl)*ZC(path.beads[tslice], ptcl); 
            r = sqrt(r);
            result = result + V_HeAr(r);
        }
    }

    return result / path.numTimeSlices;
}

double KineticEnergy(Path path) 
{
    double tot = 0.0;
    double norm = 1.0/(4.0*lam * path.tau*path.tau);

    for (size_t tslice = 0; tslice < path.numTimeSlices; ++tslice) {
        size_t tslicep1 = (tslice + 1) % path.numTimeSlices;

        for (size_t ptcl = 0; ptcl < path.numParticles; ++ptcl) {
            double dx = XC(path.beads[tslicep1], ptcl) - XC(path.beads[tslice], ptcl);
            double dy = YC(path.beads[tslicep1], ptcl) - YC(path.beads[tslice], ptcl);
            double dz = ZC(path.beads[tslicep1], ptcl) - ZC(path.beads[tslice], ptcl);
            tot = tot - norm*(dx*dx + dy*dy + dz*dz);
        }
    }

    // @NOTE:
    // note the 3.0/2.0 factor in the constant part of the energy estimator
    return 3.0/2.0*path.numParticles/path.tau + tot/path.numTimeSlices;
}

int COM_Move(Path path, size_t ptcl)
{
    assert(ptcl < path.numParticles);
    int result = 0;

    double delta = 0.1;
    double shiftx = delta*RMAX_CUBE*(-1.0 + 2.0*mt_drand());
    double shifty = delta*RMAX_CUBE*(-1.0 + 2.0*mt_drand());
    double shiftz = delta*RMAX_CUBE*(-1.0 + 2.0*mt_drand());


    double oldAction = 0.0;
    for (size_t tslice = 0; tslice < path.numTimeSlices; ++tslice) {
        oldAction = oldAction + PotentialAction(path, tslice);
    } 
    
    double *oldbeadsx = (double*) arena_alloc(&arena, path.numTimeSlices * sizeof(double));
    memset(oldbeadsx, 0.0, path.numTimeSlices * sizeof(double));
    double *oldbeadsy = (double*) arena_alloc(&arena, path.numTimeSlices * sizeof(double));
    memset(oldbeadsy, 0.0, path.numTimeSlices * sizeof(double));
    double *oldbeadsz = (double*) arena_alloc(&arena, path.numTimeSlices * sizeof(double));
    memset(oldbeadsz, 0.0, path.numTimeSlices * sizeof(double));

    for (size_t tslice = 0; tslice < path.numTimeSlices; ++tslice) {
        oldbeadsx[tslice] = XC(path.beads[tslice], ptcl);
        oldbeadsy[tslice] = YC(path.beads[tslice], ptcl);
        oldbeadsz[tslice] = ZC(path.beads[tslice], ptcl);
        XC(path.beads[tslice], ptcl) = XC(path.beads[tslice], ptcl) + shiftx;
        YC(path.beads[tslice], ptcl) = YC(path.beads[tslice], ptcl) + shifty;
        ZC(path.beads[tslice], ptcl) = ZC(path.beads[tslice], ptcl) + shiftz;
    }
    
    // reject move that moves COM in such that any particles goes beyond RMAX_CUBE
    bool within_cube = true;

    for (size_t tslice = 0; tslice < path.numTimeSlices; ++tslice) {
        double r = XC(path.beads[tslice], ptcl)*XC(path.beads[tslice], ptcl) + 
                   YC(path.beads[tslice], ptcl)*YC(path.beads[tslice], ptcl) + 
                   ZC(path.beads[tslice], ptcl)*ZC(path.beads[tslice], ptcl); 
        r = sqrt(r);

        if (r > RMAX_CUBE) within_cube = false;
    }

    if (!within_cube) 
    { 
        for (size_t tslice = 0; tslice < path.numTimeSlices; ++tslice) {
            XC(path.beads[tslice], ptcl) = oldbeadsx[tslice];
            YC(path.beads[tslice], ptcl) = oldbeadsy[tslice];
            ZC(path.beads[tslice], ptcl) = oldbeadsz[tslice];
        }

        return_defer(0); 
    }

    double newAction = 0.0;
    for (size_t tslice = 0; tslice < path.numTimeSlices; ++tslice) {
        newAction = newAction + PotentialAction(path, tslice);
    }
    
    // accept the move, or reject and restore the bead positions
    double u = mt_drand();
    double alpha = exp(-(newAction - oldAction));

    if (u < alpha) {
        return_defer(1);
    } else {
        for (size_t tslice = 0; tslice < path.numTimeSlices; ++tslice) {
            XC(path.beads[tslice], ptcl) = oldbeadsx[tslice];
            YC(path.beads[tslice], ptcl) = oldbeadsy[tslice];
            ZC(path.beads[tslice], ptcl) = oldbeadsz[tslice];
        }

        return_defer(0);
    }

defer:
    arena_reset(&arena);
    return result;
}

int Staging_Move(Path path, size_t ptcl)
// http://link.aps.org/doi/10.1103/PhysRevB.31.4234
{
    // NOTE: how should we set the stage length?  
    size_t stage_len = path.numTimeSlices/2; 
    assert(stage_len < path.numTimeSlices);

    int result = 0;

    size_t alpha_start = mt_lrand() % path.numTimeSlices;
    size_t alpha_end = (alpha_start + stage_len) % path.numTimeSlices; 

    double *oldbeadsx = (double*) arena_alloc(&arena, path.numTimeSlices * sizeof(double));
    memset(oldbeadsx, 0.0, path.numTimeSlices * sizeof(double));
    double *oldbeadsy = (double*) arena_alloc(&arena, path.numTimeSlices * sizeof(double));
    memset(oldbeadsy, 0.0, path.numTimeSlices * sizeof(double));
    double *oldbeadsz = (double*) arena_alloc(&arena, path.numTimeSlices * sizeof(double));
    memset(oldbeadsz, 0.0, path.numTimeSlices * sizeof(double));
    
    double oldAction = 0.0;

    for (size_t i = 1; i < stage_len; ++i) {
        size_t tslice = (alpha_start + i) % path.numTimeSlices;
        oldbeadsx[i - 1] = XC(path.beads[tslice], ptcl);
        oldbeadsy[i - 1] = YC(path.beads[tslice], ptcl);
        oldbeadsz[i - 1] = ZC(path.beads[tslice], ptcl);
        oldAction = oldAction + PotentialAction(path, tslice);
    }

    double newAction = 0.0;
    for (size_t i = 1; i < stage_len; ++i) {
        size_t tslice = (alpha_start + i) % path.numTimeSlices;
        size_t tslicem1 = (tslice - 1) % path.numTimeSlices;

        double tau1 = (stage_len - i) * path.tau;

        double avx = (tau1*XC(path.beads[tslicem1], ptcl) + path.tau*XC(path.beads[alpha_end], ptcl))/(path.tau+tau1);
        double avy = (tau1*YC(path.beads[tslicem1], ptcl) + path.tau*YC(path.beads[alpha_end], ptcl))/(path.tau+tau1);
        double avz = (tau1*ZC(path.beads[tslicem1], ptcl) + path.tau*ZC(path.beads[alpha_end], ptcl))/(path.tau+tau1);

        double sigma2 = 2.0*lam / (1.0/path.tau + 1.0/tau1);

        XC(path.beads[tslice], ptcl) = avx + sqrt(sigma2)*generate_normal(1.0);
        YC(path.beads[tslice], ptcl) = avy + sqrt(sigma2)*generate_normal(1.0);
        ZC(path.beads[tslice], ptcl) = avz + sqrt(sigma2)*generate_normal(1.0);
        newAction = newAction + PotentialAction(path, tslice); 
    }

    double u = mt_drand();
    double alpha = exp(-(newAction - oldAction));
    
    if (u < alpha) {
        return_defer(1);
    } else {
        for (size_t i = 1; i < stage_len; ++i) {
            size_t tslice = (alpha_start + i) % path.numTimeSlices;
            XC(path.beads[tslice], ptcl) = oldbeadsx[i - 1];
            YC(path.beads[tslice], ptcl) = oldbeadsy[i - 1];
            ZC(path.beads[tslice], ptcl) = oldbeadsz[i - 1];
        }

        return_defer(0);
    }

defer:
    arena_reset(&arena);
    return result;

}

void run_PIMC(Path path, size_t numSteps, bool collect_positions)
{
    AcceptanceRate acc = {0};

    size_t equilSkip = 20000;
    size_t observableSkip = 400;
    printf("Total MC steps: %zu\n", numSteps);
    printf("Equilibration skip: %zu\n", equilSkip);
    printf("Collecting 1 out of %zu steps\n", observableSkip); 

    //EnergyTrace trace = {0};

    for (size_t step = 0; step < numSteps; ++step) 
    {
        for (size_t i = 0; i < path.numParticles; ++i) {
            size_t ptcl = mt_lrand() % path.numParticles;
            acc.CenterOfMass += COM_Move(path, ptcl);
        }

        for (size_t i = 0; i < path.numParticles; ++i) {
            size_t ptcl = mt_lrand() % path.numParticles;
            acc.Staging += Staging_Move(path, ptcl);
        }

       if ((step % observableSkip == 0) && (step > equilSkip)) {
            //da_append(&trace, Energy(path)); 
            if (collect_positions) {
                assert(path.numParticles == 1);
                size_t ptcl = 0;

                for (size_t tslice = 0; tslice < path.numTimeSlices; ++tslice) {
                    double r = XC(path.beads[tslice], ptcl)*XC(path.beads[tslice], ptcl) + 
                               YC(path.beads[tslice], ptcl)*YC(path.beads[tslice], ptcl) + 
                               ZC(path.beads[tslice], ptcl)*ZC(path.beads[tslice], ptcl); 
                    r = sqrt(r);
                    
                    if (r < RMAX_COLLECT) {
                        da_append(&trace.positions, r);

                        double Uval = V_HeAr(r);
                        double svc_est = exp(-path.beta*Uval) - 1.0;
                        da_append(&trace.svcs, svc_est);
                    }
                }
            }
       } 
    }

    printf("--------------------------------------\n");
    printf("Acceptance Ratios:\n");
    printf("Center of Mass : %.3lf\n", (double) acc.CenterOfMass/numSteps/path.numParticles);
    printf("Staging        : %.3lf\n", (double) acc.Staging/numSteps/path.numParticles);
    printf("--------------------------------------\n");
}

void load_resources()
{
    int fileSize = 0;
    unsigned char* fileData = LoadFileData("resources/Alegreya-Regular.ttf", &fileSize);

    font.baseSize = FONT_SIZE_LOAD;
    font.glyphCount = 95;
    font.glyphs = LoadFontData(fileData, fileSize, FONT_SIZE_LOAD, 0, 95, FONT_SDF);
    Image atlas = GenImageFontAtlas(font.glyphs, &font.recs, 95, FONT_SIZE_LOAD, 4, 0);
    font.texture = LoadTextureFromImage(atlas);
    UnloadImage(atlas);
    UnloadFileData(fileData);
}

void draw_histogram(const char* title, gsl_histogram *h) 
{
    size_t nbins = h->n;
    double ymax = gsl_histogram_max_val(h);
    double xmin = h->range[0];
    double xmax = h->range[nbins];
   
    size_t samples_count = 0;
    for (size_t i = 0; i < nbins; ++i) {
        samples_count += (int) h->bin[i];
    }
    
    double mean = gsl_histogram_mean(h);

    SetConfigFlags(FLAG_WINDOW_RESIZABLE);
    SetConfigFlags(FLAG_MSAA_4X_HINT);
    
    size_t factor = 80;
    InitWindow(16*factor, 9*factor, title);
    SetExitKey(KEY_Q);
   
    load_resources();

    SetTargetFPS(60);

    while (!WindowShouldClose()) {
        BeginDrawing();
        ClearBackground(COLOR_BACKGROUND);
        
        int screen_width = GetScreenWidth();
        int screen_height = GetScreenHeight();

        int simul_sz = (int) (0.8 * fminf(screen_width, screen_height));

        Rectangle world = {
            .x = 0.1*screen_width, 
            .y = 0.1*screen_height, 
            .width = simul_sz, 
            .height = simul_sz
        };
            
        double rect_width = world.width / nbins;

        for (size_t i = 0; i < nbins; ++i) {
            double height = gsl_histogram_get(h, i)/ymax * world.height; 
            double factor = 0.9;

            Rectangle r = {
                .x = world.x + i*rect_width,
                .y = world.y + world.height - factor*height, 
                .width = rect_width,
                .height = factor*height, 
            };

            DrawRectangleLinesEx(r, 2.0, YELLOW);
        }

        DrawRectangleLinesEx(world, 3.0, LIGHTGRAY);

        
        int font_size = 24;
        {
            Vector2 mean_st = {
                .x = world.x + (mean - xmin)/(xmax - xmin)*world.width, 
                .y = world.y,
            };
            Vector2 mean_end = {
                .x = world.x + (mean - xmin)/(xmax - xmin)*world.width, 
                .y = world.y + world.height,
            };
            DrawLineEx(mean_st, mean_end, 2.0, RED);

            const char *buffer = TextFormat("%.3lf", mean);
            Vector2 text_len = MeasureTextEx(font, buffer, font_size, 0);
            Vector2 text_pos = {
                world.x + (mean - xmin)/(xmax - xmin)*world.width - 0.25*text_len.x,
                world.y + world.height + 0.25 * text_len.y,
            };
            DrawTextEx(font, buffer, text_pos, font_size, 0, RED);
        }

        {
            const char *buffer = TextFormat("%.3lf", xmin);
            Vector2 text_len = MeasureTextEx(font, buffer, font_size, 0);
            Vector2 text_pos = {
                world.x - 0.25 * text_len.x,
                world.y + world.height + 0.25 * text_len.y,
            };
            DrawTextEx(font, buffer, text_pos, font_size, 0, WHITE);
        }
        {
            const char *buffer = TextFormat("%.3lf", xmax);
            Vector2 text_len = MeasureTextEx(font, buffer, font_size, 0);
            Vector2 text_pos = {
                world.x + world.width - 0.25 * text_len.x,
                world.y + world.height + 0.25 * text_len.y,
            };
            DrawTextEx(font, buffer, text_pos, font_size, 0, WHITE);
        }
        {
            const char *buffer = TextFormat("Samples: %zu", samples_count);
            Vector2 text_pos = {
                world.x + world.width + 50,
                world.y,
            };
            DrawTextEx(font, buffer, text_pos, font_size, 0, WHITE); 
        }

        EndDrawing();
    }

    CloseWindow();

}

double minval(double *arr, size_t count)
{
    double c = FLT_MAX;
    for (size_t i = 0; i < count; ++i) {
        if (arr[i] < c) c = arr[i];
    }

    return c;
}

double maxval(double *arr, size_t count) 
{
    double c = FLT_MIN;
    for (size_t i = 0; i < count; ++i) {
        if (arr[i] > c) c = arr[i];
    }

    return c;
}

int main()
{
    uint32_t seed = mt_goodseed();
    mt_seed32(seed);

    double T = 300.0; // K
    double beta = 1.0/(Boltzmann_Hartree * T); 

    path.numParticles = 1;
    path.numTimeSlices = 4;
    path.tau = beta/path.numTimeSlices;
    path.beta = beta;
    alloc_beads(&path);

    printf("Simulation parameters:\n");
    printf("Number of Particles   = %zu\n", path.numParticles);
    printf("Number of Time Slices = %zu\n", path.numTimeSlices);
    printf("beta                  = %.3lf\n", beta);
    printf("tau                   = %.3lf\n", path.tau);

    for (size_t i = 0; i < path.numTimeSlices; ++i) {
        for (size_t j = 0; j < 3*path.numParticles; ++j) {
            path.beads[i][j] = RMAX_CUBE * 0.5*(-1.0 + 2.0*mt_drand()); 
        }
    }
   
    size_t MC_steps = 30 * 1000 * 1000;
    run_PIMC(path, MC_steps, true);

    double Lambda = /* Planck */ 2.0*M_PI / sqrt(2.0*M_PI * mu / beta); // a.u.^3
    printf("Lambda: %.5f\n", Lambda);

    double V = RMAX_COLLECT*RMAX_COLLECT*RMAX_COLLECT;
    double Lambda3 = Lambda*Lambda*Lambda; 
    double F = Lambda3 * pow(pow(path.numTimeSlices, 3.0/2.0) / Lambda3, path.numTimeSlices);
    printf("F = %.5f\n", F);

    double Coeff = -0.5 * V * ALU*ALU*ALU * AVOGADRO * 1e6; 
    for (size_t i = 0; i < trace.svcs.count; ++i) {
        trace.svcs.items[i] *=  F * Coeff;
    }

    printf("Collected SVC vals: %zu\n", trace.svcs.count);
    
    double svc_min = minval(trace.svcs.items, trace.svcs.count);
    double svc_max = maxval(trace.svcs.items, trace.svcs.count);
    printf("Minimum value: %lf\n", svc_min);
    printf("Maximum value: %lf\n", svc_max);

    gsl_histogram *p_histogram;

    size_t nbins = 100;
    p_histogram = gsl_histogram_alloc(nbins);
    gsl_histogram_set_ranges_uniform(p_histogram, -50.0, 50.0);


    for (size_t i = 0; i < trace.positions.count; ++i) {
        gsl_histogram_increment(p_histogram, trace.svcs.items[i]);
    }

    gsl_histogram_fprintf(stdout, p_histogram, "%g", "%g");

    draw_histogram("SVC histogram", p_histogram); 

    gsl_histogram_free(p_histogram);
    return 0;
}

// T = 90 K => svc = -0.5 (1.18)
// T = 300 K -> svc = -0.3 (18.5)
