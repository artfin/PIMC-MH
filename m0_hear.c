#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <float.h>

#include <mpi.h>

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


#define HEAR_IMPLEMENTATION
#include "HeAr.h"

#include "raylib.h"

#define COLOR_BACKGROUND GetColor(0x181818FF)

#define FONT_SIZE_LOAD 160 
Font font = {0};

#include "common.h"

#define PROTOCOL_IMPLEMENTATION 
#include "protocol.h"

// ---------------------------------------------------------
#define m_He (4.00260325413 * RAMTOAMU)
#define m_Ar (39.9623831237 * RAMTOAMU)
#define mu m_He * m_Ar / (m_He + m_Ar)  // a.u.
                                         
#define lam 1.0/(2.0*mu) // hbar^2/2m

#define COORD_SAMPLE_MAX  4.0 // a.u -- we sample coordinates within this cube ???
#define RMIN_COLLECT 4.0  // a.u.
#define RMAX_COLLECT 30.0 // a.u.
#define COORD_MAX    50.0 // a.u.
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
} M0_Trace;

typedef struct {
    double *items;
    size_t count;
    size_t capacity;
} PositionTrace;

typedef struct {
    PositionTrace positions;
    M0_Trace m0s;
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

void sample_beads(Path *path)
{
    for (size_t i = 0; i < path->numTimeSlices; ++i) {
        for (size_t j = 0; j < 3*path->numParticles; ++j) {
            path->beads[i][j] = COORD_SAMPLE_MAX * 0.5*(-1.0 + 2.0*mt_drand()); 
        }
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

    double delta = 1.0e-2; 
    double shiftx = delta*COORD_SAMPLE_MAX*(-1.0 + 2.0*mt_drand());
    double shifty = delta*COORD_SAMPLE_MAX*(-1.0 + 2.0*mt_drand());
    double shiftz = delta*COORD_SAMPLE_MAX*(-1.0 + 2.0*mt_drand());

    double oldAction = 0.0;
    for (size_t tslice = 0; tslice < path.numTimeSlices; ++tslice) {
        oldAction = oldAction + PotentialAction(path, tslice);
    } 
    
    double *oldbeadsx = (double*) arena_alloc(&arena, path.numTimeSlices * sizeof(double));
    double *oldbeadsy = (double*) arena_alloc(&arena, path.numTimeSlices * sizeof(double));
    double *oldbeadsz = (double*) arena_alloc(&arena, path.numTimeSlices * sizeof(double));
    memset(oldbeadsx, 0.0, path.numTimeSlices * sizeof(double));
    memset(oldbeadsy, 0.0, path.numTimeSlices * sizeof(double));
    memset(oldbeadsz, 0.0, path.numTimeSlices * sizeof(double));

    for (size_t tslice = 0; tslice < path.numTimeSlices; ++tslice) {
        oldbeadsx[tslice] = XC(path.beads[tslice], ptcl);
        oldbeadsy[tslice] = YC(path.beads[tslice], ptcl);
        oldbeadsz[tslice] = ZC(path.beads[tslice], ptcl);
        XC(path.beads[tslice], ptcl) = XC(path.beads[tslice], ptcl) + shiftx;
        YC(path.beads[tslice], ptcl) = YC(path.beads[tslice], ptcl) + shifty;
        ZC(path.beads[tslice], ptcl) = ZC(path.beads[tslice], ptcl) + shiftz;
    }
    
    // reject move that moves COM in such that any particles goes outside of the cube
    bool within_cube = true;

    for (size_t tslice = 0; tslice < path.numTimeSlices; ++tslice) {
        if (abs(XC(path.beads[tslice], ptcl) > COORD_MAX) ||
            abs(YC(path.beads[tslice], ptcl) > COORD_MAX) ||
            abs(ZC(path.beads[tslice], ptcl) > COORD_MAX)) 
        {
            within_cube = false;
            break;
        }
    }

    if (!within_cube) {
        printf("resampling the chain\n");
        sample_beads(&path);
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

        double attx = avx + sqrt(sigma2)*generate_normal(1.0); 
        double atty = avy + sqrt(sigma2)*generate_normal(1.0);
        double attz = avz + sqrt(sigma2)*generate_normal(1.0);

        bool allow = (fabs(attx) < COORD_MAX) && (fabs(atty) < COORD_MAX) && (fabs(attz) < COORD_MAX); 

        if (allow) {
            XC(path.beads[tslice], ptcl) = attx;
            YC(path.beads[tslice], ptcl) = atty;
            ZC(path.beads[tslice], ptcl) = attz; 
        
            newAction = newAction + PotentialAction(path, tslice); 
        }
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
   
    return_defer(1);

defer:
    arena_reset(&arena);
    return result;

}

void pimc_driver(MPI_Context ctx, Path path, size_t numSteps, int sockfd, bool collect_positions)
{
    AcceptanceRate acc = {0};

    size_t equilSkip = 20000;
    size_t observableSkip = 400;
    printf("Total MC steps: %zu\n", numSteps);
    printf("Equilibration skip: %zu\n", equilSkip);
    printf("Collecting 1 out of %zu steps\n", observableSkip); 

    size_t send_size = 1000;
    assert(send_size <= 1000 && " NOTE: keep the send_size under 100 for now. In the local network we encounter problems with sending larger packets\n");

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
            if (collect_positions) {
                assert(path.numParticles == 1);

                double m0_est = 0.0;
                size_t c = 0;

                for (size_t tslice = 0; tslice < path.numTimeSlices; ++tslice) {
                    double r = XC(path.beads[tslice], 0)*XC(path.beads[tslice], 0) + 
                               YC(path.beads[tslice], 0)*YC(path.beads[tslice], 0) + 
                               ZC(path.beads[tslice], 0)*ZC(path.beads[tslice], 0); 
                    r = sqrt(r);
                    
                    if ((r > RMIN_COLLECT) && (r < RMAX_COLLECT)) {
                        da_append(&trace.positions, r);

                        double dipval = dip_HeAr(r);
                        m0_est += dipval*dipval; 
                        c = c + 1;
                    }
                }
                
                if (c > 0) {
                    m0_est = m0_est / c;
                    da_append(&trace.m0s, m0_est);
                }

                if (trace.positions.count == send_size) {
                    if (ctx.rank > 0) {
                        MPI_Send(trace.positions.items, send_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
                        trace.positions.count = 0;
                    } else {
                        printf("Sending packets...\n");
                        sendFloat64Array(sockfd, trace.positions.items, trace.positions.count); 
                        trace.positions.count = 0;
                        
                        MPI_Status status = {0}; 
                        for (int i = 1; i < ctx.size; ++i) {
                            MPI_Recv(trace.positions.items, send_size, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                            trace.positions.count = send_size;
                            sendFloat64Array(sockfd, trace.positions.items, trace.positions.count);
                            trace.positions.count = 0;
                        }
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

/*
void draw_histogram(const char* title, gsl_histogram *h, bool log_scale) 
{
    size_t nbins = h->n;
    double ymax = gsl_histogram_max_val(h);
    if (log_scale) ymax = log(ymax);

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
            double height; 
            if (log_scale) { 
                height = log(gsl_histogram_get(h, i))/ymax * world.height;
            } else {
                height = gsl_histogram_get(h, i)/ymax * world.height;
            }
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

            const char *buffer = TextFormat("%.3e", mean);
            Vector2 text_len = MeasureTextEx(font, buffer, font_size, 0);
            Vector2 text_pos = {
                world.x + (mean - xmin)/(xmax - xmin)*world.width - 0.25*text_len.x,
                world.y + world.height + 0.75 * text_len.y,
            };
            DrawTextEx(font, buffer, text_pos, font_size, 0, RED);
        }

        {
            const char *buffer = TextFormat("%.3e", xmin);
            Vector2 text_len = MeasureTextEx(font, buffer, font_size, 0);
            Vector2 text_pos = {
                world.x - 0.25 * text_len.x,
                world.y + world.height + 0.25 * text_len.y,
            };
            DrawTextEx(font, buffer, text_pos, font_size, 0, WHITE);
        }
        {
            const char *buffer = TextFormat("%.3e", xmax);
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
*/

typedef struct {
    double min;
    double max;
    double mean;
    double std; 
} Stats;  

Stats getStats(double *arr, size_t count) 
{
    Stats s = {0};
    
    s.max = FLT_MIN;
    for (size_t i = 0; i < count; ++i) {
        if (arr[i] > s.max) s.max = arr[i];
    }

    s.min = FLT_MAX;
    for (size_t i = 0; i < count; ++i) {
        if (arr[i] < s.min) s.min = arr[i];
    }
    
    for (size_t i = 0; i < count; ++i) {
        s.mean = s.mean + arr[i];
    }

    s.mean /= count;

    double r = 0;
    for (size_t i = 0; i < count; ++i) {
        r = r + (arr[i] - s.mean) * (arr[i] - s.mean); 
    }

    r = r / (count - 1);
    s.std = sqrt(r);

    return s;
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    
    MPI_Context ctx = {0}; 
    MPI_Comm_size(MPI_COMM_WORLD, &ctx.size);
    MPI_Comm_rank(MPI_COMM_WORLD, &ctx.rank);

    uint32_t seed = mt_goodseed();
    mt_seed32(seed);
        
    double T = 300.0; // K
    double beta = 1.0/(Boltzmann_Hartree * T); 

    path.numParticles = 1;
    path.numTimeSlices = 1;
    path.tau = beta/path.numTimeSlices;
    path.beta = beta;
    alloc_beads(&path);

    // Source: M0 at 295K from diploma
    double refVal = 5.27e-05; 
    int blockSize = 1; // TODO: ignore for now

    printf("Simulation parameters:\n");
    printf("Number of Particles   = %zu\n", path.numParticles);
    printf("Number of Time Slices = %zu\n", path.numTimeSlices);
    printf("beta                  = %.3lf\n", beta);
    printf("tau                   = %.3lf\n", path.tau);

    sample_beads(&path);

    int sockfd = 0; 
    if (ctx.rank == 0) { 
        sockfd = initClient();
        printf("connection established at socket = %d\n", sockfd);

        if (sockfd > 0) {
            sendFloat64(sockfd, path.beta);
            sendInt32(sockfd, (int) path.numTimeSlices);
            sendInt32(sockfd, ctx.size);
            sendFloat64(sockfd, refVal);
            recvInt32(sockfd, &blockSize);
        } else {
            fprintf(stderr, "ERROR: client could not connect to server\n");
            fprintf(stderr, "Continuing calculation without communicating with the server\n\n");
        }
    }
   
    size_t MC_steps =  1 * 1000 * 1000 * 1000;
    pimc_driver(ctx, path, MC_steps, sockfd, true);

    double Lambda = /* Planck */ 2.0*M_PI / sqrt(2.0*M_PI * mu / beta); // a.u.^3
    printf("Lambda: %.5f\n", Lambda);

    double V = 4.0*M_PI/3.0*RMAX_COLLECT*RMAX_COLLECT*RMAX_COLLECT;
    const double ZeroCoeff = 0.00361479637 / (4.0 * M_PI);

    for (size_t i = 0; i < trace.m0s.count; ++i) {
        trace.m0s.items[i] *= V * ZeroCoeff;
    }

    printf("Collected M0 estimates: %zu\n", trace.m0s.count);
    
    Stats m0 = getStats(trace.m0s.items, trace.m0s.count);
    printf("Minimum: %.10lf\n", m0.min);
    printf("Maximum: %.10lf\n", m0.max);
    printf("Mean: %.3e\n", m0.mean); 

    /*
    gsl_histogram *p_histogram;

    size_t nbins = 100;
    p_histogram = gsl_histogram_alloc(nbins);
    //gsl_histogram_set_ranges_uniform(p_histogram, 0.0, 30.0);
    gsl_histogram_set_ranges_uniform(p_histogram, 0.0, 1.0e-3);

    for (size_t i = 0; i < trace.m0s.count; ++i) {
        gsl_histogram_increment(p_histogram, trace.m0s.items[i]);
    }

    gsl_histogram_fprintf(stdout, p_histogram, "%g", "%g");

    draw_histogram("M0 histogram", p_histogram, true); 

    gsl_histogram_free(p_histogram);
    */

    MPI_Finalize();

    return 0;
}

// M0 (true) = 5.27e-6
// T = 300 K -> M0 (est) = 5.47e-6, N = 150 mln
