// Based on the notes by Adrian Del Maestro 
// "Path Integral Monte Carlo and the Worm algorithm in the Spatial Continuum"

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

#include "raylib.h"

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

void subcmd_visualize(const char *program_path, int argc, char **argv);
void subcmd_run(const char *program_path, int argc, char **argv);

Subcmd subcmds[] = {
    DEFINE_SUBCMD(visualize, "Visualize the chain modifications"),
    DEFINE_SUBCMD(run, "Run the PIMC calculation of the harmonic oscillator"),
};
#define SUBCMDS_COUNT (sizeof(subcmds)/sizeof(subcmds[0]))

typedef struct {
    double *items;
    size_t count;
    size_t capacity;
} EnergyTrace;

// ---------------------------------------------------------
double lam = 0.5; // hbar^2/2m k_B
// ---------------------------------------------------------
    
Path path = {0};

double V(double x) { return 0.5*x*x; }

double SHOExact(double beta) {
    //return 0.5/tanh(0.5/T);
    return 0.5/tanh(0.5*beta);
}

double PotentialAction(Path path, size_t tslice) 
{
    assert(tslice < path.numTimeSlices);

    double pot = 0.0;
    for (size_t ptcl = 0; ptcl < path.numParticles; ++ptcl) {
        pot = pot + V(path.beads[tslice][ptcl]);
    }

    return path.tau * pot;
}

double PotentialEnergy(Path path) {
    double result = 0.0;

    for (size_t tslice = 0; tslice < path.numTimeSlices; ++tslice) {
        for (size_t ptcl = 0; ptcl < path.numParticles; ++ptcl) {
            double x = path.beads[tslice][ptcl];
            result = result + V(x);
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
            double dx = path.beads[tslicep1][ptcl] - path.beads[tslice][ptcl];
            tot = tot - norm*dx*dx;
        }
    }

    return 0.5*path.numParticles/path.tau + tot/path.numTimeSlices;
}

double Energy(Path path) {
    return PotentialEnergy(path) + KineticEnergy(path);
}

typedef struct {
    double mean;
    double std;
} Stats;

Stats getStats(EnergyTrace t) 
{
    Stats s = {0};

    for (size_t i = 0; i < t.count; ++i) {
        s.mean = s.mean + t.items[i];
    }

    s.mean /= t.count;

    double r = 0;
    for (size_t i = 0; i < t.count; ++i) {
        r = r + (t.items[i] - s.mean) * (t.items[i] - s.mean); 
    }

    r = r / (t.count - 1);
    s.std = sqrt(r);

    return s;
}

int compar(const void *a, const void *b) {
    double *x = (double *) a;
    double *y = (double *) b;
    if (*x < *y) { 
        return -1;
    } else if (*x > *y) {
        return 1; 
    }
    return 0;
}

Stats getStatsEx(EnergyTrace t, size_t binSize)
// TODO: make sure that binning actually improves the std 
{
    qsort(t.items, t.count, sizeof(t.items[0]), compar);    

    int numBins = ceil((float) t.count / binSize);
    printf("numBins: %d\n", numBins);

    double *bins = (double*) arena_alloc(&arena, numBins * sizeof(double));
    memset(bins, 0.0, numBins * sizeof(double)); 

    int c = 0; // bin cursor
    for (int i = 0; i < t.count; ++i) 
    {
        if ((i > 0) && (i % binSize == 0)) {
            bins[c] /= binSize;
            c++;
        }

        bins[c] = bins[c] + t.items[i];
    }

    if (t.count % binSize != 0) {
        bins[c] /= (t.count % binSize);
    }

    //printf("Bins: ");
    //for (int i = 0; i < numBins; ++i) {
    //    printf("%.2lf ", bins[i]);
    //}
    //printf("\n");

    Stats s = {0};

    for (size_t i = 0; i < numBins; ++i) {
        s.mean = s.mean + bins[i];
    }

    s.mean /= numBins;

    double r = 0;
    for (size_t i = 0; i < numBins; ++i) {
        r = r + (bins[i] - s.mean) * (bins[i] - s.mean); 
    }

    r = r / (numBins - 1);
    s.std = sqrt(r);

    return s; 
}


int COM_Move(Path path, size_t ptcl)
{
    assert(ptcl < path.numParticles);
    int result = 0;

    double delta = 0.75;
    double shift = delta*(-1.0 + 2.0*mt_drand());

    double oldAction = 0.0;
    for (size_t tslice = 0; tslice < path.numTimeSlices; ++tslice) {
        oldAction = oldAction + PotentialAction(path, tslice);
    } 
    
    double *oldbeads = (double*) arena_alloc(&arena, path.numTimeSlices * sizeof(double));
    memset(oldbeads, 0.0, path.numTimeSlices * sizeof(double));

    for (size_t tslice = 0; tslice < path.numTimeSlices; ++tslice) {
        oldbeads[tslice] = path.beads[tslice][ptcl];
        path.beads[tslice][ptcl] = path.beads[tslice][ptcl] + shift;
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
            path.beads[tslice][ptcl] = oldbeads[tslice];
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
    size_t stage_len = path.numTimeSlices/2; 
    assert(stage_len < path.numTimeSlices);

    int result = 0;

    size_t alpha_start = mt_lrand() % path.numTimeSlices;
    size_t alpha_end = (alpha_start + stage_len) % path.numTimeSlices; 

    double *oldbeads = (double*) arena_alloc(&arena, (stage_len - 1)*sizeof(double));
    memset(oldbeads, 0.0, (stage_len - 1)*sizeof(double));
    
    double oldAction = 0.0;

    for (size_t i = 1; i < stage_len; ++i) {
        size_t tslice = (alpha_start + i) % path.numTimeSlices;
        oldbeads[i - 1] = path.beads[tslice][ptcl];
        oldAction = oldAction + PotentialAction(path, tslice);
    }

    double newAction = 0.0;
    for (size_t i = 1; i < stage_len; ++i) {
        size_t tslice = (alpha_start + i) % path.numTimeSlices;
        size_t tslicem1 = (tslice - 1) % path.numTimeSlices;

        double tau1 = (stage_len - i) * path.tau;
        double avex = (tau1*path.beads[tslicem1][ptcl] + path.tau*path.beads[alpha_end][ptcl])/(path.tau+tau1);
        double sigma2 = 2.0*lam / (1.0/path.tau + 1.0/tau1);

        path.beads[tslice][ptcl] = avex + sqrt(sigma2)*generate_normal(1.0);
        //printf("avex: %.5lf; sigma2: %.5lf => pos = %.5lf\n", avex, sigma2, path.beads[tslice][ptcl]);
        newAction = newAction + PotentialAction(path, tslice); 
    }

    double u = mt_drand();
    double alpha = exp(-(newAction - oldAction));
    
    if (u < alpha) {
        return_defer(1);
    } else {
        for (size_t i = 1; i < stage_len; ++i) {
            size_t tslice = (alpha_start + i) % path.numTimeSlices;
            path.beads[tslice][ptcl] = oldbeads[i - 1];
        }

        return_defer(0);
    }

defer:
    arena_reset(&arena);
    return result;

}

EnergyTrace run_PIMC(Path path, size_t numSteps)
{
    AcceptanceRate acc = {0};

    int equilSkip = 20000;
    int observableSkip = 400;
    printf("Total MC steps: %zu\n", numSteps);
    printf("Equilibration skip: %d\n", equilSkip);
    printf("Collecting 1 out of %d steps\n", observableSkip); 

    EnergyTrace trace = {0};

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
            da_append(&trace, Energy(path)); 
       } 
    }

    printf("--------------------------------------\n");
    printf("Acceptance Ratios:\n");
    printf("Center of Mass : %.3lf\n", (double) acc.CenterOfMass/numSteps/path.numParticles);
    printf("Staging        : %.3lf\n", (double) acc.Staging/numSteps/path.numParticles);
    printf("--------------------------------------\n");

    return trace;
}

void alloc_beads(Path *path)
{
    path->beads = (double**) malloc(path->numTimeSlices * sizeof(double*));
    for (size_t i = 0; i < path->numTimeSlices; ++i) {
        path->beads[i] = (double*) malloc(path->numParticles * sizeof(double));
        memset(path->beads[i], 0.0, path->numParticles*sizeof(double));
    }
}



#define COLOR_BACKGROUND GetColor(0x181818FF)

void update_draw_frame()
{
    static int simulation_step = 0; 

    if (IsKeyPressed(KEY_SPACE)) {
        printf("Simulation step: %d\n", simulation_step);

        if (simulation_step % 2 == 0) {
            for (size_t ptcl = 0; ptcl < path.numParticles; ++ptcl) {
                COM_Move(path, ptcl);
            }
        } else {
            for (size_t ptcl = 0; ptcl < path.numParticles; ++ptcl) {
                Staging_Move(path, ptcl);
            }
        }

        simulation_step++;
    }

    BeginDrawing();

    ClearBackground(COLOR_BACKGROUND);
   
    Color colors[] = {
        GetColor(0xF2AF29FF),
        GetColor(0xB52A2AFF),
        GetColor(0x7DD181FF),
    };
    #define COLORS_COUNT (sizeof(colors)/sizeof(colors[0]))

    int screen_width = GetScreenWidth();
    int screen_height = GetScreenHeight();

    int simul_sz = (int) (0.8 * fminf(screen_width, screen_height));

    Rectangle world = {
        .x = 0.1*screen_width, 
        .y = 0.1*screen_height, 
        .width = simul_sz, 
        .height = simul_sz
    };

    DrawRectangleLinesEx(world, 3.0, LIGHTGRAY);
    double GRID_MIN = -2.5;
    double GRID_MAX =  2.5;   

    double tau = world.height / (path.numTimeSlices - 1);

    double prevx, prevy;
    for (size_t ptcl = 0; ptcl < path.numParticles; ++ptcl) {
        assert(ptcl < COLORS_COUNT);
        Color c = colors[ptcl];

        for (size_t i = 0; i < path.numTimeSlices; ++i) {
            double bead = path.beads[i][ptcl];
            if ((bead < GRID_MIN) || (bead > GRID_MAX)) {
                printf("%.5lf\n", bead);
            }

            double screenx = world.x + (bead - GRID_MIN) / (GRID_MAX - GRID_MIN) * simul_sz;
            double screeny = world.y + tau * i; 
            DrawCircle(screenx, screeny, 6.0, c);

            if (i > 0) {
                DrawLine(screenx, screeny, prevx, prevy, c);
            }

            prevx = screenx;
            prevy = screeny;
        }
    }

    EndDrawing(); 
}

void subcmd_visualize(const char *program_path, int argc, char **argv)
{
    SetConfigFlags(FLAG_WINDOW_RESIZABLE);
    SetConfigFlags(FLAG_MSAA_4X_HINT);
    
    size_t factor = 80;
    InitWindow(16*factor, 9*factor, "PIMC");
    SetExitKey(KEY_Q);

    double T = 1.0; // K
    path.numParticles = 1;
    path.numTimeSlices = 16;
    path.tau = 1.0/(T * path.numTimeSlices);
    alloc_beads(&path);
    
    for (size_t i = 0; i < path.numTimeSlices; ++i) {
        for (size_t j = 0; j < path.numParticles; ++j) {
            path.beads[i][j] = 0.5 * (-1.0 + 2.0*mt_drand()); 
        }
    }

    SetTargetFPS(60);
    while (!WindowShouldClose()) {
        update_draw_frame();
    }

    CloseWindow();
}

void subcmd_run(const char *program_path, int argc, char **argv)
{
    uint32_t seed = mt_goodseed();
    mt_seed32(seed);

    //double T = 1.0; // K
    double beta = 10.0;

    path.numParticles = 1;
    path.numTimeSlices = 32;
    path.tau = beta/path.numTimeSlices;
    alloc_beads(&path);

    printf("Simulation parameters:\n");
    printf("Number of Particles   = %zu\n", path.numParticles);
    printf("Number of Time Slices = %zu\n", path.numTimeSlices);
    printf("beta                  = %.3lf\n", beta);
    printf("tau                   = %.3lf\n", path.tau);

    for (size_t i = 0; i < path.numTimeSlices; ++i) {
        for (size_t j = 0; j < path.numParticles; ++j) {
            path.beads[i][j] = 0.5 * (-1.0 + 2.0*mt_drand()); 
        }
    }
    
    size_t MC_steps = 20 * 1000 * 1000;
    EnergyTrace t = run_PIMC(path, MC_steps);

    int binSize = 500;
    Stats s = getStatsEx(t, binSize);
    printf("Collected %zu values\n", t.count); 
    printf("(PIMC) Energy = %.5f +/- %.5f\n", s.mean, s.std);

    double en_exact = SHOExact(beta);
    printf("(Exact) Energy = %.5f\n", en_exact);

    double err = fabsf(s.mean - en_exact) / en_exact;
    printf("Error: %.2f%%\n", err*100.0);
}


int main(int argc, char *argv[])
{
    char *program_path = shift(&argc, &argv);

    if (argc <= 0) {
        usage(program_path, subcmds, SUBCMDS_COUNT);
        fprintf(stderr, "ERROR: no subcommand is provided\n");
        exit(1);
    }

    const char *subcmd_id = shift(&argc, &argv);
    Subcmd *subcmd = find_subcmd_by_id(subcmds, SUBCMDS_COUNT, subcmd_id);

    if (subcmd != NULL) {
        subcmd->run(program_path, argc, argv);
    } else {
        usage(program_path, subcmds, SUBCMDS_COUNT);
        fprintf(stderr, "ERROR: unknown subcommand  `%s`\n", subcmd_id);
        exit(1);
    }

    return 0;
}

int main2()
{
    EnergyTrace t = {0};
    da_append(&t, 6.0);
    da_append(&t, 1.0);
    da_append(&t, 2.0);
    da_append(&t, 3.0);
    da_append(&t, 4.0);

    Stats s = getStatsEx(t, 2);
    printf("m = %.5lf +/- %.5lf\n", s.mean, s.std);

    return 0;
}
