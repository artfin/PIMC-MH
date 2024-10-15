// Based on the notes by Adrian Del Maestro 
// "Path Integral Monte Carlo and the Worm algorithm in the Spatial Continuum"

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <unistd.h>
#include <inttypes.h>


#ifndef NO_MPI
#include <mpi.h>
#endif

// TODO: can the code for gsl_histogram be extracted and brought in the repo? 
#ifndef NO_VISUALIZE
#include <gsl/gsl_histogram.h>
#include "raylib.h"

#define FONT_SIZE_LOAD 160 
Font font = {0};

#define COLOR_BACKGROUND GetColor(0x181818FF)
#endif


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

#define MAX_FLOAT_TEXTLEN 100

#define ARENA_IMPLEMENTATION
#include "arena.h"
static Arena arena = {0};

#define COMMON_IMPLEMENTATION
#include "common.h"

void subcmd_visualize(MPI_Context ctx, int argc, char **argv);
void subcmd_run(MPI_Context ctx, int argc, char **argv);
void subcmd_server(MPI_Context ctx, int argc, char **argv);

MPI_Subcmd subcmds[] = {
#ifndef NO_VISUALIZE
    DEFINE_SUBCMD(visualize, "Visualize the chain modifications"),
#endif
    DEFINE_SUBCMD(run, "Run the PIMC calculation of the harmonic oscillator"),
};
#define SUBCMDS_COUNT (sizeof(subcmds)/sizeof(subcmds[0]))

typedef struct {
    double *items;
    size_t count;
    size_t capacity;
} EnergyTrace;

typedef struct {
    double *items;
    size_t count;
    size_t capacity;
} PositionTrace;

typedef struct {
    EnergyTrace energies;
    PositionTrace positions;
} PIMC_Trace;

PIMC_Trace trace = {0};

//#define USE_UNIX_SOCKET
#define PROTOCOL_IMPLEMENTATION 
#include "protocol.h"

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
    for (size_t i = 0; i < t.count; ++i) 
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

    for (int i = 0; i < numBins; ++i) {
        s.mean = s.mean + bins[i];
    }

    s.mean /= numBins;

    double r = 0;
    for (int i = 0; i < numBins; ++i) {
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

void pimc_driver(MPI_Context ctx, Path path, size_t numSteps, int sockfd, bool collect_positions)
{
    AcceptanceRate acc = {0};

    size_t equilSkip = 20000;
    size_t observableSkip = 400;

    size_t send_size = 100;
    assert(send_size <= 100 && " NOTE: keep the send_size under 100 for now. In the local network we encounter problems with sending larger packets\n");

    printf("Total MC steps: %zu\n", numSteps);
    printf("Equilibration skip: %zu\n", equilSkip);
    printf("Collecting 1 out of %zu steps\n", observableSkip); 

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
            da_append(&trace.energies, Energy(path));

            if (collect_positions) {
                assert(path.numParticles == 1);
                for (size_t tslice = 0; tslice < path.numTimeSlices; ++tslice) {
                    da_append(&trace.positions, path.beads[tslice][0]);
                }
            }

            // if (sockfd <= 0) continue;
            // TODO: if the socket is not open I still want to be able to calculate the mean value
            // we should calculate statistics on the fly instead of accumulating EVERY energy value 

            if (trace.energies.count == send_size) {
#ifndef NO_MPI
                if (ctx.rank > 0) {
                    MPI_Send(trace.energies.items, send_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
                    trace.energies.count = 0;
                } else 
#endif // NO_MPI
                {
                    sendFloat64Array(sockfd, trace.energies.items, trace.energies.count); 
                    trace.energies.count = 0;
                    printf("Sending energies...\n");

#ifndef NO_MPI
                    MPI_Status status = {0}; 
                    for (int i = 1; i < ctx.size; ++i) {
                        MPI_Recv(trace.energies.items, send_size, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                        trace.energies.count = send_size;
                        sendFloat64Array(sockfd, trace.energies.items, trace.energies.count);
                        trace.energies.count = 0;
                    }
#endif // NO_MPI
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

void alloc_beads(Path *path)
{
    path->beads = (double**) malloc(path->numTimeSlices * sizeof(double*));
    for (size_t i = 0; i < path->numTimeSlices; ++i) {
        path->beads[i] = (double*) malloc(path->numParticles * sizeof(double));
        memset(path->beads[i], 0.0, path->numParticles*sizeof(double));
    }
}

void dealloc_beads(Path *path)
{
    for (size_t i = 0; i < path->numTimeSlices; ++i) {
        free(path->beads[i]);
    }

    free(path->beads);
}

#ifndef NO_VISUALIZE
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

void subcmd_visualize(MPI_Context ctx, int argc, char **argv)
{
    (void) ctx;
    (void) argc;
    (void) argv;

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
#endif

void subcmd_run(MPI_Context ctx, int argc, char **argv)
{
    bool opt_client = false; 
    bool opt_position_histogram = false;

    if (argc > 0) {
        char *opt = shift(&argc, &argv);

        if (strcmp(opt, "--client") == 0) {
            fprintf(stdout ,"> start communicating to the server\n");
            opt_client = true;
        } else if (strcmp(opt, "--position_histogram") == 0) {
            fprintf(stdout, "-- Collecting positions in the histogram\n");
            opt_position_histogram = true; 
        } else {
            // TODO: show available options
            fprintf(stderr, "ERROR: unknown option passed to `run` subcommand: %s\n", opt);
            exit(1); 
        }
        
        // TODO: accepting only one option for now
        assert(argc == 0);
    }
    
    uint32_t seed = mt_goodseed();
    mt_seed32(seed);
    
    path.numParticles = 1;
    path.numTimeSlices = 64;
    path.beta = 10.0;
    path.tau = path.beta/path.numTimeSlices;
    
    double en_exact = SHOExact(path.beta);

    printf("ctx.rank = %d, size = %d\n", ctx.rank, ctx.size);

    int sockfd = 0;
    if (opt_client && (ctx.rank == 0)) {
        sockfd = initClient();
        printf("client: connection established at socket = %d\n", sockfd);

        if (sockfd > 0) {
            sendFloat64(sockfd, path.beta);
            sendInt32(sockfd, (int) path.numTimeSlices);
            sendInt32(sockfd, ctx.size);
            sendFloat64(sockfd, en_exact);
        } else {
            fprintf(stderr, "ERROR: client could not connect to server\n");
            fprintf(stderr, "Continuing calculation without communicating with the server\n\n");
        }
    }

    alloc_beads(&path);
    
    if (ctx.rank == 0) {
        printf("Simulation parameters:\n");
        printf("Number of Particles   = %zu\n", path.numParticles);
        printf("Number of Time Slices = %zu\n", path.numTimeSlices);
        printf("beta                  = %.3lf\n", path.beta);
        printf("tau                   = %.3lf\n", path.tau);
    }

    for (size_t i = 0; i < path.numTimeSlices; ++i) {
        for (size_t j = 0; j < path.numParticles; ++j) {
            path.beads[i][j] = 0.5 * (-1.0 + 2.0*mt_drand()); 
        }
    }
    
    size_t MC_steps = 100 * 1000 * 1000;
    pimc_driver(ctx, path, MC_steps, sockfd, true);

    int binSize = 500;
    Stats s = getStatsEx(trace.energies, binSize);
    printf("Collected %zu values\n", trace.energies.count); 
    
    // NOTE: the error of the mean is actually the samples std divided by the sqrt(number-of-samples)  
    double en_err = s.std/sqrt(trace.energies.count); 
    printf("(PIMC) Energy = %.5f +/- %.5f\n", s.mean, en_err); 
    
    printf("(Exact) Energy = %.5f\n", en_exact);
    
    double err = fabs(s.mean - en_exact) / en_exact;
    printf("Error: %.2f%%\n", err*100.0);
   
    /* 
    if (opt_position_histogram) {
        gsl_histogram *p_histogram;

        size_t nbins = 50;
        p_histogram = gsl_histogram_alloc(nbins);
        gsl_histogram_set_ranges_uniform(p_histogram, -1.0, 1.0);

        for (size_t i = 0; i < trace.positions.count; ++i) {
            gsl_histogram_increment(p_histogram, trace.positions.items[i]);
        }

        double q2_mean = 0.0;
        for (size_t i = 0; i < trace.positions.count; ++i) {
            q2_mean = q2_mean + trace.positions.items[i]*trace.positions.items[i];
        } 
        q2_mean /= trace.positions.count;

        // <q^2> = hbar/(2*m*omega) * coth(lambda/2), where lambda = beta*hbar*omega
        double q2_exact = 0.5/tanh(0.5*path.beta);  

        printf("(PIMC) <q^2> = %.5f\n", q2_mean); 
        printf("(Exact) <q^2> = %.5f\n", q2_exact); 
        err = fabs(q2_mean - q2_exact) / q2_exact;
        printf("Error: %.2f%%\n", err*100.0);

        // draw_histogram("Position histogram", p_histogram);
        // gsl_histogram_free(p_histogram);
    }
    */

    free(trace.energies.items);
    free(trace.positions.items);
    dealloc_beads(&path);    
    close(sockfd);
}

int main(int argc, char *argv[])
{
#ifndef NO_MPI
    MPI_Init(&argc, &argv);

    MPI_Context ctx = {0}; 
    MPI_Comm_size(MPI_COMM_WORLD, &ctx.size);
    MPI_Comm_rank(MPI_COMM_WORLD, &ctx.rank);
#else
    MPI_Context ctx = {
        .rank = 0,
        .size = 1,
    };
#endif

    char *program_path = shift(&argc, &argv);

    if (argc <= 0) {
        usage(program_path, (Subcmd*) subcmds, SUBCMDS_COUNT);
        fprintf(stderr, "ERROR: no subcommand is provided\n");
        exit(1);
    }

    const char *subcmd_id = shift(&argc, &argv);
    MPI_Subcmd *subcmd = (MPI_Subcmd*) find_subcmd_by_id((Subcmd*) subcmds, SUBCMDS_COUNT, subcmd_id);

    if (subcmd != NULL) {
        subcmd->run(ctx, argc, argv);
    } else {
        usage(program_path, (Subcmd*) subcmds, SUBCMDS_COUNT);
        fprintf(stderr, "ERROR: unknown subcommand  `%s`\n", subcmd_id);
        exit(1);
    }

    arena_free(&arena);

#ifndef NO_MPI
    MPI_Finalize();
#endif

    return 0;
}

/*
int main()
{
    gsl_histogram *h = gsl_histogram_calloc_uniform(5, 0, 5);
    gsl_histogram_increment(h, 0.5);
    gsl_histogram_increment(h, 1.5);

    gsl_histogram_fprintf(stdout, h, "%f", "%f");

    printf("------------------\n");

    h = gsl_histogram_extend(h);
    gsl_histogram_fprintf(stdout, h, "%f", "%f");

    return 0;
}
*/
