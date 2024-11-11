#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <float.h>
#include <time.h>

#ifndef NO_MPI
#include <mpi.h>
#endif

#include <gsl/gsl_histogram.h>
#include <gsl/gsl_statistics.h>

#define DA_INIT_CAP 256
#define da_append(da, item)                                                          \
    do {                                                                             \
        if ((da)->count >= (da)->capacity) {                                         \
            (da)->capacity = (da)->capacity == 0 ? DA_INIT_CAP : (da)->capacity*2;   \
            (da)->items = realloc((da)->items, (da)->capacity*sizeof(*(da)->items)); \
            assert((da)->items != NULL && "Buy more RAM lol");                       \
        }                                                                            \
                                                                                     \
        (da)->items[(da)->count++] = (item);                                         \
    } while (0)

#define return_defer(value) do { result = (value); goto defer; } while (0)

#define ALU               5.29177210903e-11 // SI: m
#define AMU               9.1093837015e-31 // SI: kg
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

#define PROTOCOL_IMPLEMENTATION 
#include "protocol.h"

// ---------------------------------------------------------
// #define m_He (4.00260325413 * RAMTOAMU)
// #define m_Ar (39.9623831237 * RAMTOAMU)
// #define mu m_He * m_Ar / (m_He + m_Ar)  // a.u.

#include "morse.h"
#define mu 12500.0 
                                         
#define lam 1.0/(2.0*mu) // hbar^2/2m

#define COORD_SAMPLE_MIN 4.0
#define COORD_SAMPLE_MAX 30.0 // a.u -- we sample coordinates within this cube ???
#define RMIN_COLLECT 4.0  // a.u.
#define RMAX_COLLECT 30.0 // a.u.
#define COORD_MAX    40.0 // a.u.
// ---------------------------------------------------------

typedef struct {
    int rank;
    int size;
} MPI_Context;

#define XC(coords, i) coords[3*i + 0]
#define YC(coords, i) coords[3*i + 1]
#define ZC(coords, i) coords[3*i + 2]

typedef struct {
    double *beads;
    size_t numTimeSlices;
    double tau;
    double beta;
    size_t steps_since_birth;
} Path;

typedef struct {
    size_t CenterOfMass;
    size_t Staging;
} AcceptanceRate;


typedef struct {
    double *items;
    size_t count;
    size_t capacity;
} DA;

typedef struct {
    DA positions;
    DA m0s;
    DA energies;
    DA necklace_sizes; 
} PIMC_Trace;

PIMC_Trace trace = {0};

typedef struct {
    void (*run)(MPI_Context ctx, int argc, char **argv);
    const char *id;
    const char *description;
} MPI_Subcmd;

void subcmd_visualize_necklace(MPI_Context ctx, int argc, char **argv);
void subcmd_visualize_ensemble(MPI_Context ctx, int argc, char **argv);
void subcmd_run(MPI_Context ctx, int argc, char **argv);

char* shift(int *argc, char ***argv);
void usage(const char *program_path, MPI_Subcmd *subcmds, size_t subcmds_count);
MPI_Subcmd *find_subcmd_by_id(MPI_Subcmd *subcmds, size_t subcmds_count, const char *id); 

void load_resources();
int msleep(long msec);


#define DEFINE_SUBCMD(name, desc) \
    {                             \
        .run = subcmd_##name,     \
        .id  = #name,             \
        .description = desc       \
    } 

MPI_Subcmd subcmds[] = {
    DEFINE_SUBCMD(visualize_necklace, "Visualize the changes single necklace upon COM & Staging actions"),
    DEFINE_SUBCMD(visualize_ensemble, "Visualize an ensemble of necklaces"), 
    DEFINE_SUBCMD(run, "Run the PIMC calculation of the harmonic oscillator"),
};
#define SUBCMDS_COUNT (sizeof(subcmds)/sizeof(subcmds[0]))

void alloc_beads(Path *path)
{
    path->beads = (double*) malloc(path->numTimeSlices * 3*sizeof(double));
    assert(path->beads != NULL);

    memset(path->beads, 0.0, path->numTimeSlices * 3*sizeof(double));
}

double sample_bead()
{
    while (true) {
        double c = COORD_SAMPLE_MAX * 0.5*(-1.0 + 2.0*mt_drand());
        if (fabs(c) > COORD_SAMPLE_MIN) return c; 
    }

}

void sample_beads(Path *path)
{
    path->steps_since_birth = 0;

    for (size_t i = 0; i < 3*path->numTimeSlices; ++i) {
        path->beads[i] = sample_bead(); 
    }
}

double PotentialAction(Path path, size_t tslice) 
{
    assert(tslice < path.numTimeSlices);

    double r2 = XC(path.beads, tslice)*XC(path.beads, tslice) + 
                YC(path.beads, tslice)*YC(path.beads, tslice) + 
                ZC(path.beads, tslice)*ZC(path.beads, tslice); 
    double r = sqrt(r2);

    // return path.tau * V_HeAr(r);
    return path.tau * morse(r);
}

double Path_PotentialAction(Path path) 
{
    double act = 0.0;

    for (size_t tslice = 0; tslice < path.numTimeSlices; ++tslice) {
        double r2 = XC(path.beads, tslice)*XC(path.beads, tslice) + 
                    YC(path.beads, tslice)*YC(path.beads, tslice) + 
                    ZC(path.beads, tslice)*ZC(path.beads, tslice); 
        double r = sqrt(r2);
        // act += V_HeAr(r);
        act += morse(r);
    }

    return path.tau*act; 
}


double Path_KineticAction(Path path) 
{
    double norm = 1.0/(4.0*lam * path.tau);
    double act = 0.0;

    for (size_t tslice = 0; tslice < path.numTimeSlices; ++tslice) {
        size_t tslicep1 = (tslice + 1) % path.numTimeSlices;

        double dx = XC(path.beads, tslicep1) - XC(path.beads, tslice);
        double dy = YC(path.beads, tslicep1) - YC(path.beads, tslice);
        double dz = ZC(path.beads, tslicep1) - ZC(path.beads, tslice);
        act += dx*dx + dy*dy + dz*dz;
    }

    return norm*act;
}

double EnergyEstimator(Path path) 
{
    double kin = 0.0;
    double norm = 1.0/(4.0*lam * path.tau*path.tau);

    for (size_t tslice = 0; tslice < path.numTimeSlices; ++tslice) {
        size_t tslicep1 = (tslice + 1) % path.numTimeSlices;

        double dx = XC(path.beads, tslicep1) - XC(path.beads, tslice);
        double dy = YC(path.beads, tslicep1) - YC(path.beads, tslice);
        double dz = ZC(path.beads, tslicep1) - ZC(path.beads, tslice);
        kin += norm*(dx*dx + dy*dy + dz*dz);
    }

    double pot = 0.0;
    for (size_t tslice = 0; tslice < path.numTimeSlices; ++tslice) {
        double r2 = XC(path.beads, tslice)*XC(path.beads, tslice) + 
                    YC(path.beads, tslice)*YC(path.beads, tslice) + 
                    ZC(path.beads, tslice)*ZC(path.beads, tslice);
        double r = sqrt(r2);
        // pot += V_HeAr(r);  
        pot += morse(r);  
    }

    // @NOTE:
    // note the 3.0/2.0 factor in the constant part of the energy estimator (D = 3)
    return 3.0/2.0/path.tau - kin/path.numTimeSlices + pot/path.numTimeSlices;
}

void apply_burnin(Path *path, size_t burnin_len);

int COM_Move(Path *path)
{
    int result = 0;

    double delta = 0.05; 
    double shiftx = delta*COORD_SAMPLE_MAX*0.5*(-1.0 + 2.0*mt_drand());
    double shifty = delta*COORD_SAMPLE_MAX*0.5*(-1.0 + 2.0*mt_drand());
    double shiftz = delta*COORD_SAMPLE_MAX*0.5*(-1.0 + 2.0*mt_drand());

    double oldAction = 0.0;
    for (size_t tslice = 0; tslice < path->numTimeSlices; ++tslice) {
        oldAction = oldAction + PotentialAction(*path, tslice);
    } 
    
    double *oldbeads = (double*) arena_alloc(&arena, path->numTimeSlices * 3*sizeof(double));
    memset(oldbeads, 0.0, path->numTimeSlices * 3*sizeof(double));

    for (size_t tslice = 0; tslice < path->numTimeSlices; ++tslice) {
        XC(oldbeads, tslice) = XC(path->beads, tslice);
        YC(oldbeads, tslice) = YC(path->beads, tslice);
        ZC(oldbeads, tslice) = ZC(path->beads, tslice);
        XC(path->beads, tslice) = XC(path->beads, tslice) + shiftx;
        YC(path->beads, tslice) = YC(path->beads, tslice) + shifty;
        ZC(path->beads, tslice) = ZC(path->beads, tslice) + shiftz;
    }

    double newAction = 0.0;
    for (size_t tslice = 0; tslice < path->numTimeSlices; ++tslice) {
        newAction = newAction + PotentialAction(*path, tslice);
    }
    
    // accept the move, or reject and restore the bead positions
    double u = mt_drand();
    double alpha = exp(-(newAction - oldAction));

    if (u < alpha) {
        bool within = true;

        for (size_t tslice = 0; tslice < path->numTimeSlices; ++tslice) {
            if ((XC(path->beads, tslice) > COORD_MAX) || (XC(path->beads, tslice) < -COORD_MAX) || 
                (YC(path->beads, tslice) > COORD_MAX) || (YC(path->beads, tslice) < -COORD_MAX) || 
                (ZC(path->beads, tslice) > COORD_MAX) || (ZC(path->beads, tslice) < -COORD_MAX)) 
            {
                within = false;
                break;
            } 
        }

        if (!within) {
            // printf("Resampling necklace\n");
            sample_beads(path);
            apply_burnin(path, 30);
        }

        path->steps_since_birth++;
        return_defer(1);
    } else {
        for (size_t tslice = 0; tslice < path->numTimeSlices; ++tslice) {
            XC(path->beads, tslice) = XC(oldbeads, tslice);
            YC(path->beads, tslice) = YC(oldbeads, tslice);
            ZC(path->beads, tslice) = ZC(oldbeads, tslice);
        }

        return_defer(0);
    }

defer:
    arena_reset(&arena);
    return result;
}

int Staging_Move(Path *path)
// http://link.aps.org/doi/10.1103/PhysRevB.31.4234
{
    // NOTE: how should we set the stage length?  
    size_t stage_len = path->numTimeSlices/2; 
    assert(stage_len < path->numTimeSlices);

    int result = 0;

    size_t alpha_start = mt_lrand() % path->numTimeSlices;
    size_t alpha_end = (alpha_start + stage_len) % path->numTimeSlices; 

    double *oldbeads = (double*) arena_alloc(&arena, path->numTimeSlices * 3*sizeof(double));
    
    double oldAction = 0.0;
    for (size_t i = 1; i < stage_len; ++i) {
        size_t tslice = (alpha_start + i) % path->numTimeSlices;

        XC(oldbeads, i - 1) = XC(path->beads, tslice);
        YC(oldbeads, i - 1) = YC(path->beads, tslice);
        ZC(oldbeads, i - 1) = ZC(path->beads, tslice);
        oldAction = oldAction + PotentialAction(*path, tslice);
    }

    double newAction = 0.0;
    for (size_t i = 1; i < stage_len; ++i) {
        size_t tslice = (alpha_start + i) % path->numTimeSlices;
        size_t tslicem1 = (tslice - 1) % path->numTimeSlices;

        double tau1 = (stage_len - i) * path->tau;

        double avx = (tau1*XC(path->beads, tslicem1) + path->tau*XC(path->beads, alpha_end))/(path->tau+tau1);
        double avy = (tau1*YC(path->beads, tslicem1) + path->tau*YC(path->beads, alpha_end))/(path->tau+tau1);
        double avz = (tau1*ZC(path->beads, tslicem1) + path->tau*ZC(path->beads, alpha_end))/(path->tau+tau1);

        double sigma2 = 2.0*lam / (1.0/path->tau + 1.0/tau1);
        
        XC(path->beads, tslice) = avx + sqrt(sigma2)*generate_normal(1.0);
        YC(path->beads, tslice) = avy + sqrt(sigma2)*generate_normal(1.0);
        ZC(path->beads, tslice) = avz + sqrt(sigma2)*generate_normal(1.0); 

        newAction = newAction + PotentialAction(*path, tslice); 
    }

    double u = mt_drand();
    double alpha = exp(-(newAction - oldAction));
  
    if (u < alpha) {
        path->steps_since_birth++;
        return_defer(1);
    } else {
        for (size_t i = 1; i < stage_len; ++i) {
            size_t tslice = (alpha_start + i) % path->numTimeSlices;
            XC(path->beads, tslice) = XC(oldbeads, i - 1);
            YC(path->beads, tslice) = YC(oldbeads, i - 1);
            ZC(path->beads, tslice) = ZC(oldbeads, i - 1);
        }

        return_defer(0);
    }
   
    return_defer(1);

defer:
    arena_reset(&arena);
    return result;

}

int Simple_MH(Path path)
{
    size_t tslice = mt_lrand() % path.numTimeSlices;

    double oldbead[3] = {0};
    oldbead[0] = XC(path.beads, tslice);
    oldbead[1] = YC(path.beads, tslice);
    oldbead[2] = ZC(path.beads, tslice);
    
    double oldAction = Path_KineticAction(path) + Path_PotentialAction(path); 

    double sigma = 10.0*path.tau*lam;
    XC(path.beads, tslice) = XC(path.beads, tslice) + sqrt(sigma)*generate_normal(1.0);
    YC(path.beads, tslice) = YC(path.beads, tslice) + sqrt(sigma)*generate_normal(1.0);
    ZC(path.beads, tslice) = ZC(path.beads, tslice) + sqrt(sigma)*generate_normal(1.0);
    
    double newAction = Path_KineticAction(path) + Path_PotentialAction(path); 

    double diffAction = newAction - oldAction;
    if (diffAction < 0) {
        return 1;
    }

    double u = mt_drand();
    double alpha = exp(-diffAction);
    if (u < alpha) {
        return 1;
    } else {
        XC(path.beads, tslice) = oldbead[0];
        YC(path.beads, tslice) = oldbead[1];
        ZC(path.beads, tslice) = oldbead[2];
        return 0;
    }
}

void apply_burnin(Path *path, size_t burnin_len) {
    
    for (size_t steps = 0; steps < burnin_len; ) {
        steps += COM_Move(path);

        if (path->numTimeSlices > 1) {
            steps += Staging_Move(path);
        }
    }
}

double compute_necklace_size(Path path)
{
    double com[3] = {0};
    for (size_t tslice = 0; tslice < path.numTimeSlices; ++tslice) {
        com[0] += XC(path.beads, tslice) / path.numTimeSlices; 
        com[1] += YC(path.beads, tslice) / path.numTimeSlices; 
        com[2] += ZC(path.beads, tslice) / path.numTimeSlices; 
    }

    double necklace_size = 0.0; 
    for (size_t tslice = 0; tslice < path.numTimeSlices; ++tslice) {
        double dx = XC(path.beads, tslice) - com[0];
        double dy = YC(path.beads, tslice) - com[1];
        double dz = ZC(path.beads, tslice) - com[2];

        necklace_size = dx*dx + dy*dy + dz*dz;
    }

    return sqrt(necklace_size)/path.numTimeSlices;
}

void gather_and_send_to_server(MPI_Context ctx, int sockfd, const char *data_name, double *send_items, size_t *send_count, size_t *packets_sent)
{
    size_t packet_size = 1000;
    assert(packet_size <= 1000 && " NOTE: keep the packet_size under 100 for now. In the local network we encounter problems with sending larger packets\n");
    
    if (*send_count >= packet_size) {
#ifndef NO_MPI
        if (ctx.rank > 0) {
            MPI_Send(send_items, packet_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            *send_count = 0;
        } else {
#endif // NO_MPI

            sendNamedFloat64Array(sockfd, data_name, send_items, *send_count);
            (*packets_sent)++;
            *send_count = 0;
#ifndef NO_MPI
            MPI_Status status = {0}; 
            for (int i = 1; i < ctx.size; ++i) {
                MPI_Recv(send_items, packet_size, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                *send_count = packet_size;
            
                sendNamedFloat64Array(sockfd, data_name, send_items, *send_count); 
                (*packets_sent)++;
                *send_count = 0;
            }

            if (*packets_sent % 100 == 0) { 
                printf("Sent %zu packets...\n", *packets_sent);
            }
        } 
#endif // NO_MPI
    }
}

void print_necklace(Path path) {
    printf("Steps since birth: %zu\n", path.steps_since_birth);

    for (size_t tslice = 0; tslice < path.numTimeSlices; ++tslice) {
        printf("x = %.3f y = %.3f z = %.3f\n", XC(path.beads, tslice), YC(path.beads, tslice), ZC(path.beads, tslice));
    }
}

void pimc_driver(MPI_Context ctx, Path path, size_t numSteps, int sockfd)
{
    (void) ctx;

    AcceptanceRate acc = {0};

    size_t observableSkip = 80;
    printf("Total MC steps: %zu\n", numSteps);
    printf("Collecting 1 out of %zu steps\n", observableSkip); 

    size_t packets_sent = 0;

    for (size_t step = 0; step < numSteps; ++step) {
        acc.CenterOfMass += COM_Move(&path);

        if (path.numTimeSlices > 1) {
            acc.Staging += Staging_Move(&path);
        }

        if (path.steps_since_birth > 0 && (path.steps_since_birth % observableSkip == 0)) {
            double com[3] = {0}; 
            for (size_t tslice = 0; tslice < path.numTimeSlices; ++tslice) {
                com[0] += XC(path.beads, tslice) / path.numTimeSlices; 
                com[1] += YC(path.beads, tslice) / path.numTimeSlices; 
                com[2] += ZC(path.beads, tslice) / path.numTimeSlices; 
            }

            double rcom = sqrt(com[0]*com[0] + com[1]*com[1] + com[2]*com[2]); 

            if ((rcom > RMIN_COLLECT) && (rcom < RMAX_COLLECT)) {
                double m0_est = 0.0;
                
                for (size_t tslice = 0; tslice < path.numTimeSlices; ++tslice) {
                    double r2 = XC(path.beads, tslice)*XC(path.beads, tslice) + 
                                YC(path.beads, tslice)*YC(path.beads, tslice) + 
                                ZC(path.beads, tslice)*ZC(path.beads, tslice); 
                    double r = sqrt(r2);

                    double dipval = dip_HeAr(r);
                    m0_est = m0_est + dipval*dipval/path.numTimeSlices; 
                }

                // da_append(&trace.positions, rcom);
                // da_append(&trace.m0s, m0_est);
                
                double en = EnergyEstimator(path);
                da_append(&trace.energies, en);
                

                double necklace_size = compute_necklace_size(path); 
                if (necklace_size > 0.2) {
                    print_necklace(path);
                    msleep(1000);
                }
                da_append(&trace.necklace_sizes, necklace_size);
            }

            gather_and_send_to_server(ctx, sockfd, "Necklace size", trace.necklace_sizes.items, &trace.necklace_sizes.count, &packets_sent);
            gather_and_send_to_server(ctx, sockfd, "Energy",        trace.energies.items, &trace.energies.count, &packets_sent);
        } 
    }

    printf("--------------------------------------\n");
    printf("Acceptance Ratios:\n");
    printf("Center of Mass : %.3lf\n", (double) acc.CenterOfMass/numSteps);
    printf("Staging        : %.3lf\n", (double) acc.Staging/numSteps);
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
Stats getStats(double *arr, size_t count); 

void subcmd_run(MPI_Context ctx, int argc, char **argv)
{
    (void) argc;
    (void) argv;

    double T = 600.0; // K
    double beta = 1.0/(Boltzmann_Hartree * T); 

    Path path = {0};
    path.numTimeSlices = 16; 
    path.tau = beta/path.numTimeSlices;
    path.beta = beta;
    alloc_beads(&path);
    
    // double coeff = pow(4.0 * M_PI * lam * path.tau, -1.5);

    // Source: M0 at 295K from diploma
    // double refVal = 5.27e-05;

    // double refVal = 22.491;   // mean(R) obtained with HEP: int(R*exp(-H/kT))/int(exp(-H/kT)) 
    // double refVal = 17.525;   // mean(R) obtained with HEP: int(R*exp(-V/kT))/int(exp(-V/kT))
    // double refVal = 2.14e-6;  // mean(mu^2) obtained with HEP: int(mu^2 exp(-V/kT))/int(exp(-V/kT)) 
    // double refVal = 1.655e-07;   // mean(mu^2) obtained with HEP [300K]: int(mu^2 exp(-H/kT))/int(exp(-H/kT)) 
    // double refVal = 2.18387e-07; // mean(mu^2) obtained with HEP [400K]: int(mu^2 exp(-H/kT))/int(exp(-H/kT)) 

    // double refVal = 7.11076e-04; // mean(E) obtained with HEP [300 K] 
    double refVal = -9.05971e-03; // mean(E) for MORSE 
    int blockSize = 1; // TODO: ignore for now

    if (ctx.rank == 0) {
        printf("Simulation parameters:\n");
        printf("Number of Time Slices = %zu\n", path.numTimeSlices);
        printf("beta                  = %.3lf\n", path.beta);
        printf("tau                   = %.3lf\n", path.tau);
    }

    sample_beads(&path);

    int sockfd = 0; 
    if (ctx.rank == 0) { 
        sockfd = initClient();
        printf("connection established at socket = %d\n", sockfd);

        if (sockfd > 0) {
            sendFloat64(sockfd, path.beta);
            sendInt32(sockfd, (int) path.numTimeSlices);
            sendInt32(sockfd, ctx.size);
            sendNamedFloat64(sockfd, "Energy", refVal);
            recvInt32(sockfd, &blockSize);
        } else {
            fprintf(stderr, "ERROR: client could not connect to server\n");
            fprintf(stderr, "Continuing calculation without communicating with the server\n\n");
        }
    }
   
    size_t MC_steps =  1 * 1000 * 1000 * 1000;
    pimc_driver(ctx, path, MC_steps, sockfd);
    
    // Stats m0 = getStats(trace.m0s.items, trace.m0s.count);
    // printf("Minimum: %.10lf\n", m0.min);
    // printf("Maximum: %.10lf\n", m0.max);
    // printf("Mean: %.3e\n", m0.mean); 
}

void update_necklace_frame(Path path)
{
    static int simulation_step = 0; 

    if (IsKeyPressed(KEY_SPACE)) {
        printf("Simulation step: %d\n", simulation_step);

        if (simulation_step % 2 == 0) {
            COM_Move(&path);
        } else {
            Staging_Move(&path);
        }

        simulation_step++;
    }

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

    DrawRectangleLinesEx(world, 3.0, LIGHTGRAY);
    double GRID_MIN = -10.0;
    double GRID_MAX =  10.0;   

    double tau = world.height / (path.numTimeSlices - 1);

    double prevx, prevy;

    Color c = GetColor(0xF2AF29FF);

    for (size_t tslice = 0; tslice < path.numTimeSlices; ++tslice) {
        double bead = ZC(path.beads, tslice);

        if ((bead < GRID_MIN) || (bead > GRID_MAX)) {
            printf("%.5lf\n", bead);
        }

        double screenx = world.x + (bead - GRID_MIN) / (GRID_MAX - GRID_MIN) * simul_sz;
        double screeny = world.y + tau * tslice; 
        DrawCircle(screenx, screeny, 6.0, c);

        if (tslice > 0) {
            DrawLine(screenx, screeny, prevx, prevy, c);
        }

        prevx = screenx;
        prevy = screeny;
    }

    EndDrawing(); 
}

void subcmd_visualize_necklace(MPI_Context ctx, int argc, char **argv)
{
    (void) ctx;
    (void) argc;
    (void) argv;

    SetConfigFlags(FLAG_WINDOW_RESIZABLE);
    SetConfigFlags(FLAG_MSAA_4X_HINT);
    
    size_t factor = 80;
    InitWindow(16*factor, 9*factor, "Visualize necklace");
    SetExitKey(KEY_Q);

    double T = 50.0; // K
    double beta = 1.0/(Boltzmann_Hartree * T); 

    Path path = {0};
    path.numTimeSlices = 4;
    path.tau = beta/path.numTimeSlices;
    path.beta = beta;
    
    alloc_beads(&path);
    // sample_beads(&path);
    for (size_t tslice = 0; tslice < path.numTimeSlices; ++tslice) {
        XC(path.beads, tslice) = 7.0 + mt_drand();
        YC(path.beads, tslice) = 7.0 + mt_drand();
        ZC(path.beads, tslice) = 7.0 + mt_drand();
    }
     
    SetTargetFPS(60);
    while (!WindowShouldClose()) {
        update_necklace_frame(path);
    }

    CloseWindow();
}
    
#define ENSEMBLE_SIZE 50 
Path ensemble[ENSEMBLE_SIZE] = {0};

void update_ensemble_frame()
{
#define CIRCLE_SIZE 3.0

    static int simulation_step = 0;
    static size_t acc = 0;

    for (size_t i = 0; i < ENSEMBLE_SIZE; ++i) {
        acc += COM_Move(&ensemble[i]);
        acc += Staging_Move(&ensemble[i]);
        // acc += Simple_MH(ensemble[i]);
    }

    simulation_step++;

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

    DrawRectangleLinesEx(world, 3.0, LIGHTGRAY);
   
    DrawCircleLines(world.x+world.width/2, world.y+world.height/2, simul_sz/(2.0*COORD_MAX)*RMAX_COLLECT, LIGHTGRAY);
    DrawCircleLines(world.x+world.width/2, world.y+world.height/2, simul_sz/(2.0*COORD_MAX)*RMIN_COLLECT, RED);

    Color color = GetColor(0xF2AF29FF);
    Color color_out = GetColor(0xA83E32FF);

    double GRID_MIN = -COORD_MAX;
    double GRID_MAX = COORD_MAX;

    int inside = 0, outside = 0;

    for (size_t i = 0; i < ENSEMBLE_SIZE; ++i)
    { 
        double com[3] = {0};
        for (size_t tslice = 0; tslice < ensemble[i].numTimeSlices; ++tslice) {
            com[0] += XC(ensemble[i].beads, tslice);
            com[1] += YC(ensemble[i].beads, tslice);
            com[2] += ZC(ensemble[i].beads, tslice);
       
            double screen_slicex = world.x + (XC(ensemble[i].beads, tslice) - GRID_MIN) / (GRID_MAX - GRID_MIN) * simul_sz;
            double screen_slicey = world.y + (YC(ensemble[i].beads, tslice) - GRID_MIN) / (GRID_MAX - GRID_MIN) * simul_sz;
            DrawCircle(screen_slicex, screen_slicey, 2.0, LIGHTGRAY);
        }
        com[0] /= ensemble[i].numTimeSlices;
        com[1] /= ensemble[i].numTimeSlices;
        com[2] /= ensemble[i].numTimeSlices;

        double screenx = world.x + (com[0] - GRID_MIN) / (GRID_MAX - GRID_MIN) * simul_sz;
        double screeny = world.y + (com[1] - GRID_MIN) / (GRID_MAX - GRID_MIN) * simul_sz;

        if ((com[0] < GRID_MIN) || (com[0] > GRID_MAX) || (com[1] < GRID_MIN) || (com[1] > GRID_MAX)) {
            DrawCircle(screenx, screeny, CIRCLE_SIZE, color_out);
            outside++;

            //for (size_t tslice = 0; tslice < ensemble[i].numTimeSlices; ++tslice) {
            //    printf("necklace[%zu] = %.5f\n", tslice, ensemble[i].beads[tslice]);
            //}

        } else { 
            double rc = sqrt(com[0]*com[0] + com[1]*com[1] + com[2]*com[2]);
            if (rc < 4.5) {
                DrawCircle(screenx, screeny, CIRCLE_SIZE, RED);
            } else {
                DrawCircle(screenx, screeny, CIRCLE_SIZE, color);
            }
            inside++;
        }
    }
        
    int font_size = 24;
    Vector2 text_pos = { .x = world.x + world.width + 50, .y = world.y };
    {
        const char *buffer = TextFormat("Simulation step: %d", simulation_step);
        DrawTextEx(font, buffer, text_pos, font_size, 0, WHITE);
        text_pos.y += font_size;
    }
    {
        const char *buffer = TextFormat("Necklaces inside: %d", inside);
        DrawTextEx(font, buffer, text_pos, font_size, 0, WHITE);
        text_pos.y += font_size;
    }
    {
        const char *buffer = TextFormat("Necklaces outside: %d", outside);
        DrawTextEx(font, buffer, text_pos, font_size, 0, WHITE);
        text_pos.y += font_size;
    }
    {
        const char *buffer = TextFormat("Acceptance rate: %.3f%%", (float)acc/simulation_step/ENSEMBLE_SIZE);
        DrawTextEx(font, buffer, text_pos, font_size, 0, WHITE);
        text_pos.y += font_size;
    }

    EndDrawing();

    msleep(100);
}

void subcmd_visualize_ensemble(MPI_Context ctx, int argc, char **argv)
{
    (void) ctx;
    (void) argc;
    (void) argv;
    
    SetConfigFlags(FLAG_WINDOW_RESIZABLE);
    SetConfigFlags(FLAG_MSAA_4X_HINT);
    
    size_t factor = 80;
    InitWindow(16*factor, 9*factor, "M0-HeAr visualize");
    SetExitKey(KEY_Q);

    load_resources();

    double T = 600.0; // K
    double beta = 1.0/(Boltzmann_Hartree*T); 

    for (size_t i = 0; i < ENSEMBLE_SIZE; ++i) {
        ensemble[i].numTimeSlices = 16;
        ensemble[i].tau = beta/ensemble[i].numTimeSlices;
        ensemble[i].beta = beta;
        
        alloc_beads(&ensemble[i]);
        // sample_beads(&ensemble[i]);
        
        for (size_t tslice = 0; tslice < ensemble[i].numTimeSlices; ++tslice) {
            XC(ensemble[i].beads, tslice) = 6.0 + mt_drand();
            YC(ensemble[i].beads, tslice) = 6.0 + mt_drand();
            ZC(ensemble[i].beads, tslice) = 6.0 + mt_drand();
        }
    }
    
    SetTargetFPS(60);
    while (!WindowShouldClose()) {
        update_ensemble_frame();
    }

    CloseWindow();
}


int main(int argc, char *argv[])
{
    MPI_Context ctx = {0};
#ifndef NO_MPI 
    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &ctx.size);
    MPI_Comm_rank(MPI_COMM_WORLD, &ctx.rank);
#else 
    ctx.size = 1;
    ctx.rank = 0;
#endif // NO_MPI 
   

    uint32_t seed = mt_goodseed();
    mt_seed32(seed);
    
    char *program_path = shift(&argc, &argv);

    if (argc <= 0) {
        usage(program_path, subcmds, SUBCMDS_COUNT);
        fprintf(stderr, "ERROR: no subcommand is provided\n");
        exit(1);
    }

    const char *subcmd_id = shift(&argc, &argv);
    MPI_Subcmd *subcmd = find_subcmd_by_id(subcmds, SUBCMDS_COUNT, subcmd_id);

    if (subcmd != NULL) {
        subcmd->run(ctx, argc, argv);
    } else {
        usage(program_path, subcmds, SUBCMDS_COUNT);
        fprintf(stderr, "ERROR: unknown subcommand  `%s`\n", subcmd_id);
        exit(1);
    }

#ifndef NO_MPI 
    MPI_Finalize();
#endif // NO_MPI
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


    return 0;
}

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

char* shift(int *argc, char ***argv)
{
    assert(*argc > 0);
    char *result = *argv[0];

    *argc -= 1;
    *argv += 1;

    return result; 
}

void usage(const char *program_path, MPI_Subcmd *subcmds, size_t subcmds_count)
{
    fprintf(stderr, "Usage: %s [subcommand]\n", program_path);
    fprintf(stderr, "Subcommands:\n");

    int width = 0;
    for (size_t k = 0; k < subcmds_count; ++k) {
        int len = strlen(subcmds[k].id);
        if (width < len) width = len;
    }

    for (size_t k = 0; k < subcmds_count; ++k) {
        fprintf(stderr, "    %-*s - %s\n", width, subcmds[k].id, subcmds[k].description);
    }
}

MPI_Subcmd *find_subcmd_by_id(MPI_Subcmd *subcmds, size_t subcmds_count, const char *id) 
{
    for (size_t k = 0; k < subcmds_count; ++k) {
        if (strcmp(subcmds[k].id, id) == 0) {
            return &subcmds[k];
        }
    }

    return NULL; 
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

int msleep(long msec)
{
    struct timespec ts = { 
        .tv_sec = msec / 1000,
        .tv_nsec = (msec % 1000) * 1000000,
    };

    int res;
    do {
        res = nanosleep(&ts, &ts);
    } while (res && errno == EINTR);

    return res;
}
// M0 (true) = 5.27e-6
// T = 300 K -> M0 (est) = 5.47e-6, N = 150 mln
