#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

#include <gsl/gsl_statistics.h>

#define ALU               5.29177210903e-11 // SI: m
#define RAMTOAMU          1822.888485332
#define Hartree           4.3597447222071e-18 // SI: J 
#define Boltzmann         1.380649e-23 // SI: J * K^(-1)
#define Boltzmann_Hartree Boltzmann/Hartree // a.u. 
#define AVOGADRO          6.022140857 * 1e23 // mol^(-1)

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

void subcmd_estimate_distribution_function(const char *program_path, int argc, char **argv);
void subcmd_estimate_energy(const char *program_path, int argc, char **argv);

Subcmd subcmds[] = {
    DEFINE_SUBCMD(estimate_energy, "collect the energy estimator during the PIMC simulation"),
    DEFINE_SUBCMD(estimate_distribution_function, "collect the diagonal term of the density matrix during the PIMC simulation"),
};
#define SUBCMDS_COUNT (sizeof(subcmds)/sizeof(subcmds[0]))

typedef struct {
    double *items;
    size_t count;
    size_t capacity;
} EnergyTrace;

// ---------------------------------------------------------
#define m_He (4.00260325413 * RAMTOAMU)
#define m_Ar (39.9623831237 * RAMTOAMU)
#define mu m_He * m_Ar / (m_He + m_Ar)  // a.u.
                                         
#define lam 0.5 // 1.0/(2.0*mu) // hbar^2/2m

#define RMAX 1.0 // a.u.
//#define TT 300.0 // K
// ---------------------------------------------------------

#define XC(coords, i) coords[3*i + 0]
#define YC(coords, i) coords[3*i + 1]
#define ZC(coords, i) coords[3*i + 2]

Path path = {0};

double V(double r) { 
    return 0.5*r*r;
}

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
        pot = pot + V(r);
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
            result = result + V(r);
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

double Energy(Path path) {
    return PotentialEnergy(path)  + KineticEnergy(path);
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
    double shiftx = delta*RMAX*(-1.0 + 2.0*mt_drand());
    double shifty = delta*RMAX*(-1.0 + 2.0*mt_drand());
    double shiftz = delta*RMAX*(-1.0 + 2.0*mt_drand());

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


EnergyTrace run_PIMC(Path path, size_t numSteps)
{
    AcceptanceRate acc = {0};

    size_t equilSkip = 20000;
    size_t observableSkip = 400;
    printf("Total MC steps: %zu\n", numSteps);
    printf("Equilibration skip: %zu\n", equilSkip);
    printf("Collecting 1 out of %zu steps\n", observableSkip); 

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

void subcmd_estimate_distribution_function(const char *program_path, int argc, char **argv)
{
    (void) program_path;
    (void) argc;
    (void) argv;
    assert(false);
}

void subcmd_estimate_energy(const char *program_path, int argc, char **argv)
{
    (void) program_path;
    (void) argc;
    (void) argv;

    uint32_t seed = mt_goodseed();
    mt_seed32(seed);

    double beta = 1.0;
        
    path.numParticles = 1;
    path.numTimeSlices = 16;
    path.tau = beta/path.numTimeSlices;
    alloc_beads(&path);

    printf("Simulation parameters:\n");
    printf("Number of Particles   = %zu\n", path.numParticles);
    printf("Number of Time Slices = %zu\n", path.numTimeSlices);
    printf("beta                  = %.3lf\n", beta);
    printf("tau                   = %.3lf\n", path.tau);

    for (size_t i = 0; i < path.numTimeSlices; ++i) {
        for (size_t j = 0; j < 3*path.numParticles; ++j) {
            path.beads[i][j] = RMAX * 0.5*(-1.0 + 2.0*mt_drand()); 
        }
    }
   
    size_t MC_steps = 10 * 1000 * 1000;
    EnergyTrace t = run_PIMC(path, MC_steps);

    double mean = gsl_stats_mean(t.items, 1, t.count);        
    double std = gsl_stats_sd_m(t.items, 1, t.count, mean);

    double omega = 1.0;
    double exact = 1.5/tanh(0.5*beta);

    printf("Collected %zu values\n", t.count); 
    printf("(PIMC) Energy = %.5f +/- %.5f\n", mean, std);
    printf("(exact) Energy = %.5e\n", exact);
    printf("Error: %.2f%%\n", (mean - exact)/exact*100);
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

