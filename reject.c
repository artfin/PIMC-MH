// naive version of proposal:
//   P = 2 => rel.diff: 2.13% (60mln)
//   P = 2 => rel.diff: 1.29% (120mln)
//   P = 2 => rel.diff: 1.13% (240mln)
//   P = 2 => rel.diff: 1.06% (480mln)
// 
// 300K: P = 2
//   SIGMA = 0.20
// acceptance rate: 89% 
// UPDATED proposal, uniform margin (0.8): rel.diff: 0.80% (24mln)
//                   normal  margin: 1.8885e-02, rel.diff: 0.90% (400mln)  
//                                   1.8896e-02, rel.diff: 0.95% (1600mln)  
//                                   1.8920e-02, rel.diff: 1.08% (6400mln)
//                                   1.8907e-02, rel.diff: 1.02% (25.6bln)
//                                   1.8912e-02, rel.diff: 1.04% (128bln) -- 0.013% error
// 100K: P = 2
//   SIGMA = 0.35
//                    normal margin: 7.6177e-03 (30bln) -- 0.02%
//
// 50 K: P = 2
//   SIGMA = 0.50
//                    normal margin: 5.3194e-03 (14bln) -- 0.02% 
//    
// -------------------------------------------------------------------------------------
// 300K: P = 4
//   SIGMA = 0.26
// acceptance rate: 1.03%
// UPDATED proposal, normal margin:  1.8990e-02, rel.diff: 1.46% (2000mln)
//                                   1.8970e-02, rel.diff: 1.35% (7040mln) -- 0.05% error 
//
// 100K: P = 4
//   SIGMA = 0.47
//                   normal magin:   7.718137e-03 (1260mln) -- 0.1% error
//                                   7.702722e-03 (4680mln) -- 0.05% error
//
// 50K: P = 4
//   SIGMA = 0.65
//                   normal magin:   5.4429e-03 (4500mln) -- 0.04% error
//
// -------------------------------------------------------------------------------------
// 300K: P = 8
//   SIGMA = 0.30 -- probably will work 
// acceptance rate: 3e-6%
            
// Pooled variance: https://en.wikipedia.org/wiki/Pooled_variance 

#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <assert.h>
#include <stdlib.h>
#include <inttypes.h>
#include <string.h>
#include <mpi.h>

#define HEAR_IMPLEMENTATION
#include "HeAr.h"

#define Hartree           4.3597447222071e-18 // SI: J 
#define Boltzmann         1.380649e-23 // SI: J * K^(-1)
#define Boltzmann_Hartree Boltzmann/Hartree // a.u. 

#define RAMTOAMU 1822.888485332
#define m_He (4.00260325413 * RAMTOAMU)
#define m_Ar (39.9623831237 * RAMTOAMU)
#define mu m_He * m_Ar / (m_He + m_Ar)  // a.u.

#define Rmin 6.597835932201344e+00

#define NBEADS 16

// Sampling scheme outline:
// The construction of the proposal distribution is the following:
//   first, we sample the zeroth bead with uniform distribution in the volume cube
//   then we select sample the next P-1 beads with normal distribution with 
//   variance set to SIGMA^2. 
// Then, in the 'reject' function we compute the ratio of target distribution density
//    exp(-sum_{i = 1}^{P}(m*P/(2*beta*hbar^2)*(x[i+1] - x[i])^2 + beta/P*V(x_i))
// to proposal distribution density:
//    exp(-(x[i]-x[0])^2/2*sigma^2) * exp(-beta*Vmin/P) 
// In rejection scheme, this ratio gives us the probability to choose the current point as the point
// of target distribution. If ratio turns out to be >= 1, it means that parameter SIGMA needs to
// be increased.   

#define SIGMA 0.25 // the width of normally distributed for 'dx' around zeroth bead
#define MARGIN 0.8 // the width of uniform distribution for 'dx' around zeroth bead
#define RMIN 4.0   
#define XMAX 30.0
#define VOLUME (2*XMAX)*(2*XMAX)*(2*XMAX)
#define Temperature 50.0 
#define beta 1.0/(Boltzmann_Hartree * Temperature)

#define MCMC_KSI_SIGMA 0.02
#define MCMC_KSI_BURNIN 1000
#define MCMC_KSI_SKIP 20

#define REFVALUE (4.0/3.0*M_PI*30.0*30.0*30.0 * 1.655e-07) // mean(mu^2) for HeAr [300K]: int(mu^2 exp(-H/kT))/int(exp(-H/kT)) 
//double mu2_ref_value = 4.0/3.0*M_PI*30.0*30.0*30.0 * 1.32216e-07; // mean(mu^2) for HeAr [10K]: int(mu^2 exp(-H/kT))/int(exp(-H/kT))

typedef enum {
    CALCULATION_M0,
    CALCULATION_M1,
} CALCULATION_TYPE;

CALCULATION_TYPE ct = CALCULATION_M0;

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

double generate_normal(mt_state* mts, double sigma) 
/*
 * Generate normally distributed variable using Box-Muller method
 */
{
    double U = mts_drand(mts);
    double V = mts_drand(mts);
    return sigma * sqrt(-2 * log(U)) * cos(2.0 * M_PI * V);
}

void sample(mt_state* mts, double *x)
{   
    x[0] = XMAX * (2.0*mts_drand(mts) - 1.0); // [-XMAX, XMAX]
    x[1] = XMAX * (2.0*mts_drand(mts) - 1.0);
    x[2] = XMAX * (2.0*mts_drand(mts) - 1.0);

    for (size_t i = 1; i < NBEADS; ++i) {
        //x[3*i + 0] = x[0] + MARGIN * (2.0*mt_drand() - 1.0); 
        //x[3*i + 1] = x[1] + MARGIN * (2.0*mt_drand() - 1.0); 
        //x[3*i + 2] = x[2] + MARGIN * (2.0*mt_drand() - 1.0);
        x[3*i + 0] = x[0] + generate_normal(mts, SIGMA); 
        x[3*i + 1] = x[1] + generate_normal(mts, SIGMA); 
        x[3*i + 2] = x[2] + generate_normal(mts, SIGMA); 
    }
}


double density(double *ksi) {
    double C = mu*NBEADS/(2.0*beta);
    double s = -ksi[3*0]*ksi[3*0] - ksi[3*(NBEADS-2)]*ksi[3*(NBEADS-2)];
    for (size_t i = 0; i < NBEADS - 2; ++i) {
        s -= (ksi[3*(i + 1)] - ksi[3*i]) * (ksi[3*(i + 1)] - ksi[3*i]); 
    }

    return exp(C*s); 
    
    // return exp(C * (-ksi[3*0]*ksi[3*0] - (ksi[3*2]-ksi[3*1])*(ksi[3*2]-ksi[3*1]) - ksi[3*2]*ksi[3*2] - (ksi[3*1]-ksi[3*0])*(ksi[3*1]-ksi[3*0])));
}

bool make_step(mt_state* mts, double *ksi)
{
    double oldksi[3*(NBEADS - 1)];
    for (size_t i = 0; i < NBEADS - 1; ++i) {
        oldksi[3 * i] = ksi[3 * i];
    }

    for (size_t i = 0; i < NBEADS-1; ++i) {
        ksi[3 * i] = ksi[3 * i] + generate_normal(mts, MCMC_KSI_SIGMA); 
    }

    double old_dens = density(oldksi);
    double new_dens = density(ksi);

    double alpha = new_dens/old_dens;
    double u = mts_drand(mts);

    if (u > alpha) {
        // reject 
        for (size_t i = 0; i < NBEADS - 1; ++i) {
            ksi[3 * i] = oldksi[3 * i];
        }
        return false; 
    }

    return true;
}

void sample_with_mcmc(mt_state *mts, double *x, double *ksi, bool first)
{
    if (first) {
        size_t acc = 0;

        memset(ksi, 0, 3*(NBEADS-1)*sizeof(double));
    
        for (size_t i = 0; i < MCMC_KSI_BURNIN; ++i) {
            acc += make_step(mts, &ksi[0]);
            acc += make_step(mts, &ksi[1]);
            acc += make_step(mts, &ksi[2]);
        }

        printf("MCMC acceptance rate: %.3e%%\n", (float)acc/3.0/MCMC_KSI_BURNIN * 100.0); 
    } else { 
        x[0] = XMAX * (2.0*mts_drand(mts) - 1.0); // [-XMAX, XMAX]
        x[1] = XMAX * (2.0*mts_drand(mts) - 1.0);
        x[2] = XMAX * (2.0*mts_drand(mts) - 1.0);

        for (size_t i = 0; i < MCMC_KSI_SKIP; ++i) {
            make_step(mts, &ksi[0]);
            make_step(mts, &ksi[1]);
            make_step(mts, &ksi[2]);
        }

        for (size_t i = 1; i < NBEADS; ++i) {
           x[3*i + 0] = ksi[3*(i-1) + 0] + x[0]; 
           x[3*i + 1] = ksi[3*(i-1) + 1] + x[1]; 
           x[3*i + 2] = ksi[3*(i-1) + 2] + x[2]; 
        }
    }
}

bool reject(mt_state* mts, double *x) 
{
    double u = mts_drand(mts); // uniform [0, 1]
   
    double Vs = 0.0; 
    for (size_t i = 0; i < NBEADS; ++i) {
        double r = sqrt(x[3*i + 0]*x[3*i + 0] + x[3*i + 1]*x[3*i + 1] + x[3*i + 2]*x[3*i + 2]);
        if (r < RMIN) return false;

        Vs += beta*V_HeAr(r)/NBEADS;
    }

    double C = mu*NBEADS/(2.0*beta);
    double kin = 0.0;
    for (size_t i = 1; i < NBEADS; ++i) {
        kin = kin + C*((x[3*i + 0] - x[3*(i-1) + 0])*(x[3*i + 0] - x[3*(i-1) + 0]) + \
                       (x[3*i + 1] - x[3*(i-1) + 1])*(x[3*i + 1] - x[3*(i-1) + 1]) + \
                       (x[3*i + 2] - x[3*(i-1) + 2])*(x[3*i + 2] - x[3*(i-1) + 2]));
    }
    kin = kin + C*((x[0] - x[3*(NBEADS-1) + 0])*(x[0] - x[3*(NBEADS-1) + 0]) + \
                   (x[1] - x[3*(NBEADS-1) + 1])*(x[1] - x[3*(NBEADS-1) + 1]) + \
                   (x[2] - x[3*(NBEADS-1) + 2])*(x[2] - x[3*(NBEADS-1) + 2]));

    double desired = exp(-(kin + Vs));  
    
    // double proposal = 1.0; // 1.0/(2.0*XMAX)*1.0/(2.0*MARGIN) * 1.0/(2.0*XMAX)*1.0/(2.0*MARGIN) * 1.0/(2.0*XMAX)*1.0/(2.0*MARGIN);
    double proposal = 1.0; // 1.0/(2.0*SIGMA*SIGMA)*1.0/(2.0*SIGMA*SIGMA)*1.0/(2.0*SIGMA*SIGMA);
    for (size_t i = 1; i < NBEADS; ++i) {
        proposal *= exp(-(x[3*i + 0] - x[0])*(x[3*i + 0] - x[0]) / (2.0*SIGMA*SIGMA));
        proposal *= exp(-(x[3*i + 1] - x[1])*(x[3*i + 1] - x[1]) / (2.0*SIGMA*SIGMA));
        proposal *= exp(-(x[3*i + 2] - x[2])*(x[3*i + 2] - x[2]) / (2.0*SIGMA*SIGMA));
    }
  
    double Vmin = V_HeAr(Rmin); 
    double M = exp(-beta*Vmin);
    double alpha = desired/M/proposal;

    if (alpha > 1.0) {
        printf("x[0] = %.5e, x[1] = %.5e, x[2] = %.5e\n", x[0], x[1], x[2]);
        printf("x[3] = %.5e, x[4] = %.5e, x[5] = %.5e\n", x[3], x[4], x[5]);
        printf("proposal = %.10e\n", proposal);
        printf("desired = %.10e\n", desired);
        printf("M = %.10e\n", M);

        assert(false);
    }

    return u < alpha;
} 

bool reject_with_mcmc(mt_state* mts, double *x) 
{
    double u = mts_drand(mts); // uniform [0, 1]
   
    double Vs = 0.0; 
    for (size_t i = 0; i < NBEADS; ++i) {
        double r = sqrt(x[3*i + 0]*x[3*i + 0] + x[3*i + 1]*x[3*i + 1] + x[3*i + 2]*x[3*i + 2]);
        if (r < RMIN) return false;

        Vs += beta*V_HeAr(r)/NBEADS;
    }

    double desired_proposal_ratio = exp(-Vs);  
    
    double Vmin = V_HeAr(Rmin); 
    double M = exp(-beta*Vmin);
    double alpha = desired_proposal_ratio/M;

    if (alpha > 1.0) {
        printf("x[0] = %.5e, x[1] = %.5e, x[2] = %.5e\n", x[0], x[1], x[2]);
        printf("x[3] = %.5e, x[4] = %.5e, x[5] = %.5e\n", x[3], x[4], x[5]);
        printf("desired_proposal_ratio = %.10e\n", desired_proposal_ratio);
        printf("M = %.10e\n", M);

        assert(false);
    }

    return u < alpha;
} 


void make_estimate(mt_state *mts, size_t npoints, double *mean, double *var)
// This function will return the mean and variance of sampled sequence.
// The variance of mean is 'npoints' times smaller: sigma_mu^2 = sigma^2/N, 
// this division is performed by caller of the function. 
{
    int wrank, wsize;
    MPI_Comm_size(MPI_COMM_WORLD, &wsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);

    size_t accumulated = 0;
    size_t attempted = 0; 
    
    double *x = malloc(3*NBEADS * sizeof(double)); 
    
    FILE *fd = fopen("samples.txt", "w");
    
    for (attempted = 0; accumulated < npoints; attempted++) {
        sample(mts, x);

        if (reject(mts, x)) {
            double fval = 0.0;
            for (size_t i = 0; i < NBEADS; ++i) {
                double r = sqrt(x[3*i + 0]*x[3*i + 0] + x[3*i + 1]*x[3*i + 1] + x[3*i + 2]*x[3*i + 2]);

                switch (ct) {
                    case CALCULATION_M0: {
                        fval = fval + VOLUME*dip_HeAr(r)*dip_HeAr(r)/NBEADS;
                        break;
                    }
                    case CALCULATION_M1: {
                        double dip = dip_HeAr(r);
                        double dip_deriv = ddip_HeAr(r);
                        fval = fval + VOLUME*(dip_deriv*dip_deriv + 2.0/r/r*dip*dip)/NBEADS; 
                        break;
                    }
                } 
            }

            double diff = fval - *mean;
            *mean += diff / (accumulated + 1.0);
            *var += diff * diff * (accumulated / (accumulated + 1.0));
            accumulated++;
            
            // fprintf(fd, "%.5e %.5e %.5e %.5e\n", x[3] - x[0], x[6] - x[3], x[9] - x[6], x[12] - x[9]); 

            //if ((accumulated > 0) && (accumulated % 10 == 0)) {
            //    printf("[%d] Accumulated %zu/%zu, acceptance rate: %.2e%%\n", 
            //            wrank, accumulated, npoints, (float)accumulated/attempted*100.0);
            //    printf("[%d] Integral: %.6e +/- %.6e, reference value: %.10e, rel.diff: %.2f%%\n\n", 
            //           wrank, VOLUME * *mean, VOLUME*sqrt(*var / accumulated / (accumulated - 1.0)), 
            //           mu2_ref_value, (VOLUME * *mean - mu2_ref_value)/mu2_ref_value*100.0);
            // }
        } 
    }

    *var = *var/(accumulated-1.0);

    free(x);
    fclose(fd); 
}

void make_estimate_with_mcmc_proposal(mt_state *mts, size_t npoints, double *mean, double *var)
// This function will return the mean and variance of sampled sequence.
// The variance of mean is 'npoints' times smaller: sigma_mu^2 = sigma^2/N, 
// this division is performed by caller of the function. 
{
    int wrank, wsize;
    MPI_Comm_size(MPI_COMM_WORLD, &wsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);

    size_t accumulated = 0;
    size_t attempted = 0; 
    
    FILE *fd = fopen("samples.txt", "w");

    double x[3*NBEADS]; 
    double ksi[3*(NBEADS-1)]; 
    sample_with_mcmc(mts, x, ksi, true);
    
    for (attempted = 0; accumulated < npoints; attempted++) {
        sample_with_mcmc(mts, x, ksi, false);

        if (reject_with_mcmc(mts, x)) {
            double fval = 0.0;
            for (size_t i = 0; i < NBEADS; ++i) {
                double r = sqrt(x[3*i + 0]*x[3*i + 0] + x[3*i + 1]*x[3*i + 1] + x[3*i + 2]*x[3*i + 2]);

                switch (ct) {
                    case CALCULATION_M0: {
                        fval = fval + VOLUME*dip_HeAr(r)*dip_HeAr(r)/NBEADS;
                        break;
                    }
                    case CALCULATION_M1: {
                        double dip = dip_HeAr(r);
                        double dip_deriv = ddip_HeAr(r);
                        fval = fval + VOLUME*(dip_deriv*dip_deriv + 2.0/r/r*dip*dip)/NBEADS; 
                        break;
                    }
                } 
            }

            double diff = fval - *mean;
            *mean += diff / (accumulated + 1.0);
            *var += diff * diff * (accumulated / (accumulated + 1.0));
            accumulated++;
            
            //fprintf(fd, "%.5e %.5e %.5e\n", x[3] - x[0], x[6] - x[3], x[9] - x[6]); 

            if ((accumulated > 0) && (accumulated % 100000 == 0)) {
                printf("[%d] Accumulated %zu/%zu, acceptance rate: %.2e%%\n", 
                        wrank, accumulated, npoints, (float)accumulated/attempted*100.0);
            //    printf("[%d] Integral: %.6e +/- %.6e, reference value: %.10e, rel.diff: %.2f%%\n\n", 
            //           wrank, VOLUME * *mean, VOLUME*sqrt(*var / accumulated / (accumulated - 1.0)), 
            //           mu2_ref_value, (VOLUME * *mean - mu2_ref_value)/mu2_ref_value*100.0);
            }
        } 
    }

    *var = *var/(accumulated-1.0);

    fclose(fd); 
}


int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    int wrank, wsize;
    MPI_Comm_size(MPI_COMM_WORLD, &wsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
    
    mt_state mts = {0};
    uint32_t seed = mts_goodseed(&mts);
    printf("Initialized process (%d) with seed = %" PRIu32 "\n", wrank, seed); 
    //int seeds[10] = {
    //    5856844, 53218983, 91452731, 43166694, 49914920, 31614602
    //}; 
    //assert(wsize <= 10);
    //mts_seed32(&mts, seeds[wrank]);
    //printf("[%d] INFO: seed = %d\n", wrank, seeds[wrank]);

    double total_mean = 0.0; 
    double total_var = 0.0;
    size_t total_samples_count = 0; // @todo: check that we don't overflow this at any point
    size_t packets = 0; 

    //size_t ACC_MAX = 400000000;
    size_t iters = 100;
    size_t packet_samples_count = 1000000; 
           
    if (wrank == 0) { 
        printf("INFO: the requested calculation needs to collect = %zu, SIZE_MAX = %zu\n", 
                iters*packet_samples_count*wsize, SIZE_MAX);
        assert(iters*packet_samples_count*wsize < SIZE_MAX);

        printf("NBEADS = %d, Temperature = %.3f\n", NBEADS, Temperature);
    }

    for (size_t iter = 0; iter < iters; ++iter) {
        double packet_mean = 0.0; 
        double packet_var = 0.0; 
        make_estimate_with_mcmc_proposal(&mts, packet_samples_count, &packet_mean, &packet_var);

        if (wrank == 0) {
            double packet_stdev = sqrt(packet_var/packet_samples_count);
            printf("[0] Adding packet (%zu/%zu) from %d = (%.5e +/- %.5e)\n", 
                    packets, iters*wsize, 0, packet_mean, packet_stdev);
            packets++;

            total_var = ((total_samples_count-1)*total_var + (packet_samples_count-1)*packet_var)/(total_samples_count-1 + packet_samples_count-1);
            total_mean = (total_mean*total_samples_count + packet_samples_count*packet_mean) / (total_samples_count + packet_samples_count);
            total_samples_count += packet_samples_count;

            MPI_Status status;
            for (int i = 1; i < wsize; ++i) {
                MPI_Recv(&packet_mean, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
                MPI_Recv(&packet_var, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);

                double packet_stdev = sqrt(packet_var/packet_samples_count);
                printf("[0] Adding packet (%zu/%zu) from %d = (%.5e +/- %.5e)\n", 
                        packets, iters*wsize, status.MPI_SOURCE, packet_mean, packet_stdev);
                packets++;

                total_var = ((total_samples_count-1)*total_var + (packet_samples_count-1)*packet_var)/(total_samples_count-1 + packet_samples_count-1);
                total_mean = (total_mean*total_samples_count + packet_samples_count*packet_mean) / (total_samples_count + packet_samples_count);
                total_samples_count += packet_samples_count;
            }

            double total_stdev = sqrt(total_var/total_samples_count);
            printf("[0] Aggregate estimate %zu/%zu: %.6e +/- %.6e (%.2f%%)\n", 
                    packets*packet_samples_count, iters*wsize*packet_samples_count,
                    total_mean, total_stdev, total_stdev/total_mean*100.0);
        } else {
            MPI_Send(&packet_mean, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&packet_var, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
    } 

    if (wrank == 0) {
        double total_stdev = sqrt(total_var);
        printf("---------------------------------------------------------------\n"); 
        //printf("FINAL: Accumulated %zu/%zu, acceptance rate: %.2e%%\n", accumulated, ACC_MAX, (float)accumulated/attempted*100.0);
        printf("FINAL: Integral: %.6e +/- %.6e (%.2f%%), reference value: %.10e, rel.diff: %.2f%%\n", 
                total_mean, total_stdev, total_stdev/total_mean*100.0, REFVALUE, (total_mean - REFVALUE)/REFVALUE*100.0);
    }

    MPI_Finalize();

    return 0;
}
