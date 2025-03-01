// ????????
// P = 2 => rel.diff: 2.13% (60mln)
// P = 2 => rel.diff: 1.29% (120mln)
// P = 2 => rel.diff: 1.13% (240mln)
// P = 2 => rel.diff: 1.06% (480mln)
// 
// P = 2 
// UPDATED proposal, uniform margin: rel.diff: 0.80% (24mln)
//                   normal  margin: 1.8885e-02, rel.diff: 0.9% (400mln)  
//                                   1.8896e-02, rel.diff: 0.95% (1600mln)  
//                                   1.8920e-02, rel.diff: 1.08% (6400mln) 
/*
 * NOTES on the parameter values:
 *  for P = 4 at 300K
 *   SIGMA = 0.26
 *   MARGIN = 0.8
 *  for P = 2 at 300K
 *   SIGMA = 0.20
 *   MARGIN = 0.8
 */

#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <assert.h>
#include <stdlib.h>
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

#define NBEADS 2

#define SIGMA 0.20
#define MARGIN 0.8
#define RMIN 4.0
#define XMAX 30.0
#define VOLUME (2*XMAX)*(2*XMAX)*(2*XMAX)
#define Temperature 300.0 
#define beta 1.0/(Boltzmann_Hartree * Temperature)
    
double mu2_ref_value = 4.0/3.0*M_PI*30.0*30.0*30.0 * 1.655e-07; // mean(mu^2) for HeAr [300K]: int(mu^2 exp(-H/kT))/int(exp(-H/kT))
//double mu2_ref_value = 4.0/3.0*M_PI*30.0*30.0*30.0 * 1.32216e-07; // mean(mu^2) for HeAr [10K]: int(mu^2 exp(-H/kT))/int(exp(-H/kT))


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

void make_estimate(mt_state *mts, size_t npoints, double *m_integral, double *q_integral)
{
    int wrank, wsize;
    MPI_Comm_size(MPI_COMM_WORLD, &wsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);

    size_t accumulated = 0;
    size_t attempted = 0; 
    
    double *x = malloc(3*NBEADS * sizeof(double)); 
    
    
    for (attempted = 0; accumulated < npoints; attempted++) {
        sample(mts, x);

        if (reject(mts, x)) {
            double fval = 0.0;
            for (size_t i = 0; i < NBEADS; ++i) {
                double r = sqrt(x[3*i + 0]*x[3*i + 0] + x[3*i + 1]*x[3*i + 1] + x[3*i + 2]*x[3*i + 2]);
                fval = fval + VOLUME*dip_HeAr(r)*dip_HeAr(r)/NBEADS;
            }

            double diff = fval - *m_integral;
            *m_integral += diff / (accumulated + 1.0);
            *q_integral += diff * diff * (accumulated / (accumulated + 1.0));
            accumulated++;
            
            //fprintf(fd, "%.5e\n", x[3] - x[0]); 

            //if ((accumulated > 0) && (accumulated % 1000000 == 0)) {
            //    printf("[%d] Accumulated %zu/%zu, acceptance rate: %.2e%%\n", 
            //            wrank, accumulated, npoints, (float)accumulated/attempted*100.0);
            //    printf("[%d] Integral: %.6e +/- %.6e, reference value: %.10e, rel.diff: %.2f%%\n\n", 
            //           wrank, VOLUME * *m_integral, VOLUME*sqrt(*q_integral / accumulated / (accumulated - 1.0)), 
            //           mu2_ref_value, (VOLUME * *m_integral - mu2_ref_value)/mu2_ref_value*100.0);
    
            //}
        } 
    }

    *q_integral = sqrt(*q_integral/accumulated/(accumulated-1.0));

    free(x);
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    int wrank, wsize;
    MPI_Comm_size(MPI_COMM_WORLD, &wsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
    
    mt_state mts = {0};
    int seeds[10] = {
        5856844, 53218983, 91452731, 43166694, 49914920, 31614602
    }; 
    assert(wsize <= 10);
    mts_seed32(&mts, seeds[wrank]);
    printf("[%d] INFO: seed = %d\n", wrank, seeds[wrank]);

    double tm_integral = 0.0; // total mean
    double tq_integral = 0.0; // total stdev         
    size_t packets = 0; 

    // FILE *fd = fopen("samples.txt", "w");
    //size_t ACC_MAX = 400000000;
    size_t iters = 160;
    size_t npoints_per_iter = 10000000; 

    for (size_t iter = 0; iter < iters; ++iter) {
        double m_integral = 0.0; // mean
        double q_integral = 0.0; // standard deviation
        make_estimate(&mts, npoints_per_iter, &m_integral, &q_integral);

        if (wrank == 0) {
            tm_integral += m_integral;
            tq_integral += q_integral;
            packets++;

            MPI_Status status;
            for (int i = 1; i < wsize; ++i) {
                MPI_Recv(&m_integral, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
                MPI_Recv(&q_integral, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
                printf("[0] Received packet (%zu/%zu) from %d = (%.5e +/- %.5e)\n", 
                        packets, iters*wsize, status.MPI_SOURCE, m_integral, q_integral);

                tm_integral += m_integral;
                tq_integral += q_integral;
                packets++;
            }

            printf("[0] Aggregate estimate %zu/%zu: %.5e +/- %.5e\n", 
                    packets*npoints_per_iter, iters*wsize*npoints_per_iter,
                    tm_integral/packets, tq_integral/packets);
        } else {
            MPI_Send(&m_integral, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&q_integral, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
    } 

    if (wrank == 0) {
        tm_integral /= packets;
        tq_integral /= packets;
        printf("---------------------------------------------------------------\n"); 
        //printf("FINAL: Accumulated %zu/%zu, acceptance rate: %.2e%%\n", accumulated, ACC_MAX, (float)accumulated/attempted*100.0);
        printf("FINAL: Integral: %.6e +/- %.6e, reference value: %.10e, rel.diff: %.2f%%\n", 
                tm_integral, sqrt(tq_integral / packets*npoints_per_iter / (packets*npoints_per_iter - 1.0)), 
                mu2_ref_value, (tm_integral - mu2_ref_value)/mu2_ref_value*100.0);
    }

    // fclose(fd); 
    MPI_Finalize();

    return 0;
}
