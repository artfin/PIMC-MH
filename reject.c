// ????????
// P = 2 => rel.diff: 2.13% (60mln)
// P = 2 => rel.diff: 1.29% (120mln)
// P = 2 => rel.diff: 1.13% (240mln)
// P = 2 => rel.diff: 1.06% (480mln)
// 
// P = 2 
// UPDATED proposal, constant margin: rel.diff: 0.80% (24mln)
//
//
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <assert.h>
#include <stdlib.h>

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

#define NBEADS 4

#define SIGMA 0.26
#define MARGIN 0.8
#define RMIN 4.0
#define XMAX 30.0
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

double generate_normal(double sigma) 
/*
 * Generate normally distributed variable using Box-Muller method
 */
{
    double U = mt_drand();
    double V = mt_drand();
    return sigma * sqrt(-2 * log(U)) * cos(2.0 * M_PI * V);
}

void sample(double *x)
{   
    /*
    assert(NBEADS == 2);
    double sigma = sqrt(beta/mu/NBEADS);
    printf("sigma = %.5e\n", sigma);

    double n1[3];
    n1[0] = generate_normal(sigma); 
    n1[1] = generate_normal(sigma); 
    n1[2] = generate_normal(sigma); 
    
    double n2[3];
    n2[0] = generate_normal(sigma); 
    n2[1] = generate_normal(sigma); 
    n2[2] = generate_normal(sigma); 

    x[0] = 0.5 * (n1[0] + n2[0]);
    x[1] = 0.5 * (n1[1] + n2[1]);
    x[2] = 0.5 * (n1[2] + n2[2]);
    
    x[3] = 0.5 * (n1[0] - n2[0]);
    x[4] = 0.5 * (n1[1] - n2[1]);
    x[5] = 0.5 * (n1[2] - n2[2]);
    */
    /* OPTION 1 */

    x[0] = XMAX * (2.0*mt_drand() - 1.0); // [-XMAX, XMAX]
    x[1] = XMAX * (2.0*mt_drand() - 1.0);
    x[2] = XMAX * (2.0*mt_drand() - 1.0);

    for (size_t i = 1; i < NBEADS; ++i) {
        //x[3*i + 0] = x[0] + MARGIN * (2.0*mt_drand() - 1.0); 
        //x[3*i + 1] = x[1] + MARGIN * (2.0*mt_drand() - 1.0); 
        //x[3*i + 2] = x[2] + MARGIN * (2.0*mt_drand() - 1.0);
        x[3*i + 0] = x[0] + generate_normal(SIGMA); 
        x[3*i + 1] = x[1] + generate_normal(SIGMA); 
        x[3*i + 2] = x[2] + generate_normal(SIGMA); 
    }
}


bool reject(double *x) 
{
    double u = mt_drand(); // uniform [0, 1]
   
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

int main()
{
    mt_seed32(42);
    double *x = malloc(3*NBEADS * sizeof(double)); 

    double integral = 0.0; 
    size_t accumulated = 0;
    size_t attempted = 0; 
    size_t ACC_MAX = 1000000;
            
    double V = (2*XMAX)*(2*XMAX)*(2*XMAX);

    FILE *fd = fopen("samples.txt", "w");

    while (accumulated < ACC_MAX) {
        sample(x);
        attempted++;

        if (reject(x)) {
            for (size_t i = 0; i < NBEADS; ++i) {
                double r = sqrt(x[3*i + 0]*x[3*i + 0] + x[3*i + 1]*x[3*i + 1] + x[3*i + 2]*x[3*i + 2]);
                integral += dip_HeAr(r)*dip_HeAr(r)/NBEADS; 
            }
            //fprintf(fd, "%.5e\n", x[9] - x[6]); 
            accumulated++;

            if ((accumulated > 0) && (accumulated % 1000000 == 0)) {
                printf("Accumulated %zu/%zu, acceptance rate: %.2e%%\n", accumulated, ACC_MAX, (float)accumulated/attempted*100.0);
    
                printf("Integral: %.10e, reference value: %.10e, rel.diff: %.2f%%\n\n", 
                       V*integral/accumulated, mu2_ref_value, (V*integral/accumulated - mu2_ref_value)/mu2_ref_value*100.0);
    
            }
        } 
    }
    
    integral /= accumulated;
  
    printf("---------------------------------------------------------------\n"); 
    printf("FINAL: Accumulated %zu/%zu, acceptance rate: %.2e%%\n", accumulated, ACC_MAX, (float)accumulated/attempted*100.0);

    printf("FINAL: Integral: %.10e, reference value: %.10e, rel.diff: %.2f%%\n", 
            V*integral, mu2_ref_value, (V*integral - mu2_ref_value)/mu2_ref_value*100.0);
   
    fclose(fd); 

    return 0;
}
