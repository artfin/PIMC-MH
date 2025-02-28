#include <stdio.h>
#include <stdbool.h>

#define HEAR_IMPLEMENTATION
#include "HeAr.h"

#define Hartree           4.3597447222071e-18 // SI: J 
#define Boltzmann         1.380649e-23 // SI: J * K^(-1)
#define Boltzmann_Hartree Boltzmann/Hartree // a.u. 

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

#define Temperature 300.0
#define DIM 3
#define M 1000.0
#define L 100

#define SKIP 30
#define BOUNDARY 20.0
#define COLLECT_RMIN 4.0 
#define COLLECT_RMAX 15.0 

double Hamiltonian(double *x, double *p) {
    double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]); 
    return 0.5*(p[0]*p[0] + p[1]*p[1] + p[2]*p[2])/M + V_HeAr(r)/(Boltzmann_Hartree*Temperature); 
}

void periodic_bc(double *x) 
{
    if (x[0] > BOUNDARY) x[0] = x[0] - 2*BOUNDARY; 
    if (x[0] < -BOUNDARY) x[0] = x[0] + 2*BOUNDARY;

    if (x[1] > BOUNDARY) x[1] = x[1] - 2*BOUNDARY; 
    if (x[1] < -BOUNDARY) x[1] = x[1] + 2*BOUNDARY;
    
    if (x[2] > BOUNDARY) x[2] = x[2] - 2*BOUNDARY; 
    if (x[2] < -BOUNDARY) x[2] = x[2] + 2*BOUNDARY;
}

/*
void periodic_bc(double *x)
{
    double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]); 
    if (r > BOUNDARY) {
        x[0] = -x[0];
        x[1] = -x[1];
        x[2] = -x[2];
    }
}
*/

bool HMC_Move(double *x)
{
    double oldx[DIM];
    oldx[0] = x[0];
    oldx[1] = x[1];
    oldx[2] = x[2]; 

    double p[3];
    p[0] = generate_normal(M); 
    p[1] = generate_normal(M); 
    p[2] = generate_normal(M);
    
    double energy0 = Hamiltonian(x, p); 

    // leapfrog algorithm
    double deltaT = 1.0e-2; 

    for (size_t step = 0; step < L; ++step) {
        double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]); 
        double dV = dV_HeAr(r);

        double p_halfstep[3];
        p_halfstep[0] = p[0] - deltaT/2.0*dV;
        p_halfstep[1] = p[1] - deltaT/2.0*dV;
        p_halfstep[2] = p[2] - deltaT/2.0*dV;

        x[0] = x[0] + deltaT/M*p_halfstep[0];
        x[1] = x[1] + deltaT/M*p_halfstep[1];
        x[2] = x[2] + deltaT/M*p_halfstep[2];

        r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]); 
        dV = dV_HeAr(r);
        p[0] = p_halfstep[0] - deltaT/2.0*dV; 
        p[1] = p_halfstep[1] - deltaT/2.0*dV; 
        p[2] = p_halfstep[2] - deltaT/2.0*dV; 

        periodic_bc(x); 
    }
   
    double energyL = Hamiltonian(x, p); 
    
    double u = mt_drand(); // uniform [0, 1]
    double alpha = exp(-(energyL - energy0));
    if (u > alpha) {
        // rejected 
        x[0] = oldx[0]; 
        x[1] = oldx[1]; 
        x[2] = oldx[2];
        return false; 
    }

    return true;
}

void sample(double *x)
{
    x[0] = BOUNDARY * 0.5*(-1.0 + 2.0*mt_drand()); // [-BOUNDARY, BOUNDARY] 
    x[1] = BOUNDARY * 0.5*(-1.0 + 2.0*mt_drand()); // [-BOUNDARY, BOUNDARY] 
    x[2] = BOUNDARY * 0.5*(-1.0 + 2.0*mt_drand()); // [-BOUNDARY, BOUNDARY] 
}


int main()
{
    size_t accepted = 0;
    size_t attempted = 10000000;
    size_t accumulated_counter = 0;
    double integral = 0.0;

    double x[3];
    sample(x);

    FILE *fd = fopen("samples.txt", "w");

    for (size_t i = 0; i < attempted; ++i) {
        accepted += HMC_Move(x);
    
        periodic_bc(x);

        if ((i > 0) && (i % SKIP == 0)) {
            double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);

            if ((r > COLLECT_RMIN) && (r < COLLECT_RMAX)) {
                fprintf(fd, "%.5e %.5e %.5e %.5e\n", x[0], x[1], x[2], r);
                integral += dip_HeAr(r)*dip_HeAr(r); 
                accumulated_counter++;
            }     
        }

        if (i % 20000 == 0) {
            printf("Attempted steps: %zu/%zu\n", i, attempted);
        } 
    }

    integral /= accumulated_counter;

    fclose(fd);
    
    double V = 4.0/3.0*M_PI*30.0*30.0*30.0;
    double mu2_ref_value = V * 1.655e-07; // mean(mu^2) for HeAr [300K]: int(mu^2 exp(-H/kT))/int(exp(-H/kT))
    //double mu2_ref_value = V * 1.32216e-07; // mean(mu^2) for HeAr [10K]: int(mu^2 exp(-H/kT))/int(exp(-H/kT))
 
    printf("Accumulated: %zu\n", accumulated_counter);

    V = 4.0/3.0*M_PI*COLLECT_RMAX*COLLECT_RMAX*COLLECT_RMAX;
    printf("Integral: %.10e, reference value: %.10e, rel.diff: %.2f%%\n", 
            V*integral, mu2_ref_value, (V*integral - mu2_ref_value)/mu2_ref_value*100.0);

    printf("Accepted: %.2f%%\n", (float)accepted/attempted*100.0);

    return 0;
}
