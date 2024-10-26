#include <cmath>
#include <iomanip>
#include <hep/mc-mpi.hpp>

#define HEAR_IMPLEMENTATION
#include "HeAr.h"

#define DIM 6
#define IR 0
#define IPR 1
#define IPHI 2
#define IPPHI 3
#define ITHETA 4
#define IPTHETA 5

const double Rmin = 4.0;
const double Rmax = 30.0;

const double RAMTOAMU          = 1822.888485332;
const double Hartree           = 4.3597447222071e-18; // SI: J
const double Boltzmann         = 1.380649e-23; // SI: J * K^(-1)
const double Boltzmann_Hartree = Boltzmann/Hartree; // a.u.
const double HTOJ              = 4.35974417e-18;
const double ALU               = 5.29177210903e-11; // SI: m
const double AMU               = 9.10938356e-31;

const double m_He = 4.00260325413 * RAMTOAMU;
const double m_Ar = 39.9623831237 * RAMTOAMU;
const double mu   = m_He * m_Ar / (m_He + m_Ar);  // a.u.

const double Temperature = 300.0; // K

double Hamiltonian(double qp[])
{
    double sinTheta = std::sin(qp[ITHETA]);

    return qp[IPR]*qp[IPR]/2.0/mu + \
           qp[IPTHETA]*qp[IPTHETA]/2.0/mu/qp[IR]/qp[IR] + \
           qp[IPPHI]*qp[IPPHI]/2.0/mu/qp[IR]/qp[IR]/sinTheta/sinTheta + V_HeAr(qp[IR]);
}

double numerator(hep::mc_point<double> const &x)
{
    double R = x.point()[IR] * (Rmax - Rmin) + Rmin;	
	double pR =  std::tan(M_PI * (x.point()[IPR] - 0.5));
    double Phi = x.point()[IPHI] * 2.0 * M_PI;
    double pPhi = std::tan(M_PI * (x.point()[IPPHI] - 0.5));
	double Theta = x.point()[ITHETA] * M_PI;
	double pTheta =  std::tan(M_PI * (x.point()[IPTHETA] - 0.5));

    double RJ = Rmax - Rmin;
	double PRJ = M_PI * (1.0 + pR * pR);
    double PHJ = 2.0 * M_PI;
    double PPHJ = M_PI * (1.0 + pPhi * pPhi);
	double THJ  = M_PI;
	double PTHJ = M_PI * (1.0 + pTheta * pTheta);

	double Jac = RJ * PRJ * THJ * PTHJ * PHJ * PPHJ;

    double qp[DIM] = {R, pR, Theta, pTheta, Phi, pPhi};
    double en = Hamiltonian(qp);
    
    double dip = dip_HeAr(R);
    double dipsq = dip*dip;
    // double dipsq = R;

    double result = Jac * dipsq * std::exp(-en/Temperature/Boltzmann_Hartree);
	return result;
}

double denominator(hep::mc_point<double> const& x)
{
    double R = x.point()[IR] * (Rmax - Rmin) + Rmin;	
	double pR =  std::tan(M_PI * (x.point()[IPR] - 0.5));
    double Phi = x.point()[IPHI] * 2.0 * M_PI;
    double pPhi = std::tan(M_PI * (x.point()[IPPHI] - 0.5));
	double Theta = x.point()[ITHETA] * M_PI;
	double pTheta =  std::tan(M_PI * (x.point()[IPTHETA] - 0.5));

    double RJ = Rmax - Rmin;
	double PRJ = M_PI * (1.0 + pR * pR);
    double PHJ = 2.0 * M_PI;
    double PPHJ = M_PI * (1.0 + pPhi * pPhi);
	double THJ  = M_PI;
	double PTHJ = M_PI * (1.0 + pTheta * pTheta);

	double Jac = RJ * PRJ * THJ * PTHJ * PHJ * PPHJ;

    double qp[DIM] = {R, pR, Theta, pTheta, Phi, pPhi};
    double en = Hamiltonian(qp);
    
    double result = Jac * std::exp(-en/Temperature/Boltzmann_Hartree);
	return result;
}

/*
double numerator_cnfg(hep::mc_point<double> const& x)
{
    double R = x.point()[IR] * (Rmax - Rmin) + Rmin;	
    double Phi = x.point()[IPHI] * 2.0 * M_PI;
	double Theta = x.point()[ITHETA] * M_PI;
    
    (void) Phi;
    (void) Theta;

    double RJ = Rmax - Rmin;
    double PHJ = 2.0 * M_PI;
	double THJ  = M_PI;
	double Jac = RJ * THJ * PHJ;

    double V = V_HeAr(R);
    
    double dip = dip_HeAr(R);
    double dipsq = dip*dip;
    // double dipsq = R;
    
    double result = Jac * dipsq * std::exp(-V/Temperature/Boltzmann_Hartree);
	return result;
}

double denominator_cnfg(hep::mc_point<double> const& x)
{
    double R = x.point()[IR] * (Rmax - Rmin) + Rmin;	
    double Phi = x.point()[IPHI] * 2.0 * M_PI;
	double Theta = x.point()[ITHETA] * M_PI;
    
    (void) Phi;
    (void) Theta;

    double RJ = Rmax - Rmin;
    double PHJ = 2.0 * M_PI;
	double THJ  = M_PI;
	double Jac = RJ * THJ * PHJ;

    double V = V_HeAr(R);
    
    double result = Jac * std::exp(-V/Temperature/Boltzmann_Hartree);
	return result;
}
*/

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    int world_size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
    hep::mpi_vegas_callback<double>(hep::mpi_vegas_verbose_callback<double>);
		
    auto numerator_results = hep::mpi_vegas(
        MPI_COMM_WORLD,
        hep::make_integrand<double>(numerator, DIM),
        std::vector<std::size_t>(7, 5.0e6)
    );

    auto denominator_results = hep::mpi_vegas(
        MPI_COMM_WORLD,
        hep::make_integrand<double>(denominator, DIM),
        std::vector<std::size_t>(7, 5.0e6)
    );
		
    auto numerator_result   = hep::accumulate<hep::weighted_with_variance>(numerator_results.begin() + 2, numerator_results.end());
    auto denominator_result = hep::accumulate<hep::weighted_with_variance>(denominator_results.begin() + 2, denominator_results.end());

    // double V = 4.0*M_PI/3.0 * Rmax*Rmax*Rmax; 
    // double den = std::pow(2.0*M_PI*mu*Boltzmann_Hartree*Temperature, 1.5) * V; 
    
    double dipsq_mean = numerator_result.value() / denominator_result.value();

    if (rank == 0) {
        std::cout << std::fixed << std::setprecision(5);
        std::cout << "numerator: " << numerator_result.value() << std::endl;
        std::cout << "denominator: " << denominator_result.value() << std::endl;
        // std::cout << "(analytic) denominator: " << den << std::endl;

        std::cout << std::scientific;
        std::cout << "<mu^2>: " << dipsq_mean << "\n";
        // std::cout << "V*<mu^2>: " << V*dipsq_mean << "\n";
    }

    MPI_Finalize();

    return 0;
}
