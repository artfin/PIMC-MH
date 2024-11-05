#ifndef HTOCM
#define HTOCM 2.1947463136320e5  // 1 Hartree in cm-1
#endif

// Botlzmann*(T = 300K) = 208 cm-1
// => we need to consider temperatures somewhat below 100K to avoid free states
double morse_cm(double r) {
    double De = 3000.0; // cm-1
    double a = 0.5;
    // double a_ang = a/ALU/1.0e10;
    double re = 6.5; // Bohr
    // double reang = 6.5 * ALU * 1.0e10;
    // double mu = 0.17; 
    // double MU = mu * HTOCM; 
    // double nu0 = a/(2.0*M_PI) * sqrt(2.0 * De / mu);
    // double omega0 = 2.0*M_PI*nu0;

    return De*(1.0 - exp(-a*(r-re)))*(1.0 - exp(-a*(r-re))) - De;
}

double morse(double r) 
// r: Bohr -> pot: Hartree
{
    return morse_cm(r) / HTOCM;
}

