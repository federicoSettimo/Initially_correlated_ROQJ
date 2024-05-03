// Zero-discord initial state rho = p |0><0| \otimes |1><1| + (1-p) |1><1| \otimes |0><0|
#include "../roqj_pop.h"

using namespace std;
using namespace Eigen;

static Eigen::VectorXcd plus_y = (excited_state + I*ground_state)/sqrt(2.), minus_y = (excited_state - I*ground_state)/sqrt(2.);

// Parameters...
double tmin = 0., tmax = 25, dt = 0.001, Delta = .1, g = .1, p = .8, cutoff_rates = 2.;
int N_ensemble = 10000, Ncopies = 3, dimH = 2, Ntraj = 10;
bool printTraj = false;

double sgn (double x) {return x >= 0. ? 1. : -1.;}

// Functions of JC, etc
double Omega (int n) {return sqrt(Delta*Delta + 4.*n*g*g);}

complex<double> c (int n, double t) {
    return exp(.5*I*Delta*t)*(cos(.5*Omega(n)*t) - I*Delta*sin(.5*Omega(n)*t)/Omega(n));
}
complex<double> cprime (int n, double t) {
    return .5*I*Delta*c(n,t) - .5*exp(.5*I*Delta*t)*Omega(n)*sin(.5*Omega(n)*t) - .5*exp(.5*I*Delta*t)*I*Delta*cos(.5*Omega(n)*t);
}

complex<double> alpha (int n, double t) {
    return conj(c(n,t))*c(n,t);
}
complex<double> alphaprime (int n, double t) {
    return conj(cprime(n,t))*c(n,t) + conj(c(n,t))*cprime(n,t);
}

complex<double> beta (int n, double t) {
    return conj(c(n+1,t))*c(n+1,t);
}
complex<double> betaprime (int n, double t) {
    return conj(cprime(n+1,t))*c(n+1,t) + conj(c(n+1,t))*cprime(n+1,t);
}

complex<double> gamma (int n, double t) {
    return c(n,t)*c(n+1,t);
}
complex<double> gammaprime (int n, double t) {
    return cprime(n,t)*c(n+1,t) + c(n,t)*cprime(n+1,t);
}

double gamma_p_gen (int n, double t) {return real(((alpha(n,t) - 1.)*betaprime(n,t) - beta(n,t)*alphaprime(n,t))/(alpha(n,t) + beta(n,t) - 1.));}
double gamma_m_gen (int n, double t) {return real(((beta(n,t) - 1.)*alphaprime(n,t) - alpha(n,t)*betaprime(n,t))/(alpha(n,t) + beta(n,t) - 1.));}
double gamma_z_gen (int n, double t) {return real((alphaprime(n,t) + betaprime(n,t))/(4.*(alpha(n,t) + beta(n,t) - 1.))) - .5*real(gammaprime(n,t)/(gamma(n,t)));}

// Functions for rho_E^{0,1}
// They are JC but with a cutoff such that |gamma_\alpha| <= cutoff
// For rho_E^0
double gamma_p_0 (double t) {
    int n = 1;
    if (alpha(n,t) + beta(n,t) - 1. == 0.)
        return cutoff_rates;
    double g = gamma_p_gen(n,t);
    if (g >= 0.)
        return min(g, cutoff_rates);
    return max(g, -cutoff_rates);
}
double gamma_m_0 (double t) {
    int n = 1;
    if (alpha(n,t) + beta(n,t) - 1. == 0.)
        return cutoff_rates;
    double g = gamma_m_gen(n,t);
    if (g >= 0.)
        return min(g, cutoff_rates);
    return max(g, -cutoff_rates);
}
double gamma_z_0 (double t) {
    int n = 1;
    if (alpha(n,t) + beta(n,t) - 1. == 0. || gamma(n,t) == 0.)
        return cutoff_rates;
    double g = gamma_z_gen(n,t);
    if (g >= 0.)
        return min(g, cutoff_rates);
    return max(g, -cutoff_rates);
}

// And for \rho_E^1
double gamma_p_1 (double t) {
    int n = 0;
    if (alpha(n,t) + beta(n,t) - 1. == 0.)
        return cutoff_rates;
    double g = gamma_p_gen(n,t);
    if (g >= 0.)
        return min(g, cutoff_rates);
    return max(g, -cutoff_rates);
}
double gamma_m_1 (double t) {
    int n = 0;
    if (alpha(n,t) + beta(n,t) - 1. == 0.)
        return cutoff_rates;
    double g = gamma_m_gen(n,t);
    if (g >= 0.)
        return min(g, cutoff_rates);
    return max(g, -cutoff_rates);
}
double gamma_z_1 (double t) {
    int n = 0;
    if (alpha(n,t) + beta(n,t) - 1. == 0. || gamma(n,t) == 0.)
        return cutoff_rates;
    double g = gamma_z_gen(n,t);
    if (g >= 0.)
        return min(g, cutoff_rates);
    return max(g, -cutoff_rates);
}

// Maps \Phi^{0,x,y}:
double gamma_p (double t) {return p*gamma_p_0(t) + (1.-p)*gamma_p_1(t);}
double gamma_m (double t) {return p*gamma_m_0(t) + (1.-p)*gamma_m_1(t);}
double gamma_z (double t) {return p*gamma_z_0(t) + (1.-p)*gamma_z_1(t);}
// Map \Phi^z
/*double gamma_p (double t) {return gamma_p_1(t);}
double gamma_m (double t) {return gamma_m_1(t);}
double gamma_z (double t) {return gamma_z_1(t);}*/
// Dummy maps to avoid negativity since t = 0
/*double gamma_p (double t) {return 1.;}
double gamma_m (double t) {return .5;}
double gamma_z (double t) {return cos(t);}*/

MatrixXcd H (double t) {
  return 0*sigma_z;
}

MatrixXcd J (const MatrixXcd &rho, double t) {
  return gamma_p(t)*sigma_p*rho*sigma_m + gamma_m(t)*sigma_m*rho*sigma_p + gamma_z(t)*sigma_z*rho*sigma_z;
}

MatrixXcd Gamma (double t) {
  return gamma_p(t)*sigma_m*sigma_p + gamma_m(t)*sigma_p*sigma_m + gamma_z(t)*id;
}

VectorXcd Phi (const VectorXcd &psi, double t, bool jumped) {
    double gz = gamma_z(t), gp = gamma_p(t), gm = gamma_m(t);
    complex<double> ps0 = psi(1), ps1 = psi(0), ph0, ph1;

    if (jumped) { // |1>
        if (psi(1) == 0.)
            return -gamma_z(t) * excited_state;
        return -gamma_z(t) * ground_state;
    }

    if (abs(ps1) == 0.)
        return -gz * ground_state;
    else if (abs(ps1) == 1.)
        return -gz * excited_state;


    ph1 = (2.*gz*ps1*conj(ps1) + (gm*ps1*conj(ps1)*ps1*conj(ps1))/(ps0*conj(ps0)) - gp*ps0*conj(ps0))/(2.*ps1);
    ph0 = 2.*gz*ps0 - ps0*conj(ph1)/ps1;
    return ph0*ground_state + ph1*excited_state;
}

double observable (const MatrixXcd &rho) {return real((rho*sigma_x).trace());}

MatrixXcd L_t (const MatrixXcd &rho, double t) {
    return -I*comm(H(t),rho) + J(rho, t) - .5*anticomm(Gamma(t), rho);
}

int main () {
    srand(time(NULL));
    // Defining the initial operators Q_i (use only one and comment out the others)
    // For Q_0:
    /*Vector2cd phi_p = (sqrt(3)-1.)/(2.*sqrt(3.-sqrt(3)))*(I-1.)*excited_state + 1./sqrt(3.-sqrt(3))*ground_state;
    Vector2cd phi_m = -(sqrt(3)+1.)/(2.*sqrt(3.+sqrt(3)))*(I-1.)*excited_state + 1./sqrt(3.+sqrt(3))*ground_state;
    Matrix2cd X0 = .5*(id - sigma_x - sigma_y - sigma_z);
    double mu_p = .5*(1.+sqrt(3)), mu_m = .5*(1.-sqrt(3));*/

    // For Q_x:
    Vector2cd phi_p = plus_state, phi_m = minus_state;
    Matrix2cd X0 = .5*sigma_x;
    double mu_p = .5, mu_m = -.5;

    // For Q_y:
    /*Vector2cd phi_p = plus_y, phi_m = minus_y;
    Matrix2cd X0 = .5*sigma_y;
    double mu_p = .5, mu_m = -.5;*/

    // For Q_z:
    /*Vector2cd phi_p = excited_state, phi_m = ground_state;
    Matrix2cd X0 = .5*sigma_z;
    double mu_p = .5, mu_m = -.5;*/

    // Now running the simulations
    qubit_roqj_pop jump(N_ensemble, tmin, tmax, dt, Ncopies, printTraj, Ntraj, false);

    jump.set_initial_state(phi_p);
    //jump.set_initial_state(phi_m);
    jump.run();
    jump.get_observable("average_p.txt");
    jump.get_error_observable("error_p.txt");

    ofstream out_g;
    out_g.open("gammas.txt");
    for (double t = 0.; t < tmax; t += dt)
        out_g << gamma_m(t) << " " << gamma_p(t) << " " << gamma_z(t) << endl;


  return 0;
}