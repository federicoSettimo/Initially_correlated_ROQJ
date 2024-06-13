#include "../roqj_pop.h"

using namespace std;
using namespace Eigen;

static Eigen::VectorXcd plus_y = (excited_state + I*ground_state)/sqrt(2.), minus_y = (excited_state - I*ground_state)/sqrt(2.);

// Parameters...
int N_ensemble = 10000, Ncopies = 3, dimH = 2, Ntraj = 10;
double tmin = 0., tmax = M_PI, dt = 0.001;
bool printTraj = true;

MatrixXcd H (double t) {
  return (1.-cos(t))*sigma_z;
}

Matrix2cd A (double t) {return cos(t)*sigma_x-sin(t)*sigma_y;}
Matrix2cd Atilde (double t) {return sin(t)*sigma_x-(1.-cos(t))*sigma_y;}

MatrixXcd J (const MatrixXcd &rho, double t) {
  return A(t)*rho*Atilde(t) + Atilde(t)*rho*A(t);
}

MatrixXcd Gamma (double t) {
  return 2.*id*sin(t);
}

VectorXcd Phi (const VectorXcd &psi, double t, bool jumped) {
    if (jumped)
      return 0.*ground_state;
    complex<double> phi1, phi0, psi1 = psi(0), psi0 = psi(1);
    phi1 = -2.*psi0*conj(psi0)*sin(t)/psi1;
    phi0 = ((-cos(t)+I*sin(t))*(psi0*conj(phi1)*(cos(t) + I*sin(t)) + 4.*psi1*conj(psi0)*(-I + I*cos(t) + sin(t))))/psi1;
    return phi1*excited_state + phi0*ground_state;
}

double observable (const MatrixXcd &rho) {return real((rho*sigma_x).trace());}

double gamma_p (double t) {return sin(t)+2.*sin(.5*t);}
double gamma_m (double t) {return sin(t)-2.*sin(.5*t);}


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
    qubit_roqj_pop jump(N_ensemble, tmin, tmax, dt, Ncopies, printTraj, Ntraj, true);

    string index = "_p.txt";
    //string index = "_m.txt";
    if (index == "_p.txt")
        jump.set_initial_state_R(phi_p);
    else
        jump.set_initial_state_R(phi_m);
    jump.run();
    jump.get_observable("average"+index);
    jump.get_error_observable("error"+index);
    jump.get_trajectories ("trajectories"+index);

    ofstream out_g;
    out_g.open("gammas.txt");
    for (double t = 0.; t < tmax; t += dt)
        out_g << gamma_m(t) << " " << gamma_p(t) << endl;

  return 0;
}