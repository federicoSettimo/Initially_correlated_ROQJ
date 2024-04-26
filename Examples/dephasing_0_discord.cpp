// Spin-boson dephasing, initial state 0 discord \rho(0) = p |0><0| \rho_E^0 + (1-p) |1><1| \rho_E^1
// \rho_E^i thermal states
// Two environmental modes with frequencies omega_i and coupling strengths g_1
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <complex>
#include <vector>
#include <cstring>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <cmath>

using namespace std;
using namespace Eigen;

static complex<double> I(0,1), one(1,0), zero(0,0);
static Eigen::Matrix2cd sigma_x {{0,1},{1,0}};
static Eigen::Matrix2cd sigma_y {{0,-I},{I,0}};
static Eigen::Matrix2cd sigma_z {{1,0},{0,-1}};
static Eigen::Matrix2cd sigma_p {{0,1},{0,0}};
static Eigen::Matrix2cd sigma_m {{0,0},{1,0}};
static Eigen::Matrix2cd id {{1,0},{0,1}};

static Eigen::VectorXcd ground_state {{0.,1.}};
static Eigen::VectorXcd excited_state {{1.,0.}};
static Eigen::VectorXcd plus_state {{1./sqrt(2.),1./sqrt(2.)}};
static Eigen::VectorXcd minus_state {{1./sqrt(2.),-1./sqrt(2.)}};
static Eigen::VectorXcd plus_y = (excited_state + I*ground_state)/sqrt(2.), minus_y = (excited_state - I*ground_state)/sqrt(2.);

Matrix2cd comm (const Matrix2cd &A, const Matrix2cd &B) {return A*B-B*A;}

Matrix2cd anticomm (const Matrix2cd &A, const Matrix2cd &B) {return A*B+B*A;}

Matrix2cd projector (const Vector2cd &psi) {return psi*psi.adjoint();}

double coth (double x) {return cosh(x)/sinh(x);}

// Parameters defining the gamma_i functions: omega = frequencies, g = cupling strengths, beta = inverse temperature for \rho_E^1 (\rho_E^0 at T = 0)
//double omega_1 = 1., omega_2 = 1.2, g_1 = 1., g_2 = .8, beta = .4;
double omega_1 = 1.8, omega_2 = 2.3, g_1 = .5, g_2 = .4, beta = .4;

double gamma_0 (double t) { // 0 temperature environmental state
    return g_1*g_1*sin(omega_1*t)/omega_1 + g_2*g_2*sin(omega_2*t)/omega_2;
}
double gamma_1 (double t) { // \rho_E^1 as environmental state (at T = 1/\beta)
    return g_1*g_1*coth(.5*beta*omega_1)*sin(omega_1*t)/omega_1 + g_2*g_2*coth(.5*beta*omega_2)*sin(omega_2*t)/omega_2;
}

// \rho(0) = p |0><0| \rho_E^0 + (1-p) |1><1| \rho_E^1
double p = .8;
// \tilde\gamma_0 = \tilde\gamma_x = \tilde\gamma_y
// Defining the maps \Phi^{0,x,y}_t
//double gamma (double t) {return p*gamma_0(t) + (1.-p)*gamma_1(t);}
// \tilde\gamma_z, defining \phi^z_t
double gamma (double t) {return gamma_1(t);}




Matrix2cd Gamma (double t) {return gamma(t)*sigma_z*sigma_z;}

Matrix2cd J (const Matrix2cd X, double t) {return gamma(t)*sigma_z*X*sigma_z;}

Matrix2cd L (const Matrix2cd X, double t) {return J(X,t) - .5*anticomm(Gamma(t),X);}

double xBloch (const Matrix2cd X) {return real((X*sigma_x).trace());}
double yBloch (const Matrix2cd X) {return real((X*sigma_y).trace());}

int main () {
    int Nstates = 10000, Npsi_p = Nstates, Npsibar_p = 0, Npsi_m = Nstates, Npsibar_m = 0;
    double tmax = 10., dt = .001;
    
    // Defining the initial operators Q_i (use only one and comment out the others)
    // For Q_0:
    /*Vector2cd phi_p = (sqrt(3)-1.)/(2.*sqrt(3.-sqrt(3)))*(I-1.)*excited_state + 1./sqrt(3.-sqrt(3))*ground_state;
    Vector2cd barphi_p = (sqrt(3)-1.)/(2.*sqrt(3.-sqrt(3)))*(I-1.)*excited_state - 1./sqrt(3.-sqrt(3))*ground_state;
    Vector2cd phi_m = -(sqrt(3)+1.)/(2.*sqrt(3.+sqrt(3)))*(I-1.)*excited_state + 1./sqrt(3.+sqrt(3))*ground_state;
    Vector2cd barphi_m = -(sqrt(3)+1.)/(2.*sqrt(3.+sqrt(3)))*(I-1.)*excited_state - 1./sqrt(3.+sqrt(3))*ground_state;
    Matrix2cd X0 = .5*(id - sigma_x - sigma_y - sigma_z);
    double mu_p = .5*(1.+sqrt(3)), mu_m = .5*(1.-sqrt(3));*/

    // For Q_x:
    /*Vector2cd phi_p = plus_state, phi_m = minus_state;
    Vector2cd barphi_p = minus_state, barphi_m = plus_state;
    Matrix2cd X0 = .5*sigma_x;
    double mu_p = .5, mu_m = -.5;*/

    // For Q_y:
    Vector2cd phi_p = plus_y, phi_m = minus_y;
    Vector2cd barphi_p = minus_y, barphi_m = plus_y;
    Matrix2cd X0 = .5*sigma_y;
    double mu_p = .5, mu_m = -.5;

    // fileID = _mapIndex_QIndex.txt --- _0_x.txt for \Phi^0_t(Q_x)
    string fileID = "_z_y.txt";

    ofstream out_ex, out_avg, out_params, out_gamma, out_gamma_0_original, out_gamma_1_original;
    out_ex.open("exact"+fileID);
    out_avg.open("avg"+fileID);
    out_params.open("params"+fileID);
    out_gamma.open("gamma"+fileID);
    out_gamma_0_original.open("gamma_0.txt");
    out_gamma_1_original.open("gamma_1.txt");
    out_params << omega_1 << endl << g_1 << endl << omega_2 << endl << g_2 << endl << beta << endl << p << endl << tmax << endl << dt << endl;

    Vector2cd psi_p = phi_p, psibar_p = barphi_p;
    Vector2cd psi_m = phi_m, psibar_m = barphi_m;

    for (double t = 0.; t < tmax; t += dt) {
        out_gamma << gamma(t) << endl;
        out_gamma_0_original << gamma_0(t) << endl;
        out_gamma_1_original << gamma_1(t) << endl;

        out_ex << xBloch(X0) << " " << yBloch(X0) << endl;
        X0 += L(X0,t)*dt;

        // Trajectories
        Matrix2cd Xavg = (mu_p*((double)Npsi_p*projector(psi_p) +(double) Npsibar_p*projector(psibar_p)) + mu_m*((double)Npsi_m*projector(psi_m) + (double)Npsibar_m*projector(psibar_m)))/(double)Nstates;
        out_avg << xBloch(Xavg) << " " << yBloch(Xavg) << endl;

        // First with +
        int Npsi_old = Npsi_p, Npsibar_old = Npsibar_p;
        double pp = (double)Npsi_p/(double)Nstates;
        // States in psi
        for (int i = 0; i < Npsi_old; ++i) {
            double p_jump = gamma(t)*dt;
            if (gamma(t) < 0.) {
                if (pp >= .5)
                    p_jump = 0.;
                else p_jump = -dt*abs(gamma(t)*(1.-2.*pp)/pp);
            }
            double z = (double)rand()/(double)RAND_MAX;
            if (z <= p_jump) {
                Npsi_p--;
                Npsibar_p++;
            }
        }
        // States in psibar
        for (int i = 0; i < Npsibar_old; ++i) {
            double p_jump = gamma(t)*dt;
            if (gamma(t) < 0.) {
                if (pp >= .5)
                    p_jump = dt*abs(gamma(t)*(1.-2.*pp)/(1.-pp));
                else p_jump = 0.;
            }
            double z = (double)rand()/(double)RAND_MAX;
            if (z <= p_jump) {
                Npsibar_p--;
                Npsi_p++;
            }
        }

        // Now with -
        Npsi_old = Npsi_m; Npsibar_old = Npsibar_m;
        pp = (double)Npsi_m/(double)Nstates;
        // States in psi
        for (int i = 0; i < Npsi_old; ++i) {
            double p_jump = gamma(t)*dt;
            if (gamma(t) < 0.) {
                if (pp >= .5)
                    p_jump = 0.;
                else p_jump = -dt*gamma(t)*(1.-2.*pp)/pp;
            }
            double z = (double)rand()/(double)RAND_MAX;
            if (z <= p_jump) {
                Npsi_m--;
                Npsibar_m++;
            }
        }
        // States in psibar
        for (int i = 0; i < Npsibar_old; ++i) {
            double p_jump = gamma(t)*dt;
            if (gamma(t) < 0.) {
                if (pp >= .5)
                    p_jump = dt*gamma(t)*(1.-2.*pp)/(1.-pp);
                else p_jump = 0.;;
            }
            double z = (double)rand()/(double)RAND_MAX;
            if (z <= p_jump) {
                Npsibar_m--;
                Npsi_m++;
            }
        }
    }
}