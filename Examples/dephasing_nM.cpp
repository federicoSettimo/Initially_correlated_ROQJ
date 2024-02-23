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

//double gamma (double t) {return cos(M_PI*t)/(1.+2.*t*t);} // gamma_x
//double gamma (double t) {return cos(.5*M_PI*t)/(1.+.5*t*t);} // gamma_y
//double gamma (double t) {return 1.;} // gamma_z
double gamma (double t) {return sin(2.*M_PI*t)/(1.+.5*t*t);} // gamma_0

Matrix2cd Gamma (double t) {return gamma(t)*sigma_z*sigma_z;}

Matrix2cd J (const Matrix2cd X, double t) {return gamma(t)*sigma_z*X*sigma_z;}

Matrix2cd L (const Matrix2cd X, double t) {return J(X,t) - .5*anticomm(Gamma(t),X);}

double xBloch (const Matrix2cd X) {return real((X*sigma_x).trace());}
double yBloch (const Matrix2cd X) {return real((X*sigma_y).trace());}

int main () {
    int Nstates = 100000, Npsi_p = Nstates, Npsibar_p = 0, Npsi_m = Nstates, Npsibar_m = 0;
    double tmax = 3., dt = .001;
    
    Vector2cd phi_p = (sqrt(3)-1.)/(2.*sqrt(3.-sqrt(3)))*(I-1.)*excited_state + 1./sqrt(3.-sqrt(3))*ground_state;
    Vector2cd barphi_p = (sqrt(3)-1.)/(2.*sqrt(3.-sqrt(3)))*(I-1.)*excited_state - 1./sqrt(3.-sqrt(3))*ground_state;
    Vector2cd phi_m = -(sqrt(3)+1.)/(2.*sqrt(3.+sqrt(3)))*(I-1.)*excited_state + 1./sqrt(3.+sqrt(3))*ground_state;
    Vector2cd barphi_m = -(sqrt(3)+1.)/(2.*sqrt(3.+sqrt(3)))*(I-1.)*excited_state - 1./sqrt(3.+sqrt(3))*ground_state;
    Vector2cd psi_p = phi_p, psibar_p = barphi_p;
    Vector2cd psi_m = phi_m, psibar_m = barphi_m;
    Matrix2cd X = .5*(id - sigma_x - sigma_y - sigma_z);
    string fileID = "_0.txt";

    ofstream out_ex, out_avg, out_params, out_gamma;
    out_ex.open("exact"+fileID);
    out_avg.open("avg"+fileID);
    out_params.open("params"+fileID);
    out_gamma.open("gamma"+fileID);
    out_params << tmax << endl << dt << endl;

    for (double t = 0.; t < tmax; t += dt) {
        out_gamma << gamma(t) << endl;

        // Exact solution
        out_ex << xBloch(X) << " " << yBloch(X) << endl;
        X += L(X,t)*dt;

        // Trajectories
        //Matrix2cd Xavg = ((double)Npsi_p*projector(psi_p) +(double) Npsibar_p*projector(psibar_p) - (double)Npsi_m*projector(psi_m) - (double)Npsibar_m*projector(psibar_m))/(double)(2*Nstates);
        //double mu_p = .5*sqrt(3.-sqrt(3))*(1.+sqrt(3)), mu_m = .5*sqrt(3.+sqrt(3))*(1.-sqrt(3));
        double mu_p = .5*(1.+sqrt(3)), mu_m = .5*(1.-sqrt(3));
        Matrix2cd Xavg = (mu_p*((double)Npsi_p*projector(psi_p) +(double) Npsibar_p*projector(psibar_p)) + mu_m*((double)Npsi_m*projector(psi_m) + (double)Npsibar_m*projector(psibar_m)))/(double)Nstates;
        out_avg << xBloch(Xavg) << " " << yBloch(Xavg) << endl;

        // First with +
        int Npsi_old = Npsi_p, Npsibar_old = Npsibar_p;
        double p = (double)Npsi_p/(double)Nstates;
        // States in psi
        for (int i = 0; i < Npsi_old; ++i) {
            double p_jump = gamma(t)*dt;
            if (gamma(t) < 0.) {
                if (p >= .5)
                    p_jump = 0.;
                else p_jump = -dt*abs(gamma(t)*(1.-2.*p)/p);
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
                if (p >= .5)
                    p_jump = dt*abs(gamma(t)*(1.-2.*p)/(1.-p));
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
        p = (double)Npsi_m/(double)Nstates;
        // States in psi
        for (int i = 0; i < Npsi_old; ++i) {
            double p_jump = gamma(t)*dt;
            if (gamma(t) < 0.) {
                if (p >= .5)
                    p_jump = 0.;
                else p_jump = -dt*gamma(t)*(1.-2.*p)/p;
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
                if (p >= .5)
                    p_jump = dt*gamma(t)*(1.-2.*p)/(1.-p);
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