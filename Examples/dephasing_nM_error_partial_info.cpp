// Dephasing partial info, now computing the error depending on the numer of realizations.
// Using the map acting on the state for simplicity
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/MatrixFunctions>
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

double gamma (double t) {return cos(M_PI*t)/(1.+2.*t*t);} // gamma_x
//double gamma (double t) {return cos(.5*M_PI*t)/(1.+.5*t*t);} // gamma_y
//double gamma (double t) {return 1.;} // gamma_z
//double gamma (double t) {return sin(2.*M_PI*t)/(1.+.5*t*t);} // gamma_0

Matrix2cd Gamma (double t) {return gamma(t)*sigma_z*sigma_z;}

Matrix2cd J (const Matrix2cd X, double t) {return gamma(t)*sigma_z*X*sigma_z;}

Matrix2cd L (const Matrix2cd X, double t) {return J(X,t) - .5*anticomm(Gamma(t),X);}

double Fidelity (const Matrix2cd rho, const Matrix2cd sigma) {
    Matrix2cd sqrt_rho = rho.sqrt(), prod = sqrt_rho*sigma*sqrt_rho, sqrt_prod = prod.sqrt();
    return pow(real(sqrt_prod.trace()),2);
}

double TD (const Matrix2cd rho, const Matrix2cd sigma) {
    Matrix2cd X = rho-sigma;
    return .5*real((X.adjoint()*X).trace());
}

int main () {
    int Nstates = 10000;
    double tmax = 1.5, dt = .005;

    Vector2cd psi_p = plus_state, psibar_p = minus_state;
    Matrix2cd rho = projector(plus_state);

    ofstream out_gamma, out_err;
    out_gamma.open("gamma_error_partial_info.txt");
    out_err.open("error_partial_info.txt");
    out_gamma << tmax << endl << dt << endl;

    vector<bool> vec_psi_p; // If true = the state is in psi, else it's in psibar
    for (int i = 0; i < Nstates; ++i) {
        vec_psi_p.push_back(true);
    }

    for (double t = 0.; t < tmax; t += dt) {
        out_gamma << gamma(t) << endl;

        // First with +
        int Npsi_p = 0, Npsibar_p = 0;
        for (int i = 0; i < Nstates; ++i) {
            double p_jump = gamma(t)*dt;
            if (vec_psi_p[i]) {
                Npsi_p++;
                double p = (double)Npsi_p/((double)i+1.); // fraction of states in psi so far - only this info needed
                if (gamma(t) < 0.) {
                    if (p == 0.)
                        p_jump = 1.;
                    else if (p >= .5)
                        p_jump = 0.;
                    else p_jump = -dt*abs(gamma(t)*(1.-2.*p)/p);
                }
            }
            else {
                Npsibar_p++;
                double p = (double)Npsi_p/((double)i+1.); // fraction of states in psi so far
                if (gamma(t) < 0.) {
                    if (p == 1.)
                        p_jump = 1.;
                    if (p >= .5)
                        p_jump = dt*abs(gamma(t)*(1.-2.*p)/(1.-p));
                    else p_jump = 0.;;
                }
            }
            double z = (double)rand()/(double)RAND_MAX;
            if (z <= p_jump) 
                vec_psi_p[i] = !vec_psi_p[i];
            Matrix2cd rho_avg = ((double)Npsi_p*projector(psi_p) + (double)Npsibar_p*projector(psibar_p))/((double)i+1.);
            if (i >= 49 && gamma(t) < 0.)
                //out_err << i+1 << " " << t << " " << Fidelity(rho,rho_avg) << endl;
                out_err << i+1 << " " << t << " " << TD(rho,rho_avg) << endl;
        }
        // Exact solution
        rho += L(rho,t)*dt;
    }
}