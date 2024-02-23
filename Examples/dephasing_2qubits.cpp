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

Vector4cd ket_zero{{0.,0.,0.,1.}}, ket_one{{0.,0.,1.,0.}}, ket_two {{0.,1.,0.,0.}}, ket_three{{1.,0.,0.,0.}};
DiagonalMatrix<complex<double>, 4> id (1., 1., 1., 1.);
complex<double> I(0,1), one(1,0), zero(0,0);

Matrix4cd comm (const Matrix4cd &A, const Matrix4cd &B) {return A*B-B*A;}
Matrix4cd anticomm (const Matrix4cd &A, const Matrix4cd &B) {return A*B+B*A;}
Matrix4cd projector (const Vector4cd &psi) {return psi*psi.adjoint();}

Vector4cd ket_i (int i) {
    if (i == 0) return ket_zero;
    else if (i == 1) return ket_one;
    else if (i == 2) return ket_two;
    else if (i == 3) return ket_three;
    else exit(EXIT_FAILURE);
}

Vector4cd ket_plus (int i, int j, bool isX) { // plus eigenstate of the sigma_x or sigma_y matrix in the i,j subspace
    if (isX) return (ket_i(i) + ket_i(j))/sqrt(2.);
    return (ket_i(i) + I*ket_i(j))/sqrt(2.);
}

Vector4cd ket_minus (int i, int j, bool isX) { // minus eigenstate of the sigma_x or sigma_y matrix in the i,j subspace
    if (isX) return (ket_i(i) - ket_i(j))/sqrt(2.);
    return (ket_i(i) - I*ket_i(j))/sqrt(2.);
}

Matrix4cd sigma (int i, int j) {return projector(ket_i(i)) - projector(ket_i(j));}

double gamma (double t, int i, int j, bool isX) { // The gamma functions are different if considering the x or y map, however the operators are the same
    /*if (isX && i == 1 && j == 0)
        return cos(M_PI*t)/(1.+2.*t*t);
    else if (isX && i == 3 && j == 0)
        return sin(2.*M_PI*t)/(1.+.5*t*t);
    else if (isX && i == 2 && j == 0)
        return cos(.5*M_PI*t)/(1.+.5*t*t);*/
    return .1*(double)(4.*i + j);
}

Matrix4cd Gamma (double t, int i, int j, bool isX) {return gamma(t,i,j,isX)*sigma(i,j).adjoint()*sigma(i,j);}

Matrix4cd J (const Matrix4cd X, double t, int i, int j, bool isX) {return gamma(t,i,j,isX)*sigma(i,j)*X*sigma(i,j).adjoint();}

Matrix4cd L (const Matrix4cd X, double t, int i, int j, bool isX) {return J(X,t,i,j,isX) - .5*anticomm(Gamma(t,i,j,isX),X);}


// For 0
double gamma0 (double t, int i, int j) { // Some gammas must be the same - see notes
    return 1.;
}

Matrix4cd Gamma0 (double t) {
    Matrix4cd g;
    for (int i = 1; i < 4; ++i) {
        for (int j = 0; j < i; ++j) {
            g += gamma0(t,i,j)*sigma(i,j).adjoint()*sigma(i,j);
        }
    }
    return g;
}

Matrix4cd J0 (const Matrix4cd X, double t) {
    Matrix4cd JJ;
    for (int i = 1; i < 4; ++i) {
        for (int j = 0; j < i; ++j) {
            JJ += gamma0(t,i,j)*sigma(i,j)*X*sigma(i,j).adjoint();
        }
    }
    return JJ;
}

Matrix4cd L0 (const Matrix4cd X, double t) {return J0(X,t) - .5*anticomm(Gamma0(t),X);}

double ReLocalCoherence (const Matrix4cd X) {return real(ket_two.dot(X*ket_zero)) + real(ket_one.dot(X*ket_zero));} // <2|X|0> + <1|X|0>: the two terms that I would get if one of the two is |0>
double ReNonLocalCoherence (const Matrix4cd X) {return real(ket_three.dot(X*ket_zero));} // The coherent part of the maximally entangled state

int main () {
    int Nstates = 10000;
    double tmax = 3., dt = .001;

    // First: do it for Phi^0(Q_0)
    cout << "Unraveling Phi^0...\n";
    Matrix4cd X = .5*id;
    for (int i = 1; i < 4; ++i) {
        for (int j = 0; j < i; ++j) {
            for (int k = 0; k < 2; ++k) {
                bool isX = k == 0;
                X -= .5*(projector(ket_plus(i,j,isX)) - projector(ket_minus(i,j,isX)));
            }
        }
    }
    string fileID = "_0.txt";
    Vector4cd phi1u = (-0.92388 - 0.382683*I)*ket_three + (-0.707107 + 0.707107*I)*ket_two + (0.382683 + 0.92388*I)*ket_one + ket_zero,
        phi2u = (0.382683 - 0.92388*I)*ket_three + (0.707107 - 0.707107*I)*ket_two + (0.92388 - 0.382683*I)*ket_one + ket_zero,
        phi3u = (-0.382683 + 0.92388*I)*ket_three + (0.707107 - 0.707107*I)*ket_two + (-0.92388 + 0.382683*I)*ket_one + ket_zero,
        phi4u = (0.92388 + 0.382683*I)*ket_three + (-0.707107 + 0.707107*I)*ket_two + (-0.382683 - 0.92388*I)*ket_one + ket_zero;
    Vector4cd phi1 = phi1u.normalized(), phi2 = phi2u.normalized(), phi3 = phi3u.normalized(), phi4 = phi4u.normalized();
    //double mu1 = 1.7483, mu2 = -1.51367, mu3 = 1.09946, mu4 = 0.665911;
    double mu1 = 1.7483*pow(phi1u.norm(),2), mu2 = -1.51367*pow(phi2u.norm(),2), mu3 = 1.09946*pow(phi3u.norm(),2), mu4 = 0.665911*pow(phi4u.norm(),2);
    vector<Vector4cd> psis;
    vector<double> mu;
    for (int i = 0; i < Nstates/4; ++i) {
        psis.push_back(phi1);
        psis.push_back(phi2);
        psis.push_back(phi3);
        psis.push_back(phi4);
        mu.push_back(mu1);
        mu.push_back(mu2);
        mu.push_back(mu3);
        mu.push_back(mu4);
    }
    ofstream out_ex, out_avg, out_params, out_gamma;
    out_ex.open("exact"+fileID);
    out_avg.open("avg"+fileID);
    out_params.open("params"+fileID);
    out_gamma.open("gamma"+fileID);
    out_params << tmax << endl << dt << endl;
    for (double t = 0.; t < tmax; t += dt) {
        for (int i = 1; i < 4; ++i) {
            for (int j = 0; j < i; ++j) {
                out_gamma << gamma0(t,i,j) << " ";
            }
        }
        out_gamma << endl;

        // Exact solution
        out_ex << ReLocalCoherence(X) << " " << ReNonLocalCoherence(X) << endl;
        X += L0(X,t)*dt;

        Matrix4cd Xavg, K = -.5*I*Gamma0(t);
        for (int i = 0; i < Nstates; ++i) {
            Xavg += mu[i]*projector(psis[i])/(double)(Nstates);

            Matrix4cd R = J0(projector(psis[i]),t);
            double z = (double)rand()/(double)RAND_MAX, pj = real(R.trace())*dt;
            if (z <= pj) {
                ComplexEigenSolver<Matrix4cd> eigs;
                eigs.compute(R);
                Vector4cd eigval = eigs.eigenvalues();
                Matrix4cd eigvec = eigs.eigenvectors();
                if (z <= real(eigval[0])*dt)
                    psis[i] = eigvec.col(0);
                else if (z <= real(eigval[0] + eigval[1])*dt)
                    psis[i] = eigvec.col(1);
                else if (z <= real(eigval[0] + eigval[1] + eigval[2])*dt)
                    psis[i] = eigvec.col(2);
                else psis[i] = eigvec.col(3);
            }
            else 
                psis[i] -= I*K*psis[i]*dt;
            psis[i].normalize();
        }
        out_avg << ReLocalCoherence(Xavg) << " " << ReNonLocalCoherence(Xavg) << endl;
    }



    // Now all other maps: they are dephasing in the subspace |i>, |j>
    for (int i = 1; i < 4; ++i) {
        for (int j = 0; j < i; ++j) {
            for (int k = 0; k < 2; ++k) { // I have to do it for both x and y direction (diff Pauli matrices)
                bool isX = k == 0;
                Vector4cd psi_p = ket_plus(i,j,isX), psibar_p = ket_minus(i,j,isX);
                Vector4cd psi_m = psibar_p, psibar_m = psi_p;
                Matrix4cd X = .5*(projector(psi_p) - projector(psi_m));
                string fileID = "_"+to_string(i)+"_"+to_string(j);
                if (isX) fileID += "_x.txt";
                else fileID += "_y.txt";

                ofstream out_ex, out_avg, out_params, out_gamma;
                out_ex.open("exact"+fileID);
                out_avg.open("avg"+fileID);
                out_params.open("params"+fileID);
                out_gamma.open("gamma"+fileID);
                out_params << tmax << endl << dt << endl;
                if (isX) cout << "Unraveling Phi^{" << i << "," << j << "}...\n";

                int  Npsi_p = Nstates/2, Npsibar_p = 0, Npsi_m = Nstates/2, Npsibar_m = 0;

                // Now cycle on time for all maps...
                for (double t = 0.; t < tmax; t += dt) {
                    out_gamma << gamma(t,i,j,isX) << endl;

                    // Exact solution
                    out_ex << ReLocalCoherence(X) << " " << ReNonLocalCoherence(X) << endl;
                    X += L(X,t,i,j,isX)*dt;

                    // Trajectories
                    Matrix4cd Xavg = ((double)Npsi_p*projector(psi_p) +(double) Npsibar_p*projector(psibar_p) - (double)Npsi_m*projector(psi_m) - (double)Npsibar_m*projector(psibar_m))/(double)(Nstates);
                    out_avg << ReLocalCoherence(Xavg) << " " << ReNonLocalCoherence(Xavg) << endl;

                    // First with +
                    int Npsi_old = Npsi_p, Npsibar_old = Npsibar_p;
                    double p = (double)Npsi_p/(double)Nstates;
                    // States in psi
                    for (int ii = 0; ii < Npsi_old; ++ii) {
                        double p_jump = gamma(t,i,j,isX)*dt;
                        if (gamma(t,i,j,isX) < 0.) {
                            if (p >= .5)
                                p_jump = 0.;
                            else p_jump = dt*abs(gamma(t,i,j,isX)*(1.-2.*p)/p);
                        }
                        double z = (double)rand()/(double)RAND_MAX;
                        if (z <= p_jump) {
                            Npsi_p--;
                            Npsibar_p++;
                        }
                    }
                    // States in psibar
                    for (int ii = 0; ii < Npsibar_old; ++ii) {
                        double p_jump = gamma(t,i,j,isX)*dt;
                        if (gamma(t,i,j,isX) < 0.) {
                            if (p >= .5)
                                p_jump = dt*abs(gamma(t,i,j,isX)*(1.-2.*p)/(1.-p));
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
                    for (int ii = 0; ii < Npsi_old; ++ii) {
                        double p_jump = gamma(t,i,j,isX)*dt;
                        if (gamma(t,i,j,isX) < 0.) {
                            if (p >= .5)
                                p_jump = 0.;
                            else p_jump = dt*abs(gamma(t,i,j,isX)*(1.-2.*p)/p);
                        }
                        double z = (double)rand()/(double)RAND_MAX;
                        if (z <= p_jump) {
                            Npsi_m--;
                            Npsibar_m++;
                        }
                    }
                    // States in psibar
                    for (int ii = 0; ii < Npsibar_old; ++ii) {
                        double p_jump = gamma(t,i,j,isX)*dt;
                        if (gamma(t,i,j,isX) < 0.) {
                            if (p >= .5)
                                p_jump = dt*abs(gamma(t,i,j,isX)*(1.-2.*p)/(1.-p));
                            else p_jump = 0.;
                        }
                        double z = (double)rand()/(double)RAND_MAX;
                        if (z <= p_jump) {
                            Npsibar_m--;
                            Npsi_m++;
                        }
                    }
                }
            }
        }
    }

    
     return 0;
    
}