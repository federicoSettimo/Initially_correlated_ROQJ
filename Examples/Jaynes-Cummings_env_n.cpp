// Jaynes-Cummings model, particular case that env states (in B+) are eigenstates of n 
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

double one_X_one (const Matrix2cd X) {return real(X(0,0));}
double one_X_zero_re (const Matrix2cd X) {return real(X(0,1));}
double one_X_zero_im (const Matrix2cd X) {return imag(X(0,1));}

// Parameters of the model
double Delta = 1., g = .2;
double tmax = 5., dt = .001;
int Nstates = 10000;
int n_0 = 1, n_x = 2, n_y = 10, n_z = 3;

// Functions and dynamics
double Omega (int n) {return sqrt(Delta*Delta + 4.*g*g*n);}

// Eigenvalues for the unitary dynamics
complex<double> c (int n, double t) {return exp(.5*I*Delta*t)*(cos(.5*Omega(n)*t) - I*Delta*sin(.5*Omega(n)*t)/Omega(n));}
complex<double> cdot (int n, double t) {return exp(.5*I*Delta*t)*sin(.5*t*Omega(n))*(Delta*Delta - Omega(n)*Omega(n))/(2.*Omega(n));}
double absc2 (int n, double t) {return norm(c(n,t));}
complex<double> conjc (int n, double t) {return conj(c(n,t));}
complex<double> conjc_dot (int n, double t) {return exp(-.5*I*Delta*t)*sin(.5*t*Omega(n))*(Delta*Delta - Omega(n)*Omega(n))/(2.*Omega(n));}
double absc2dot (int n, double t) {return real(cdot(n,t)*conjc(n,t) + c(n,t)*conjc_dot(n,t));}

complex<double> d (int n, double t) {return -I*exp(.5*I*Delta*t)*2.*g*sin(.5*Omega(n)*t)/Omega(n);}
complex<double> ddot (int n, double t) {return exp(.5*I*Delta*t)*g*(Delta*sin(.5*Omega(n)*t) + I*Omega(n)*cos(.5*Omega(n)*t))/Omega(n);}
double absd2 (int n, double t) {return norm(d(n,t));}
complex<double> conjd (int n, double t) {return conj(d(n,t));}
complex<double> conjd_dot (int n, double t) {return exp(-.5*I*Delta*t)*g*(Delta*sin(.5*Omega(n)*t) - I*Omega(n)*cos(.5*Omega(n)*t))/Omega(n);}
double absd2dot (int n, double t) {return real(ddot(n,t)*conjd(n,t) + d(n,t)*conjd_dot(n,t));}

// For sigma_x
complex<double> kappa (int n, double t) {return c(n+1,t)*c(n,t);}
complex<double> kappadot (int n, double t) {return cdot(n+1,t)*c(n,t) + c(n+1,t)*cdot(n,t);}
double b_x (int n, double t) {return -.5*imag(kappadot(n,t)/kappa(n,t));}
double gamma_x (int n, double t) {return -.5*real(kappadot(n,t)/kappa(n,t));}

// For sigma_y
complex<double> zeta (int n, double t) {return -I*c(n+1,t)*c(n,t);}
complex<double> zetadot (int n, double t) {return -I*(cdot(n+1,t)*c(n,t) + c(n+1,t)*cdot(n,t));}
double b_y (int n, double t) {return -.5*imag(zetadot(n,t)/zeta(n,t));}
double gamma_y (int n, double t) {return -.5*real(zetadot(n,t)/zeta(n,t));}

// For sigma_z
double xi (int n, double t) {return absc2(n+1,t) - n*absd2(n,t);}
double gamma_m_z (int n, double t) {return -(xi(n,t+dt) - xi(n,t))/(xi(n,t)*dt);}

// For Q0
double beta (int n, double t) {return absc2(n+1,t) + n*absd2(n,t) - 1.;}
double a (int n, double t) {return beta(n,t) - xi(n,t);}
double adot (int n, double t) {return (a(n,t+dt)-a(n,t))/dt;}
complex<double> b (int n, double t) {return kappa(n,t)+zeta(n,t);}
double Re_bdotb (int n, double t) {return real((b(n,t+dt)-b(n,t))/(dt*b(n,t)));}
double Im_bdotb (int n, double t) {return imag((b(n,t+dt)-b(n,t))/(dt*b(n,t)));}



int main () {
    ofstream out_ex, out_avg, out_params, out_gamma;
    out_params.open("params.txt");
    out_params << tmax << endl << dt << endl << Delta << endl << g << endl << n_0 << endl << n_x << endl << n_y << endl << n_z << endl;


    // Q_x
    out_ex.open("exact_x.txt");
    out_avg.open("avg_x.txt");
    out_gamma.open("gamma_x.txt");
    Vector2cd psi_p = plus_state, psibar_p = minus_state;
    Vector2cd psi_m = minus_state, psibar_m = plus_state;
    Matrix2cd X_ex;

    int Npsi_p = Nstates, Npsibar_p = 0, Npsi_m = Nstates, Npsibar_m = 0;

    for (double t = 0.; t < tmax; t += dt) {
        // Exact
        Matrix2cd X_ex;
        complex<double> k = kappa(n_x,t);
        X_ex(0,1) = k;
        X_ex(1,0) = conj(k);
        out_ex << one_X_one(X_ex) << " " << one_X_zero_re(X_ex) << " " << one_X_zero_im(X_ex) << endl;

        // Unraveling
        Matrix2cd Xavg = ((double)Npsi_p*projector(psi_p) +(double) Npsibar_p*projector(psibar_p) - (double)Npsi_m*projector(psi_m) - (double)Npsibar_m*projector(psibar_m))/(double)(Nstates);
        out_avg << one_X_one(Xavg) << " " << one_X_zero_re(Xavg) << " " << one_X_zero_im(Xavg) << endl;
        double g_x_t = gamma_x(n_x,t);
        out_gamma << g_x_t << " " << b_x(n_x,t) << endl;
        Matrix2cd K = b_x(n_x,t)*sigma_z - .5*I*g_x_t*sigma_z*sigma_z;
        psi_p -= I*K*dt*psi_p; psi_p.normalize();
        psibar_p -= I*K*dt*psibar_p; psibar_p.normalize();
        psi_m -= I*K*dt*psi_m; psi_m.normalize();
        psibar_m -= I*K*dt*psibar_m; psibar_m.normalize();
        // First with +
        int Npsi_old = Npsi_p, Npsibar_old = Npsibar_p;
        double p = (double)Npsi_p/(double)Nstates;
        // States in psi
        for (int i = 0; i < Npsi_old; ++i) {
            double p_jump = g_x_t*dt;
            if (g_x_t < 0.) {
                if (p >= .5)
                    p_jump = 0.;
                else p_jump = -dt*abs(g_x_t*(1.-2.*p)/p);
            }
            double z = (double)rand()/(double)RAND_MAX;
            if (z <= p_jump) {
                Npsi_p--;
                Npsibar_p++;
            }
        }
        // States in psibar
        for (int i = 0; i < Npsibar_old; ++i) {
            double p_jump = g_x_t*dt;
            if (g_x_t < 0.) {
                if (p >= .5)
                    p_jump = dt*abs(g_x_t*(1.-2.*p)/(1.-p));
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
            double p_jump = g_x_t*dt;
            if (g_x_t < 0.) {
                if (p >= .5)
                    p_jump = 0.;
                else p_jump = -dt*g_x_t*(1.-2.*p)/p;
            }
            double z = (double)rand()/(double)RAND_MAX;
            if (z <= p_jump) {
                Npsi_m--;
                Npsibar_m++;
            }
        }
        // States in psibar
        for (int i = 0; i < Npsibar_old; ++i) {
            double p_jump = g_x_t*dt;
            if (g_x_t < 0.) {
                if (p >= .5)
                    p_jump = dt*g_x_t*(1.-2.*p)/(1.-p);
                else p_jump = 0.;;
            }
            double z = (double)rand()/(double)RAND_MAX;
            if (z <= p_jump) {
                Npsibar_m--;
                Npsi_m++;
            }
        }

    }


    // Q_y
    out_ex.close();
    out_avg.close();
    out_gamma.close();
    out_ex.open("exact_y.txt");
    out_avg.open("avg_y.txt");
    out_gamma.open("gamma_y.txt");
    psi_p = plus_y; psibar_p = minus_y;
    psi_m = minus_y; psibar_m = plus_y;

    Npsi_p = Nstates; Npsibar_p = 0; Npsi_m = Nstates; Npsibar_m = 0;

    for (double t = 0.; t < tmax; t += dt) {
        // Exact
        Matrix2cd X_ex;
        complex<double> z = zeta(n_y,t);
        X_ex(0,1) = z;
        X_ex(1,0) = conj(z);
        out_ex << one_X_one(X_ex) << " " << one_X_zero_re(X_ex) << " " << one_X_zero_im(X_ex) << endl;

        // Unraveling
        Matrix2cd Xavg = ((double)Npsi_p*projector(psi_p) +(double) Npsibar_p*projector(psibar_p) - (double)Npsi_m*projector(psi_m) - (double)Npsibar_m*projector(psibar_m))/(double)(Nstates);
        out_avg << one_X_one(Xavg) << " " << one_X_zero_re(Xavg) << " " << one_X_zero_im(Xavg) << endl;
        double g_y_t = gamma_y(n_y,t);
        out_gamma << g_y_t << " " << b_y(n_y,t) << endl;
        Matrix2cd K = b_y(n_y,t)*sigma_z - .5*I*g_y_t*sigma_z*sigma_z;
        psi_p -= I*K*dt*psi_p; psi_p.normalize();
        psibar_p -= I*K*dt*psibar_p; psibar_p.normalize();
        psi_m -= I*K*dt*psi_m; psi_m.normalize();
        psibar_m -= I*K*dt*psibar_m; psibar_m.normalize();
        // First with +
        int Npsi_old = Npsi_p, Npsibar_old = Npsibar_p;
        double p = (double)Npsi_p/(double)Nstates;
        // States in psi
        for (int i = 0; i < Npsi_old; ++i) {
            double p_jump = g_y_t*dt;
            if (g_y_t < 0.) {
                if (p >= .5)
                    p_jump = 0.;
                else p_jump = -dt*abs(g_y_t*(1.-2.*p)/p);
            }
            double z = (double)rand()/(double)RAND_MAX;
            if (z <= p_jump) {
                Npsi_p--;
                Npsibar_p++;
            }
        }
        // States in psibar
        for (int i = 0; i < Npsibar_old; ++i) {
            double p_jump = g_y_t*dt;
            if (g_y_t < 0.) {
                if (p >= .5)
                    p_jump = dt*abs(g_y_t*(1.-2.*p)/(1.-p));
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
            double p_jump = g_y_t*dt;
            if (g_y_t < 0.) {
                if (p >= .5)
                    p_jump = 0.;
                else p_jump = -dt*g_y_t*(1.-2.*p)/p;
            }
            double z = (double)rand()/(double)RAND_MAX;
            if (z <= p_jump) {
                Npsi_m--;
                Npsibar_m++;
            }
        }
        // States in psibar
        for (int i = 0; i < Npsibar_old; ++i) {
            double p_jump = g_y_t*dt;
            if (g_y_t < 0.) {
                if (p >= .5)
                    p_jump = dt*g_y_t*(1.-2.*p)/(1.-p);
                else p_jump = 0.;
            }
            double z = (double)rand()/(double)RAND_MAX;
            if (z <= p_jump) {
                Npsibar_m--;
                Npsi_m++;
            }
        }

    }



    // Q_z
    out_ex.close();
    out_avg.close();
    out_gamma.close();
    out_ex.open("exact_z.txt");
    out_avg.open("avg_z.txt");
    out_gamma.open("gamma_z.txt");
    psi_p = excited_state; psibar_p = ground_state;
    psi_m = ground_state; psibar_m = excited_state;

    Npsi_p = Nstates; Npsibar_p = 0; Npsi_m = Nstates; Npsibar_m = 0;

    for (double t = 0.; t < tmax; t += dt) {
        // Exact
        Matrix2cd X_ex;
        complex<double> x = xi(n_z,t);
        X_ex = x*sigma_z;
        out_ex << one_X_one(X_ex) << " " << one_X_zero_re(X_ex) << " " << one_X_zero_im(X_ex) << endl;

        // Unraveling
        Matrix2cd Xavg = ((double)Npsi_p*projector(psi_p) +(double) Npsibar_p*projector(psibar_p) - (double)Npsi_m*projector(psi_m) - (double)Npsibar_m*projector(psibar_m))/(double)(Nstates);
        out_avg << one_X_one(Xavg) << " " << one_X_zero_re(Xavg) << " " << one_X_zero_im(Xavg) << endl;
        double g_z_t = gamma_m_z(n_z,t);
        out_gamma << g_z_t << endl;
        // First with +
        int Npsi_old = Npsi_p, Npsibar_old = Npsibar_p;
        double p = (double)Npsi_p/(double)Nstates;
        double gm = g_z_t, gp = 0.;
        if (gm < 0.) {
            if (1.-p == 0.)
                gp= 1./dt;
            else
                gp = -p*gm/(1.-p);
            gm = 0.;
        }
        // States in |1>
        for (int i = 0; i < Npsi_old; ++i) {
            double z = (double)rand()/(double)RAND_MAX;
            if (z <= gm*dt) {
                Npsi_p--;
                Npsibar_p++;
            }
        }
        // States in |0>
        for (int i = 0; i < Npsibar_old; ++i) {
            double z = (double)rand()/(double)RAND_MAX;
            if (z <= gp*dt) {
                Npsibar_p--;
                Npsi_p++;
            }
        }

        // Now with -
        Npsi_old = Npsi_m; Npsibar_old = Npsibar_m;
        p = (double)Npsibar_m/(double)Nstates;
        gm = g_z_t; gp = 0.;
        if (gm < 0.) {
            if (1.-p == 0.)
                gp= 1./dt;
            else
                gp = -p*gm/(1.-p);
            gm = 0.;
        }
        // States in |0>
        for (int i = 0; i < Npsi_old; ++i) {
            double z = (double)rand()/(double)RAND_MAX;
            if (z <= gp*dt) {
                Npsi_m--;
                Npsibar_m++;
            }
        }
        // States in |1>
        for (int i = 0; i < Npsibar_old; ++i) {
            double z = (double)rand()/(double)RAND_MAX;
            if (z <= gm*dt) {
                Npsibar_m--;
                Npsi_m++;
            }
        }

    }



    

    // Q0
    out_ex.close();
    out_avg.close();
    out_gamma.close();
    out_ex.open("exact_0.txt");
    out_avg.open("avg_0.txt");
    out_gamma.open("gamma_0.txt");

    int N0 = 0, N1 = Nstates, Nm = Nstates, Np = 0, Nmy = Nstates, Npy = 0;
    Vector2cd plus_t = plus_state, minus_t = minus_state, plusy_t = plus_y, minusy_t = minus_y;
    for (double t = 0.; t < tmax; t += dt) {
        // Exact
        Matrix2cd X_ex = id + a(n_0,t)*sigma_z;
        X_ex(0,1) = -b(n_0,t);
        X_ex(1,0) = -conj(b(n_0,t));
        out_ex << one_X_one(X_ex) << " " << one_X_zero_re(X_ex) << " " << one_X_zero_im(X_ex) << endl;

        // Unraveling
        Matrix2cd Xavg = 2.*(-N0*projector(ground_state) - N1*projector(excited_state) + Nm*projector(minus_t) + Np*projector(plus_t) + Nmy*projector(minusy_t) + Npy*projector(plusy_t))/(double)(Nstates);
        out_avg << one_X_one(Xavg) << " " << one_X_zero_re(Xavg) << " " << one_X_zero_im(Xavg) << endl;

        // Rates
        double rk = -Re_bdotb(n_0,t), A = a(n_0,t), Adot = adot(n_0,t), beta = -.5*Im_bdotb(n_0,t);
        double gp = -.5*Adot + rk*(1.-A), gm = .5*Adot + rk*(1.+A); // For |0,1>
        double gpx = -Adot/(2.*A), gz = (Adot + 2.*A*rk)/(4.*A); // For |+->
        out_gamma << gp << " " << gm << " " << beta << " " << gz << " " << gpx << endl;

        // For the sigma_z part
        int N0old = N0, N1old = N1;
        double p = (double)N1/(double)Nstates, pj0 = gp*dt, pj1 = gm*dt;
        if (pj0 < 0. || pj1 < 0.) {
            if (p*gm - (1.-p)*gp >= 0.) {
                if (p == 0.)
                    pj1 = 1./dt;
                else
                    pj1 = (gm - (1.-p)/p*gp)*dt;
                pj0 = 0.;
            }
            else {
                if (1. - p == 0.)
                    pj0 = 1./dt;
                else
                    pj0 = (gp - p*gm/(1.-p))*dt;
                pj1 = 0.;
            }
        }
        for (int i = 0; i < N0old; ++i) {
            double z = (double)rand()/(double)RAND_MAX;     
            if (z <= pj0) {
                N0--;
                N1++;
            }
        }
        for (int i = 0; i < N1old; ++i) {
            double z = (double)rand()/(double)RAND_MAX;
            if (z <= pj1) {
                N1--;
                N0++;
            }
        }

        // For the sigma_x part
        int Npold = Np, Nmold = Nm;
        p = (double)Nm/(double)Nstates;
        double gg = gz + .5*gpx;
        double pjp = gg*dt, pjm = gg*dt;
        plus_t -= I*dt*beta*sigma_z*plus_t; plus_t.normalize();
        minus_t -= I*dt*beta*sigma_z*minus_t; minus_t.normalize();
        if (pjp < 0.) {
            if (1. - 2.*p >= 0.) {
                if (p == 0.)
                    pjm = 1./dt;
                else
                    pjm = -dt*gg*(1.-2.*p)/p;
                pjp = 0.;
            }
            else {
                if (1. - p == 0.)
                    pjp = 1./dt;
                else
                    pjp = dt*gg*(1.-2.*p)/(1.-p);
                pjm = 0.;
            }
        }
        for (int i = 0; i < Nmold; ++i) {
            double z = (double)rand()/(double)RAND_MAX;
            if (z <= pjm) {
                Nm--;
                Np++;
            }
        }
        for (int i = 0; i < Npold; ++i) {
            double z = (double)rand()/(double)RAND_MAX;
            if (z <= pjp) {
                Np--;
                Nm++;
            }
        }

        // For the sigma_y part - rates and everything is the same as the sigma_x part
        int Npyold = Npy, Nmyold = Nmy;
        plusy_t -= I*dt*beta*sigma_z*plusy_t; plusy_t.normalize();
        minusy_t -= I*dt*beta*sigma_z*minusy_t; minusy_t.normalize();
        p = (double)Nmy/(double)Nstates;
        pjp = gg*dt; pjm = gg*dt;
        if (pjp < 0.) {
            if (1. - 2.*p >= 0.) {
                if (p == 0.)
                    pjm = 1./dt;
                else
                    pjm = -dt*gg*(1.-2.*p)/p;
                pjp = 0.;
            }
            else {
                if (1. - p == 0.)
                    pjp = 1./dt;
                else
                    pjp = dt*gg*(1.-2.*p)/(1.-p);
                pjm = 0.;
            }
        }
        for (int i = 0; i < Nmyold; ++i) {
            double z = (double)rand()/(double)RAND_MAX;
            if (z <= pjm) {
                Nmy--;
                Npy++;
            }
        }
        for (int i = 0; i < Npyold; ++i) {
            double z = (double)rand()/(double)RAND_MAX;
            if (z <= pjp) {
                Npy--;
                Nmy++;
            }
        }
    }


    return 0;
}