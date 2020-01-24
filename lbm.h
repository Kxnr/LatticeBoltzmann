//lbmh
#ifndef LBM_PROTECT
#define LBM_PROTECT

#include <iostream>
#include <armadillo>
#include <string>

#define MAX_TIMES 200
#define TIME 5

// Lattice constants
#define WIDTH 1024
#define HEIGHT 1024
#define Q 9
#define C 1.0f / std::sqrt(3.0f) //dimensionless speed of sound

// Tunable Params
float l0 = 1.0;        // meters
float u0 = .1;         // m/s
float t0 = l0 / u0;    // s
float v = .0000000001; // viscosity in m**2 / s

// Derived params
float dx = 1.0f / WIDTH;
float dt = std::pow(dx, 2.0f); // gives ~ 2nd order accuracy
float steps = 1.0f / dt;       // simulation steps in one characteristic time

float Re = u0 * l0 / v;                  //Reynolds number
float vlb = dt / (std::pow(dx, 2) * Re); // lattice viscosity
float ulb = dt / dx;

float omega = 1.0f / (vlb / std::pow(C, 2.0f) + .5f); // 1 / relaxation parameter

arma::imat vel = {{0, 1, 0, -1, 0, 1, -1, -1, 1},
                  {0, 0, 1, 0, -1, 1, 1, -1, -1}};
arma::uword opp[Q] = {0, 3, 4, 1, 2, 7, 8, 5, 6};
arma::uword bt[Q] = {4, 7, 8};
arma::uword bb[Q] = {2, 5, 6};
arma::uword bl[Q] = {1, 5, 8};
arma::uword br[Q] = {3, 6, 7};

arma::vec w = {4.0f / 9, 1.0f / 9, 1.0f / 9, 1.0f / 9, 1.0f / 9, 1.0f / 36, 1.0f / 36, 1.0f / 36, 1.0f / 36};

void equilibrium(arma::mat rho, arma::mat u, arma::cube &feq);

void sumpop(arma::cube pop, arma::mat &rho);

template <size_t S>
void bounceback(arma::cube fin, double (*indices)[S][2], arma::cube &fout);

void collide(arma::cube fin, arma::cube feq, arma::cube &fout);

void stream(arma::cube fin, arma::cube &fout);

void velocity(arma::cube fin, arma::mat rho, arma::cube &u);

#endif