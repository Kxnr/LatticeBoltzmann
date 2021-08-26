#pragma once

#include <iostream>
#include <armadillo>
#include <string>

#define MAX_TIMES 200
#define TIME 5

// TODO: load lattice description from file
const struct LatticeConstants {
    float c = 1.0f / std::sqrt(3.0f); // dimensionless speed of sound
    float l0 = 1.0;         // meters
    float u0 = .1;          // m/s
    float t0 = l0 / u0;     // s
    float v = .0000000001;  // viscosity in m**2 / s
};

// Lattice constants
#define C 1.0f / std::sqrt(3.0f) //dimensionless speed of sound

// Tunable Params
const float l0 = 1.0;         // meters
const float u0 = .1;          // m/s
constexpr float t0 = l0 / u0; // s
const float v = .0000000001;  // viscosity in m**2 / s

// Derived params
constexpr float dx = 1.0f / WIDTH;
constexpr float dt = std::pow(dx, 2.0f); // gives ~ 2nd order accuracy
constexpr float steps = 1.0f / dt;         // simulation steps in one characteristic time

constexpr float Re = u0 * l0 / v;                     // Reynolds number
constexpr float vlb = dt / (std::pow(dx, 2) * Re); // lattice viscosity
constexpr float ulb = dt / dx;

constexpr float omega = 1.0f / (vlb / std::pow(C, 2.0f) + .5f); // 1 / relaxation parameter


//TODO: these are properties of a particular DnQm simulation
arma::imat vel = {{0, 1, 0, -1, 0, 1, -1, -1, 1},
                  {0, 0, 1, 0, -1, 1, 1, -1, -1}};
arma::uword opp[Q] = {0, 3, 4, 1, 2, 7, 8, 5, 6};

// vectors that can carry across boundaries
arma::uword bt[Q] = {4, 7, 8};
arma::uword bb[Q] = {2, 5, 6};
arma::uword bl[Q] = {1, 5, 8};
arma::uword br[Q] = {3, 6, 7};

arma::vec w = {4.0f / 9, 1.0f / 9, 1.0f / 9, 1.0f / 9, 1.0f / 9, 1.0f / 36, 1.0f / 36, 1.0f / 36, 1.0f / 36};

// TODO: boundary conditions
// TODO: computational backend

// necessary operations: dot product, subtraction, addition, indexing

template <size_t...Dims>
class Lattice {
    /*
     * Class to represent the Boltzmann Lattice. The first dimensions represent spatial dimensions,
     * and the last dimension represents the distribution over velocities. All operations are defined
     * either as uniform over space or uniform over velocities--this is not a general tensor class.
     */

    public:
        Lattice(Vector velocities);

        virtual int dims();
        virtual auto size();

        virtual auto sum(int axis);
        virtual auto stream(); // updates lattice according to boundary conditions and velocities

        // operators
        virtual Lattice operator*(Lattice& a, Lattice& b); // Second operand must match all dimensions, spatial dimensions, or velocity dimension
        virtual Lattice operator*(Lattice& a, double& b);
        virtual Lattice operator/(Lattice& a, Lattice& b);
        virtual Lattice operator/(Lattice& a, double& b);
        virtual Lattice operator+(Lattice& a, Lattice& b);
        virtual Lattice operator+(Lattice& a, double& b);
        virtual Lattice operator-(Lattice& a, Lattice& b);
        virtual Lattice operator-(Lattice& a, double& b);
        virtual Lattice operator%(Lattice& a, Lattice& b); // Hadamard/Schur product
};

template <size_t H, size_t W, size_t Q>
class 2DArmadilloLattice : Lattice<H, W, Q> {

};

class LatticeSimulator {
    public:
        LatticeSimulator(Lattice init, LatticeConstants params);

        Lattice lattice;

        LatticeSimulator step();

    private:
        const Lattice vel;
        Lattice equilibrium;

//        void equilibrium(arma::mat rho, arma::mat u, arma::cube &feq);
//
//        template <size_t S>
//        void bounceback(arma::cube fin, double (*indices)[S][2], arma::cube &fout);
//
//        void collide(arma::cube fin, arma::cube feq, arma::cube &fout);
//
//        void stream(arma::cube fin, arma::cube &fout);
//
//        void velocity(arma::cube fin, arma::mat rho, arma::cube &u);

        // combines equilibrium and collide steps
        virtual LatticeSimulator equilibrate();

        virtual LatticeSimulator bounceback();

        virtual LatticeSimulator stream();

        virtual void save();

        LatticeSimulator step();

};

class D2Q9LatticeSimulator : LatticeSimulator {

};