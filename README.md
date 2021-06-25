# Lattice Boltzmann Simulation

## Description

Implementation of a Lattice Boltzmann simulation, originally used for a Physics 114 - Statistical Mechanics course at Swarthmore College to demonstrate both how the logic behind Ising models can be extended to more general systems and how endemic characteristics of simluations can create counterintuitive results--such as negative absolute temperature. For 2D Lattice Bolzmann simulations, one such result is that a random flow field will decay into large scale structures. This despite the random initial condition seeming to be much more statistically likely than any state with persistent large scale structure. Work done with the supervision of Prof. Michael 'Doc' Brown at Swarthmore College.

## Setup

### Python

#### Dependencies

* python >= 3
* numpy
* matplotlib

#### Executing program

```
python lbm2.py
```

### C++

#### Dependencies

* C++11 compatible compiler
* make
* armadillo

#### Installation

* Install armadillo
* Clone this repository
* Edit link instructions in Makefile as necessary for your installation of armadillo
* Run `make all`

#### Executing program

Output is produced in the form of text files, which represent the simulation grid at each time point. These values can be read by the `fviz.py` script to produce pressure, velocity, and vorticity diagrams.

```
./main.out
```

### Rust

Coming Soon

### Help

TODO: parameter description

## Authors

* [Connor Keane](kxnr.me)

## License

This project is licensed under the GPL3 License - see the LICENSE.md file for details

## Acknowledgments

* Original python implementation provided by [flowkit.com](FlowKit) under GPL3 and is replicated with modifications here
* Much more in depth work into this phenomemon published [http://plasma.physics.swarthmore.edu/brownpapers/BrownJPP96.pdf](here)
