

# Scatter2

Thie is an implementation of quantum scattering calculation in two dimension system in a Cartesian grid basis. The basic algorithm can be found in [J. Chem. Phys. **142**, 084110 (2015)](http://dx.doi.org/10.1063/1.4908032) and [J. Phys. Chem. A 2020 **124** (37), 7355-7372](http://dx.doi.org/10.1021/acs.jpca.0c04562). Please cite the two papers if you have used the code.


## Installation

Requirements:
- Armadillo > 9.500.2: Which can be found in [http://arma.sourceforge.net](http://arma.sourceforge.net).
- CXXOPT library: Which can be found in [https://github.com/jarro2783/cxxopts](https://github.com/jarro2783/cxxopts). Put cxxopts.hpp in the project root directory.
- A compiler supporting C++11 and OpenMP.

To build the code, run (should takes less than 1 min)

    make

## Usage

Arguments available to the code:

- N: Number of the grid;
- Ny: Number of the grid in the y direction (by default = N);
- dx: Spacing of the grid;
- n_stencil: Differential stencil used for the kinetic matrix. It is half of the actual stencil used, e.g. n_stencil=2 indicates a 3-stencil, n_stencil=3 indicates a 5-stencil, etc.
- n_ext: Extra grid for the boundary wave basis. By default is 3. Usually set to = n_stencil;
- n_bs: Limit of bound state used. Default is 1.
- in_bs: Incoming bound state. Default is 0.
- in_el: Incoming electronic state. 
- L: Length of the model in x direction. In atomic units.
- m: Nuclear mass. 
- E: Overall incoming energy (including potential energy).
- startpos: Incident position. Either "up" or "down". For angular potentials, specifying "up" indicates incoming from left.
- potential: Name of the potential used. Can be the return of `get_name()` of any potential included in the [angularpot.hpp](angularpot.hpp) and [linearpot.hpp](linearpot.hpp). To repeat our published work, only "2statelx" and "angular" are needed.
- param: Params for the model Hamiltonian. The format can be found in the example.
- output: The suffix of the output.

The argument list can also be found by specifying `--help`.

Running the program will give three outputs:

- prob-[suffix]: Describes the (normalized) scattering probabilities for each bound states and (unnormalized) scattering amplitudes for each bound states.
- bound-[suffix]: The wavefunction of bound states, in a grid basis.
- state-[suffix]: The wavefunction solution. Will output as a (N x Nel) row and N column matrix. The number in (a*N + i)'th row and j'th column indicates the wavefunction value at electronic state |a> and nuclear position (i, j).

## Example

[example/bifurcation.sh](example/bifurcation.sh) generates the two-dimensional scattering result in a bifurcation Hamiltonian described in [this paper](https://arxiv.org/abs/2008.02443). It should be able to run directly if the program is correctly built.

