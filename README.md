

# Scatter2

Perform 2D multistate scattering calculation in a Cartesian grid basis.


## Installation

Requirements:
- Armadillo > 9.500.2
- C++11 supported compiler (ideally C++14)
- OpenMP support
- Python3

To install, run (takes less than 1 min)

    make

## Usage

The argument list can be found by specifying `--help`.

File "linearpot.hpp" and "angularpot.hpp" contains various potential classes available for simulation. You may create your own classes of potential in these two files, which will automatically appear in the choices of "--potential" if properly set.

## Example

[example/bifurcation.sh](example/bifurcation.sh) generates the raw data for Fig.1 in [this paper](https://arxiv.org/abs/2008.02443). It may be modified to a PBS script.

## Citation

If you use this code, please consider citing

[1] [Yanze Wu, Gaohan Miao, and Joseph E. Subotnik. The Journal of Physical Chemistry A 2020 124 (37), 7355-7372](http://dx.doi.org/10.1021/acs.jpca.0c04562);

Feel free to contact us if you have any questions.


