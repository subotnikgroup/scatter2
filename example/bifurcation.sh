
# Parameters you may specify:
# E: Incoming energy
# A, eps, epsc, r, omega: Parameter of potential

# Requires a single core and 6.2 MB memory. Takes several minutes.

scatter -N 128 -L 8 --n_stencil=5 --n_ext=5 -m 1000 -E $E --in_el=0 --n_bs=20 --potential 2statelx --startpos down --output up-E$E --param $A,$mu,$lambda,$eps,$epsc,$omega,$r,1000,10 > out-E$E.txt 

# Will output three files:
# out-*: STDOUT;
# prob-*: Probability file. Describes the (normalized) scattering probability in each bound states and (unnormalized) scattering amplitude in each bound states.
# bound-*: Wavefunction of bound states.

