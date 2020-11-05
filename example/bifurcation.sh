
# Parameters you may specify:
# E: Incoming energy
# A, mu, lambda, eps, epsc, r, omega: Parameter of potential

A=0.02
mu=0.0002
lambda=0.001
eps=2.5
epsc=2.5
r=2
omega=0.01
E=0.007  # This is just one point. To repeat full Fig.1 you need to iterate E from 0.006 to 0.06.

# Requires a single core and 6.2 MB memory. Takes several minutes.

# Put scatter into your environment, or just move the script to project root.

scatter -N 128 -L 8 --n_stencil=5 --n_ext=5 -m 1000 -E $E --in_el=0 --n_bs=20 --potential 2statelx --startpos down --output up-E$E --param $A,$mu,$lambda,$eps,$epsc,$omega,$r,1000,10 > out-E$E.txt 

# Will output four files:
# out-*: STDOUT;
# prob-*: Probability file. Describes the (normalized) scattering probability in each bound states and (unnormalized) scattering amplitude in each bound states.
# bound-*: Wavefunction of bound states.
# state-*: Wavefunction of the solution.

