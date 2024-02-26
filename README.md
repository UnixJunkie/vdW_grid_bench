# vdW_grid_bench
Van der Waals 3D grid initialization benchmark

# problem statement
We want to know the interaction energy between one
atom [l_a] located at any grid point [i,j,k] and a fixed set of atoms [A].
We use the Van der Waals non-bonded interaction energy and
Universal Force Field (UFF) parameters.

# standard OCaml; before loop modification
time ./grid_init -dx 0.1 -l data/lig.pqrs -p data/prot.pqrs 
real    1m19.899s

# w/ dune's --profile=release
real    1m10.603s

# w/ V3 @@inline annotations plus dist2 optims
real    0m57.384s

# w/ optims in vdW_grid
real    0m49.096s

# w/ pow3 and V3.dist2 instead of pow6 and V3.dist in vdW_grid
real    0m41.282s
