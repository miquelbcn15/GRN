#!bin/bash

# Generates the code and results for the GRN problem

cd src
make realclean
make superall
gnuplot drawgrn.gnu
mv *.eps ../img
