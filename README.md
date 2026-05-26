# Shape Optimization Pipeline Using Parametric Reduced Order Models: An Application to Soft Swimmers

## Introduction

This repo contains the code used in the paper ``Shape Optimization Pipeline Using Parametric Reduced Order Models: An Application to Soft Swimmers'' by Mathieu Dubied, Paolo Tiso, and Robert K. Katzschmann [TODO: add link to paper once available]. 

The code uses a slightly modified version of the FEM code YetAnotherFEcode (https://github.com/jain-shobhit/YetAnotherFEcode, https://zenodo.org/records/15257935), and implements the algorithm presented in the paper.

## Requirements and installation

The main part of the code runs on Matlab, and tensors are constructed using Julia. To ensure compatibility between the two programming languages, the Mex package is used (https://juliapackages.com/p/mex). This code was tested with the following setup

1. Install Matlab2020b
2. Install version 1.9.2 of Julia from https://julialang.org/downloads/oldreleases/
3. Install Matlab.jl package in Julia from https://juliapackages.com/p/matlab.
  First, make sure to load the Julia package manager:
   ```julia
   using Pkg
   Pkg.add("MATLAB")
   ] build MATLAB
   ```
   To test the installation, run
   ```julia
   using MATLAB
   mat"2*2"
   ```
4. Install Mex.jl package (https://juliapackages.com/p/mex) in Julia.  
   First, build the package:
   ```julia
   ] build Mex
   ```
   Restart Julia and add the Julia path to Matlab. To test the installation run the following in Matlab terminal:
   ```matlab
   jl.eval('2+2')
   ```
Notes:
- If a Julia script is changed, close and re-open MATLAB to see the changes.
- For the Julia call in yetAnotherFEcode to work, you may need to install TensorOperations.jl.
If required, an error message will indicate this.
- Mathematica is used to obtain analytical tensors describing the hydrodynamic drag force, as functions. The Mathematica code does not need to be re-run as the tensors are stored as Matlab functions in `src/Hydrodynamics/3D/TET4D4-tensors`.

## Code description
To run the code, start by running `startup.m`. The examples presented in the paper are gathered in the folder `examples`, and in particular:
- `A_FOM_PROM_comparison.m`: compares the FOM to the (P)ROMs (Section 6.1 of the paper).
- `B_shape_optimisation.m`: optimizes the shape of the fish for multiple sets of search parameters (Section 6.2 of the paper).
- `C_co_optimisation.m`: co-optimizes the shape and the actuation signal's amplitude (Section 6.3 of the paper)

The algorithm and its components are implemented in the folder `src`, with the main file being `src/Optimization
/optimise_shape_3D.m`.




