# Fortran Full FEA Solver for Upsetting Process in Metal Forming
This repository provides a standalone Fortran code for simulating metal forming processes through Finite Element Analysis (FEA), specifically optimized for rigid viscoplastic materials in upsetting processes. Known as SPID (Simple Plastic Incremental Deformation), this solver provides a tailored solution for analyzing plastic deformation in metal forming.

## Features
- Material Model: Implements a rigid-viscoplastic material formulation to model plastic deformation accurately.
- Custom Solver: Independent, in-house solver for stiffness matrix assembly, boundary condition applications, and Gaussian elimination – no external FEA tools required.
- Element Handling: Supports 4-node quadrilateral elements in both 2D plane strain and axisymmetric modes.
- Boundary Conditions: Capable of handling displacement, force, and friction boundary conditions.
- Iterative Convergence Control: Uses both direct and Newton-Raphson iterative methods to handle nonlinear material behavior.
- Restart Capability: Generates a restart file (```SPID.RST```) to save simulation state, allowing for resumption of interrupted simulations.
- Post-Processing: Computes strain, stress, nodal displacements, and outputs for detailed analysis.

## Code Structure

- Main Program (```SPID.FOR```): Runs the simulation, handles input/output, and calls key subroutines for FEA operations.
- Core Subroutines:
  - STIFF: Global stiffness matrix and load vector assembly.
  - ADDBAN: Assembles element-level stiffness into the global matrix.
  - NONLIN: Manages nonlinear convergence and iterative updates.
  - DISBDY: Applies boundary conditions.
  - POTSOL: Handles post-solution updates for stress and strain evaluations.
  - RSTFIL: Creates a restart file (```SPID.RST```) to save the simulation state.

## Files and Usage 

- ```SPID.DAT```: Input file with initial parameters (node coordinates, element connectivity, boundary conditions, material properties).
- ```SPID.RST```: Optional restart file. Allows simulation resumption from a specific state, saving nodal positions, strain history, and current step data.
- ```SPID.OUT```: Main output file containing simulation results.
- ```SPID.MSG```: Log file with iteration, convergence, and error messages.

## Compile 
```sh
gfortran -o SPID SPID.FOR
./SPID
