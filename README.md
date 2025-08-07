# LSC_internship
Internship at the department of Scientific Computing Cambridge simulating matter under extreme conditions. 
## Compressible Euler Equations
Finite volume schemes modelling compressible fluids in 1 and 2 dimensions include:
- First Order schemes: FORCE, Lax-Freidrichs
- Second Order schemes (with limiting): SLIC, FLIC
- Godunov methods: Exact Riemann solver, HLL, HLLC
## Magnetohydrodynamics (MHD)
- Slope limiting scheme in 1D and 2D (SLIC)
- HLLC to approximate exact solution with 3 waves
- MUSCL-hancock in 1D and 2D with mixed GLM divergence cleaning
