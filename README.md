# LSC_internship
Internship at the department of Scientific Computing Cambridge simulating matter under extreme conditions. 
## Compressible Euler Equations
Finite volume schemes modelling compressible fluids in 1 and 2 dimensions include:
- First Order schemes: FORCE, Lax-Freidrichs
- Second Order schemes (with limiting): SLIC, FLIC
- Godunov methods: Exact Riemann solver, HLL, HLLC
## Magnetohydrodynamics (MHD)
Finite difference and finite volume schemes modelling plasma in 1 and 2 dimensions:
- Slope limiting scheme in 1D and 2D (SLIC)
- HLLC to approximate exact solution with 3 waves
- MUSCL-hancock in 1D and 2D with mixed GLM divergence cleaning
  <table>
    <tr>
      <td>
        <img src="MUSCL.gif" alt="Orzang-Tang Animation" width="300"/><br/>
        <p align="center"><em>Figure 1: Using MUSCL-Hancock to model Denisty for Orazang-Tang vortex</em></p>
      </td>
      <td>
        <img src="MUSCL2.gif" alt="Kelvin-Helmhotz Animation" width="300"/><br/>
        <p align="center"><em>Figure 2: Using MUSCL-Hancock to model Denisty for Kelvin-Helmhotz instability</em></p>
      </td>
    </tr>
  </table>
## Linear Solov'Ev equation solver
Solving the sparse linear equation to solve the Solov'Ev equation modelling a nuclear fusion reactor. The linear system was solved with the C++ Eigen library, as well as analytically.
<table>
    <tr>
      <td>
        <img src="solovev.gif" alt="Solov'Ev Solution" width="300"/><br/>
        <p align="center"><em>Figure 3: Solov'Ev solution with a linear solver</em></p>
      </td>
      <td>
        <img src="solovev_exact.gif" alt="Exact Solov'Ev Solution" width="300"/><br/>
        <p align="center"><em>Figure 4: Solov'Ev exat solution</em></p>
      </td>
    </tr>
  </table>
