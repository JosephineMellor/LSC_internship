#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <array>
#include <filesystem>

//need to be able to convert between the primitive and conservative variables both ways

std::array<double, 3> primitiveToConservative(std::array<double , 3> prim , double gamma){
    std::array<double,3> consv ;
    consv[0] = prim [0];
    consv[1] = prim [0] * prim [1];
    consv[2] = prim [2] /(gamma -1) + 0.5*prim[0]*prim[1]*prim[1];
    return consv;
}

std::array<double, 3> conservativeToPrimative(std::array<double, 3> x , double gamma){
    std::array<double,3> prim;
    prim[0] = x[0];
    prim[1] = x[1] / x[0];
    prim[2] = (gamma -1)*( x[2] - 0.5*prim[1]*prim[1]*prim[0]);
    return prim;
}


//define a general flux function 

std::array<double, 3> flux_def(std::array<double, 3> x, double gamma){
    std::array<double,3> flux;
    double rho = x[0];
    double u = x[1] / rho;
    double KE = 0.5 * rho * u * u;
    double p = (gamma - 1) * (x[2] - KE);

    flux[0] = x[1];               // mass: rho u
    flux[1] = rho * u * u + p;    // momentum: rho u² + p
    flux[2] = (x[2] + p) * u;     // energy: (E + p)u
    return flux;
}

//define a function to get the flux x is ui and y is ui+1 they are both arrays of length 3
//it will spit out an array of size 3 as well give each variable their flux

std::array<double, 3>  getFlux(std::array<double, 3>  x , std::array<double, 3>  y, double dx , double dt,double gamma){
    //impliment general flux function 
    std::array<double, 3>  f_1 = flux_def(x, gamma); //f(ui)
    std::array<double, 3>  f_2 = flux_def(y, gamma); //f(i+1)
    std::array<double, 3>  uPlusHalf; //u_{i+1/2}

    for(int i=0; i<=2; ++i){
        uPlusHalf[i] = 0.5*(x[i] + y[i]) - 0.5*(dt/dx)*(f_2[i] - f_1[i]);
    }

    std::array<double, 3> RI_flux = flux_def(uPlusHalf, gamma); //richtmyer flux
    std::array<double, 3> LF_flux; //set up 3 length array for LF flux
    std::array<double, 3> FORCE_flux; //set up 3 length array for FORCE

    for (int i=0; i<=2; ++i) {
        LF_flux[i] = (0.5 * (dx/dt) * (x[i] - y[i])) + 0.5 * (f_1[i] + f_2[i]);
        FORCE_flux[i] = 0.5 * (LF_flux[i] + RI_flux[i]);
    }

    return FORCE_flux;
}


double computeTimeStep(const std::vector<std::array<double,3>>& u , double C, double dx, double gamma) {
    double maxSpeed = 0.0;

    for (const auto& state : u) {
        double rho = state[0];
        double mom = state[1];
        double E = state[2];

        double u_val = mom / rho;
        double KE = 0.5 * rho * u_val * u_val;
        double pressure = (gamma - 1.0) * (E - KE);

        double sound_speed = std::sqrt(gamma * pressure / rho);
        double speed = std::abs(u_val) + sound_speed;

        if (speed > maxSpeed) {
            maxSpeed = speed;
        }
    }

    double dt = C * dx / maxSpeed;
     return dt;
}

int counter =0;

int main() { 
    int nCells = 100; //the distance between points is 0.01
    double x0 = 0.0;
    double x1 = 1.0;
    double tStart = 0.0; //set the start and finish time steps the same
    double tStop = 0.15;
    double C = 0.8;
    double gamma = 1.4;

    // Allocate matrices with 2 extra points for transmissive BCs
    std::vector<std::array<double,3>> u(nCells+2);
    std::vector<std::array<double,3>> uPlus1(nCells+2);
    std::vector<std::array<double,3>> flux(nCells+2);
    std::vector<double> time(nCells);
    double dx = (x1 - x0) / nCells; //the space steps 

    // Initial conditions!

    for(int i = 0; i < u.size(); i++) {
        // x 0 is at point i=1/2
        double x = x0 + (i-0.5) * dx;
        std::array<double, 3> prim;
        if(x <= 0.25) {
            prim[0] = 1; // Density
            prim[1] = 1; // Velocity
            prim[2] = 1; // Pressure
            } else {
            prim[0] = 0.1; // Density
            prim[1] = 1; // Velocity
            prim[2] = 1; // Pressure
        }

        u[i] = primitiveToConservative(prim, gamma);
    }

    // for(int i =0 ; i< u.size(); i++){
    //     std::array<double, 3> prim;
    //     prim[0] = 2;
    //     prim[1] = 1;
    //     prim[2] = 3;

    //     u[i] = primitiveToConservative(prim, gamma);
    // }


    double dt = computeTimeStep(u , C , dx, gamma); //the time steps

    double t = tStart;
    do {
        // Compute the stable time step for this iteration

        dt = computeTimeStep(u , C , dx, gamma); 
        t = t + dt;

        //You may want to manually reduce dt if this would overshoot tStop
        //Apply boundary conditions

        // Trasmissive boundary conditions
        u[0] = u[1];
        u[nCells + 1] = u[nCells];


        for(int i = 0; i < nCells+1; i++) { //Define the fluxes
            // flux[i] corresponds to cell i+1/2 
            flux[i] = getFlux( u[i], u[i+1] , dx , dt, gamma);
        }

        //the below has a i-1 flux which means we need to define a flux at 0 so make sure the above ^ starts at 0! this is because we have another edge with the number of cells (like the walls)

        for(int i = 1; i < nCells+1; i++) { //Update the data
            for(int j=0; j<=2; ++j){
                uPlus1[i][j] = u[i][j] - (dt/dx) * (flux[i][j] - flux[i-1][j]);
            }
        }
    
        // Now replace u with the updated data for the next time step

        u = uPlus1;
    } while (t < tStop);

    //still need to convert it back to primitive

    //define final results

    std::vector<std::array<double,3>> results(nCells+2);

    for(int i=0; i<= results.size() -1; ++i){
        results[i] = conservativeToPrimative(u[i], gamma);
    }


    // Output the results
    std::ofstream output("euler.dat");
    for (int i = 1; i <= nCells; ++i) {
        double x = x0 + (i - 1) * dx;
        output << x << " " << results[i][0] <<  " " << results[i][1] <<  " " << results[i][2] << std::endl;
        // std::cout << x << " " << u[i][0] <<  " " << u[i][1] <<  " " << u[i][2] << std::endl;
    }
}
