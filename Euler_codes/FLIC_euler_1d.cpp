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

std::array<double, 3>  getFluxFORCE(std::array<double, 3>  x , std::array<double, 3>  y, double dx , double dt,double gamma){
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


std::array<double, 3>  getFluxRI(std::array<double, 3>  x , std::array<double, 3>  y, double dx , double dt,double gamma){
    //impliment general flux function 
    std::array<double, 3>  f_1 = flux_def(x, gamma); //f(ui)
    std::array<double, 3>  f_2 = flux_def(y, gamma); //f(i+1)
    std::array<double, 3>  uPlusHalf; //u_{i+1/2}

    for(int i=0; i<=2; ++i){
        uPlusHalf[i] = 0.5*(x[i] + y[i]) - 0.5*(dt/dx)*(f_2[i] - f_1[i]);
    }

    std::array<double, 3> RI_flux = flux_def(uPlusHalf, gamma); //richtmyer flux
    
    return RI_flux;
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


//the inputs of this function are w =u_{i-1} x=u_{i} , y=u_{i+1} and z=u{i+2}


double phi(std::array<double, 3> w, std::array<double, 3> x, std::array<double, 3> y, std::array<double, 3> z, double C) {
    double phi_g = (1 - C) / (1 + C);

    // Calculate the slopes 
    double delta_q_i_minus_half = x[2] - w[2];
    double delta_q_i_plus_half = y[2] - x[2];
    double delta_q_i_plus_3half = z[2] - y[2];
    const double epsilon = 1e-20;

    double rL = delta_q_i_minus_half/ (delta_q_i_plus_half + epsilon);
    double rR = delta_q_i_plus_half/(delta_q_i_plus_3half+epsilon);

    // Lambda to compute limiter value based on r and phi_g
    auto computeLimiter = [phi_g](double r) -> double {
        if (r <= 0) return 0;
        if (r <= 1) return 2 * r / (1 + r);
        return phi_g + 2 * (1 - phi_g) * r / (1 + r);
    };

    double phi_rL = computeLimiter(rL);
    double phi_rR = computeLimiter(rR);

    return std::min(phi_rL, phi_rR);
}


//the inputs of this function are w =u_{i-1} x=u_{i} , y=u_{i+1} and z=u{i+2}

std::array<double , 3> FLIC_flux(std::array<double , 3> w , std::array<double , 3> x , std::array<double , 3> y  , std::array<double , 3> z , double C , double dx , double dt,double gamma){
    std::array<double , 3> FLIC_flux;
    for(int i=0 ; i<=2 ; ++i){
        FLIC_flux[i] = getFluxFORCE( x ,  y,  dx ,dt,gamma)[i] + phi(w,x,y,z,C) * (getFluxRI(x,y,dx,dt,gamma)[i] - getFluxFORCE(x,y,dx,dt,gamma)[i]);
    }
    return FLIC_flux;
}



int main() { 
    int nCells = 100; //the distance between points is 0.01
    double x0 = 0.0;
    double x1 = 1.0;
    double tStart = 0.0; //set the start and finish time steps the same
    double tStop = 0.15;
    double C = 0.6;
    double gamma = 1.4;

    // Allocate matrices with 4 extra points for transmissive BCs
    std::vector<std::array<double,3>> u(nCells+4);
    std::vector<std::array<double,3>> uPlus1(nCells+4);
    std::vector<std::array<double,3>> flux(nCells+4);

    double dx = (x1 - x0) / nCells; //the space steps 

    // Initial conditions


    for(int i = 0; i < u.size(); i++) {
        // x 0 is at point i=1/2
        double x = x0 + (i-0.5) * dx;
        std::array<double, 3> prim;
        if(x <= 0.5) {
            prim[0] = 1; // Density
            prim[1] = -2; // Velocity
            prim[2] = 0.4; // Pressure
            } else {
            prim[0] = 1; // Density
            prim[1] = 2; // Velocity
            prim[2] = 0.4; // Pressure
        }

        u[i] = primitiveToConservative(prim, gamma);
    }


    double dt = computeTimeStep(u , C , dx, gamma); //the time steps

    double t = tStart;
    do {
        // Compute the stable time step for this iteration

        dt = computeTimeStep(u , C , dx, gamma); 
        t = t + dt;
        //std::cout<<dt<<std::endl;
        //You may want to manually reduce dt if this would overshoot tStop
        //Apply boundary conditions

        // Trasmissive boundary conditions
        
        
        u[1] = u[2];
        u[0] = u[1];
        u[nCells+2] = u[nCells +1];
        u[nCells + 3] = u[nCells+2];


        for(int i = 1; i < nCells+2; i++) { //Define the fluxes
            // flux[i] corresponds to cell i+1/2 
            flux[i] = FLIC_flux( u[i-1],u[i], u[i+1] , u[i+2],C, dx , dt, gamma);
        }



        //the below has a i-1 flux which means we need to define a flux at 0 so make sure the above ^ starts at 0! this is because we have another edge with the number of cells (like the walls)

        for(int i = 2; i < nCells+2; i++) { //Update the data
            for(int j=0; j<=2; ++j){
                uPlus1[i][j] = u[i][j] - (dt/dx) * (flux[i][j] - flux[i-1][j]);
            }
        }
        
        uPlus1[1] = uPlus1[2];
        uPlus1[0] = uPlus1[1];
        uPlus1[nCells+2] = uPlus1[nCells +1];
        uPlus1[nCells + 3] = uPlus1[nCells+2];

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