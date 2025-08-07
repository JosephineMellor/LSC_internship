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


//get the HLLC intermediate 
//here, x is u_{i} and y is u_{i+1}

std::array<double , 6> uHLLC ( std::array<double , 3>x , std::array<double , 3> y , double gamma){
    double rhoL = x[0];
    double rhoR = y[0];
    double vL = x[1] / x[0];
    double vR = y[1] / y[0];
    double pL = (gamma -1)*(x[2] - 0.5*vL*vL*rhoL);
    double pR = (gamma -1)*(y[2] - 0.5*vR*vR*rhoR);
    
    double Splus = std::abs(vL) + std::sqrt(gamma *  pL / rhoL);
    if(std::abs(vR) + std::sqrt(gamma *  pR / rhoR) > Splus){
        Splus = std::abs(vR) + std::sqrt(gamma *  pR / rhoR);
    }

    double sL = -Splus;
    double sR = Splus;

    double Sstar = (pR - pL + rhoL*vL*(sL - vL) - rhoR*vR*(sR - vR)) / (rhoL*(sL - vL) - rhoR*(sR - vR));

    std::array<double , 6> uHLLCn_then_n1;

    double mL = rhoL*(sL - vL)/(sL - Sstar);
    double mR = rhoR*(sR - vR)/(sR - Sstar);

    uHLLCn_then_n1[0] = mL;
    uHLLCn_then_n1[1] = mL*Sstar;
    uHLLCn_then_n1[2] = mL * ( (x[2] / rhoL) + (Sstar - vL) * (Sstar + pL / (rhoL * (sL - vL))) );

    uHLLCn_then_n1[3] = mR;
    uHLLCn_then_n1[4] = mR*Sstar;
    uHLLCn_then_n1[5] = mR * ( (y[2] / rhoR) + (Sstar - vR) * (Sstar + pR / (rhoR * (sR - vR))) );


    return uHLLCn_then_n1;
}

std::array<double , 3> getFlux(std::array<double , 3> x, std::array<double , 3> y , double gamma){
    double rhoL = x[0];
    double rhoR = y[0];
    double vL = x[1] / x[0];
    double vR = y[1] / y[0];
    double pL = (gamma -1)*(x[2] - 0.5*vL*vL*rhoL);
    double pR = (gamma -1)*(y[2] - 0.5*vR*vR*rhoR);

    std::array<double , 6> uHLLCn_then_n1 = uHLLC(x , y, gamma);
    std::array<double , 3> uHLLCL;
    std::array<double , 3> uHLLCR;

    for(int i=0 ; i<=2 ; ++i){
        uHLLCL[i] = uHLLCn_then_n1[i];
        uHLLCR[i] = uHLLCn_then_n1[i+3];
    }

    std::array<double , 3> flux;
    
    double Splus = std::abs(vL) + std::sqrt(gamma *  pL / rhoL);
    if(std::abs(vR) + std::sqrt(gamma *  pR / rhoR) > Splus){
        Splus = std::abs(vR) + std::sqrt(gamma *  pR / rhoR);
    }

    double sL = -Splus;
    double sR = Splus;

    double Sstar = (pR - pL + rhoL*vL*(sL - vL) - rhoR*vR*(sR - vR)) / (rhoL*(sL - vL) - rhoR*(sR - vR));

    if(sL >= 0){
        flux = flux_def(x , gamma);
    }
    if(sL < 0 && Sstar >= 0){
        for(int j =0 ; j<=2; ++j){
            flux[j] = flux_def(x , gamma)[j] + sL*(uHLLCL[j] - x[j]);
        }
    }
    if(Sstar < 0 && sR >= 0){
        for(int j =0 ; j<=2; ++j){
            flux[j] = flux_def(y , gamma)[j] + sR*(uHLLCR[j] - y[j]);
        }   
    }
    if(sR <0){
        flux = flux_def(y,gamma);
    }

    return flux;
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
    double time;
    

    double dx = (x1 - x0) / nCells; //the space steps 

    // Initial conditions!
    for(int i = 0; i < u.size(); i++) {
        //x 0 is at point i=1/2
        double x = x0 + (i-0.5) * dx;
        std::array<double, 3> prim;
        if(x <= 0.25) {
            prim[0] = 1; // Density
            prim[1] = 0.75; // Velocity
            prim[2] = 1; // Pressure
            } else {
            prim[0] = 0.125; // Density
            prim[1] = 0; // Velocity
            prim[2] = 0.1; // Pressure
        }

        u[i] = primitiveToConservative(prim, gamma);
    }


    double dt = computeTimeStep(u , C , dx, gamma); //the time steps

    for(int counter =0 ; counter<=20; ++counter){
        double t = tStart;
        do {
            // Compute the stable time step for this iteration

            dt = computeTimeStep(u , C , dx, gamma); 
            t = t + dt;
            std::cout<<t<<dt<<std::endl;

            //You may want to manually reduce dt if this would overshoot tStop
            //Apply boundary conditions

            // Trasmissive boundary conditions
            u[0] = u[1];
            u[nCells + 1] = u[nCells];


            for(int i = 0; i < nCells+1; i++) { //Define the fluxes
                // flux[i] corresponds to cell i+1/2 
                flux[i] = getFlux( u[i], u[i+1] , gamma);
                //std::cout << flux[i][0] << flux[i][1] << flux[i][2]<< std::endl;
            }

            //the below has a i-1 flux which means we need to define a flux at 0 so make sure the above ^ starts at 0! this is because we have another edge with the number of cells (like the walls)

            for(int i = 1; i < nCells+1; i++) { //Update the data
                for(int j=0; j<=2; ++j){
                    uPlus1[i][j] = u[i][j] - (dt/dx) * (flux[i][j] - flux[i-1][j]);
                }
            }
        
            // Now replace u with the updated data for the next time step

            u = uPlus1;
        } while (t < tStop/(21-counter));
        std::cout << tStop/(21-counter) << std::endl;

        //still need to convert it back to primitive

        //define final results

        std::vector<std::array<double,3>> results(nCells+2);

        for(int i=0; i<= results.size() -1; ++i){
            results[i] = conservativeToPrimative(u[i], gamma);
        }


        //output
        std::string filename = "euler_" + std::to_string(counter) + ".dat";
        std::ofstream output(filename);
        for (int i = 1; i <= nCells; ++i) {
            double x = x0 + (i - 1) * dx;
            output << x << " " << results[i][0] <<  " " << results[i][1] <<  " " << results[i][2] << std::endl;
            //std::cout << x << " " << u[i][0] <<  " " << u[i][1] <<  " " << u[i][2] << std::endl;
        }
    }
}