#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <array>
#include <filesystem>

//SLIC is just FORCE but you have to put u to ubar beforehand


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

//using the minbee limiter
double minbee(double deltaMinus , double deltaPlus , double omega = 0.0){
    if (deltaPlus == 0.0) {
        return 0.0;  
    }
    double r = deltaMinus / deltaPlus;
    double xi;
    if(r <=0){
        double xi =0.0;
    }
    else if(r <=1.0){
        xi = r;
    }
    else{
        xi = std::min(1.0, 2.0/(1.0+r));
    }
    
    double Delta = 0.5 * (1.0 + omega) * deltaMinus + 0.5 * (1.0 - omega) * deltaPlus;
    
    return xi * Delta;
}

std::array<double, 3> getxiDeltas(const std::array<double, 3>& u_left, const std::array<double, 3>& u_center,  const std::array<double, 3>& u_right) {
    std::array<double, 3> slopes;
    
    for (int k = 0; k < 3; k++) {
        double backward_diff = u_center[k] - u_left[k];
        double forward_diff = u_right[k] - u_center[k];
        slopes[k] = minbee(backward_diff, forward_diff);
    }
    
    return slopes;
}

// this function takes in u and spits out uBarL and uBarR
void getUbar(const std::vector<std::array<double,3>>& u, int i,  std::array<double, 3>& uBarL,  std::array<double, 3>& uBarR) {
    
    // Compute slopes using MinBee 
    std::array<double, 3> xiDeltaL = getxiDeltas(u[i-1], u[i], u[i+1]);
    std::array<double, 3> xiDeltaR = getxiDeltas(u[i], u[i+1], u[i+2]);
    
    // get the ubars but just replace the others with them
    for (int k = 0; k < 3; k++) {
        uBarL[k] = u[i][k] + 0.5 * xiDeltaL[k];      
        uBarR[k] = u[i+1][k] - 0.5 * xiDeltaR[k];  
    }
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
    double tStop = 0.25;
    double C = 0.8;
    double gamma = 1.4;
    double omega = 0;

    // Allocate matrices with 2 extra points for transmissive BCs
    std::vector<std::array<double,3>> u(nCells+2);
    std::vector<std::array<double,3>> uPlus1(nCells+2);
    std::vector<std::array<double,3>> flux(nCells+2);
    std::vector<std::array<double,3>> uBarL(nCells+2);
    std::vector<std::array<double,3>> uBarR(nCells+2);
    std::vector<std::array<double,3>> uBarHalfL(nCells+2);
    std::vector<std::array<double,3>> uBarHalfR(nCells+2);
    double dx = (x1 - x0) / nCells; //the space steps 
    double time;

    // Initial conditions!

    for(int i = 0; i < u.size(); i++) {
        // x 0 is at point i=1/2
        double x = x0 + (i-0.5) * dx;
        std::array<double, 3> prim;
        if(x <= 0.25) {
            prim[0] = 1; // Density
            prim[1] = 0; // Velocity
            prim[2] = 1; // Pressure
            } else {
            prim[0] = 0.125; // Density
            prim[1] = 0; // Velocity
            prim[2] = 0.1; // Pressure
        }

        u[i] = primitiveToConservative(prim, gamma);
    }

    
    



    double dt = computeTimeStep(u , C , dx, gamma); //the time steps
    for(int counter =1; counter<=20; ++counter){
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

            //find ubar


            for(int i=1; i<=nCells; ++i){
                for(int j=0; j<=2; ++j){
                    double DeltaPlus = u[i+1][j] - u[i][j];
                    double DeltaMinus = u[i][j] - u[i-1][j];
                    // std::cout<< DeltaMinus << " " << DeltaPlus << std::endl;
                    double r = DeltaMinus / (DeltaPlus  + 1e-8);
                    
                    double xi_L = 2.0*r/(1+r);
                    double xi_R = 2.0/(1+r);
                    double xi;
                    double Delta = 0.5*(1+omega)*DeltaMinus + 0.5*(1-omega)*DeltaPlus;
                    // std::cout << Delta << std::endl;
                
                    if(r<=0){ xi=0;}
                    else if(r>0 && r<=1){ xi=r;}
                    else{ 
                        xi=std::fmin(1, xi_R);
                        
                    }

                    // xi = 0.2;
                    // std::cout << xi << " " << xi_R << " " << r <<  std::endl;
                    uBarL[i][j] = u[i][j] - 0.5 * xi * Delta;
                    uBarR[i][j] = u[i][j] + 0.5 * xi * Delta;
                    
                    
                }
            }

            for(int i=1; i<=nCells; ++i){
                for(int j=0; j<=2; ++j){
                    uBarHalfL[i][j] = uBarL[i][j] - 0.5*(dt/dx)*(flux_def(uBarR[i] , gamma)[j]-flux_def(uBarL[i] , gamma)[j]);
                    uBarHalfR[i][j] = uBarR[i][j] - 0.5*(dt/dx)*(flux_def(uBarR[i] , gamma)[j]-flux_def(uBarL[i] , gamma)[j]);
                    
                }
            }

            

            uBarHalfL[0] = uBarHalfL[1];
            uBarHalfL[nCells + 1] = uBarHalfL[nCells];

            uBarHalfR[0] = uBarHalfR[1];
            uBarHalfR[nCells + 1] = uBarHalfR[nCells];


            for(int i = 0; i < nCells+1; i++) { //Define the fluxes
                // flux[i] corresponds to cell i+1/2 
                flux[i] = getFlux( uBarHalfR[i], uBarHalfL[i+1] , dx , dt, gamma);
            }

            //the below has a i-1 flux which means we need to define a flux at 0 so make sure the above ^ starts at 0! this is because we have another edge with the number of cells (like the walls)

            for(int i = 1; i <= nCells+1; i++) { //Update the data
                for(int j=0; j<=2; ++j){
                    uPlus1[i][j] = u[i][j] - (dt/dx) * (flux[i][j] - flux[i-1][j]);
                }
            }
        
            // Now replace u with the updated data for the next time step

            
            u = uPlus1;
        } while (t < tStop/(21-counter));
        time = tStop/(21-counter);
        std::cout <<time << std::endl;

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
            // std::cout << x << " " << u[i][0] <<  " " << u[i][1] <<  " " << u[i][2] << std::endl;
        }
    }
}