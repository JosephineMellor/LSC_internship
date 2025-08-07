#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <array>
#include <filesystem>
#include "mutation++/mutation++.h"

//path for the Eigen and Mutation libraries: g++ -I /home/csc2301/projects/libraries/include -I /home/csc2301/projects/helloworld/internship_codes/eigen-3.4.0 -O3 SLIC2d_TeqS.cpp -o SLIC2d_TeqS -lmutation++

//export LD_LIBRARY_PATH=/home/csc2301/projects/libraries/lib:$LD_LIBRARY_PATH
//./SLIC2d_TeqS

//we will be looking at air here
Mutation::Mixture mix("air_5");
const double tolerance = 1e-5;
const int max = 100;

//the values need to be redimentionalised from reference
const double rho_ref = 1.225;     // kg/m³
const double p_ref   = 101325.0;  // Pa
const double u_ref   = 340.29;    // m/s (speed of sound at STP)


//start by converting between primitive and conservative variables

std::array<double , 4> PrimativeToConservative(std::array<double , 4> prim , double gamma){
    std::array<double , 4> consv;
    consv[0] = prim[0];
    consv[1] = prim[0] * prim[1];
    consv[2] = prim[0] * prim[2];
    consv[3] = prim[3] / (gamma -1) + ( 0.5 * prim[0] * (prim[1]*prim[1] + prim[2]*prim[2]));
    return consv;
}

std::array<double, 4> PrimativeToConservative2(std::array<double, 4> prim, Mutation::Mixture& mix) {
    std::array<double, 4> consv;
    consv[0] = prim[0];
    consv[1] = prim[0] * prim[1];
    consv[2] = prim[0] * prim[2];

    double rho = prim[0];
    double u = prim[1];
    double v = prim[2];
    double p = prim[3];

    int i = 0;
    double difference = 1.0;

    std::vector<double> species_densities(mix.nSpecies());

    double T = 300.0;

   while (std::abs(difference) > tolerance && i < max) {
        i++;

        mix.setState(species_densities.data(), &T, 1);

        difference = p - mix.P();

        // Relaxed update
        T += difference / 10.0;
        T = std::clamp(T, 100.0, 20000.0); // prevent bad values
    }

    double internalEnergy = mix.mixtureHMass() - p / rho;
    double totalEnergy = rho * (internalEnergy + 0.5 * (u*u + v*v));

    // Sanity checks
    if (std::isnan(totalEnergy) || totalEnergy <= 0.0) {
        std::cerr << " T=" << T << ", p=" << p<< ", rho=" << rho << ", internalEnergy=" << internalEnergy<< ", u=" << u << ", v=" << v << "Hmass = "<< mix.mixtureHMass()<<"\n";
        exit(1); 
    }

    consv[3] = totalEnergy;
    return consv;
}


std::array<double , 4> ConservativeToPrimative(std::array<double , 4> consv , Mutation::Mixture& mix){
    std::array < double , 4> prim;
    prim[0] = consv[0];
    prim[1] = consv[1] / consv[0];
    prim[2] = consv[2] / consv[0];

    double rho = prim[0];
    double u = prim[1];
    double v = prim[2];
    double E = consv[3];

    //isolate the internal energy

    double intenergy = (E / rho) - 0.5*(u*u + v*v);

    double T = 300.0;

    int i =0;
    double difference =1.0;
    std::vector<double> species_densities(mix.nSpecies());

    while(std::abs(difference)> tolerance && i < max ){
        i++;

        //find the values using mutation of the state with the guess temperature
        mix.setState(species_densities.data() , &T,1);

        //then use mutation to get the internal energy that would be here
        double internalEnergy = mix.mixtureHMass() - mix.P() / rho;

        //see how far away our guess was from the actual value
        difference = internalEnergy - intenergy;

        //then update the value of T accordingly

        T -= difference*100;
    }

    prim[3] = mix.P();

    return prim;
}

//define the flux function but now we split into two functions fluxX and fluxY

std::array <double, 4> fluxX_def(std::array <double ,4> x , double gamma){
    std::array <double , 4> fluxX;
    double rho = x[0];
    double v_x = x[1] / x[0];
    double v_y = x[2] / x[0];
    double KE = 0.5*rho*(v_x*v_x + v_y*v_y);
    double p = (gamma -1)*(x[3] - KE);

    fluxX[0] = x[1];
    fluxX[1] = rho*v_x*v_x + p;
    fluxX[2] = rho * v_x * v_y;
    fluxX[3] = (x[3] + p)*v_x;

    return fluxX;
} //flux for the x-direction

std::array <double, 4> fluxY_def(std::array <double ,4> x , double gamma){
    std::array <double , 4> fluxY;
    double rho = x[0];
    double v_x = x[1] / x[0];
    double v_y = x[2] / x[0];
    double KE = 0.5*rho*(v_x*v_x + v_y*v_y);
    double p = (gamma -1)*(x[3] - KE);

    fluxY[0] = x[2];
    fluxY[1] = rho*v_x*v_y;
    fluxY[2] = rho * v_y * v_y + p;
    fluxY[3] = (x[3] + p)*v_y;

    return fluxY;
} //flux for the y-direction

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

// Compute MinBee-limited slopes for reconstruction
std::array<double, 4> getxiDeltas(const std::array<double, 4>& u_left, const std::array<double, 4>& u_center,  const std::array<double, 4>& u_right) {
    std::array<double, 4> slopes;
    
    for (int k = 0; k < 4; k++) {
        double backward_diff = u_center[k] - u_left[k];
        double forward_diff = u_right[k] - u_center[k];
        slopes[k] = minbee(backward_diff, forward_diff);
    }
    
    return slopes;
}

// Reconstruct left and right states at cell interface using MinBee
void getUbar(const std::vector<std::array<double,4>>& u, int i,  std::array<double, 4>& uBarL,  std::array<double, 4>& uBarR) {
    
    // Compute slopes using MinBee 
    std::array<double, 4> xiDeltaL = getxiDeltas(u[i-1], u[i], u[i+1]);
    std::array<double, 4> xiDeltaR = getxiDeltas(u[i], u[i+1], u[i+2]);
    
    // get the ubars but just replace the others with them
    for (int k = 0; k < 3; k++) {
        uBarL[k] = u[i][k] + 0.5 * xiDeltaL[k];      
        uBarR[k] = u[i+1][k] - 0.5 * xiDeltaR[k];  
    }
}

std::array<double , 4> getFluxX(std::array<double , 4> x , std::array< double , 4> y , double dx , double dt , double gamma){
    //find the flux at u_{i} and u_{i+1}

    std::array<double , 4> f_1 = fluxX_def(x, gamma); //f(ui)
    std::array<double , 4> f_2 = fluxX_def(y, gamma); //f(ui+1)

    //find u_{i+1/2}

    std::array<double , 4> uPlusHalf;
    for(int i=0; i<=3; ++i){
        uPlusHalf[i] = 0.5*(x[i] + y[i]) - 0.5*(dt/dx)*(f_2[i] - f_1[i]);
    }

    std::array<double, 4> RI_flux = fluxX_def(uPlusHalf , gamma); //richtmyer flux
    std::array<double, 4> LF_flux; 
    std::array<double, 4> FORCE_flux;
    
    for(int i=0; i<=3; ++i){
        LF_flux[i] = (0.5*dx/dt)*(x[i] - y[i]) + 0.5*(f_1[i] + f_2[i]); // lax friedrichs flux
        FORCE_flux[i] = 0.5 * (LF_flux[i] + RI_flux[i]); //force flux
    }

    return FORCE_flux;

    
} //FORCE in x


std::array<double , 4> getFluxY(std::array<double , 4> x , std::array< double , 4> y , double dy , double dt , double gamma){
    //find the flux at u_{i} and u_{i+1}

    std::array<double , 4> f_1 = fluxY_def(x, gamma); //f(ui)
    std::array<double , 4> f_2 = fluxY_def(y, gamma); //f(ui+1)

    //find u_{i+1/2}

    std::array<double , 4> uPlusHalf;
    for(int i=0; i<=3; ++i){
        uPlusHalf[i] = 0.5*(x[i] + y[i]) - 0.5*(dt/dy)*(f_2[i] - f_1[i]);
    }

    std::array<double, 4> RI_flux = fluxY_def(uPlusHalf , gamma); //richtmyer flux
    std::array<double, 4> LF_flux; 
    std::array<double, 4> FORCE_flux;
    
    for(int i=0; i<=3; ++i){
        LF_flux[i] = (0.5*dy/dt)*(x[i] - y[i]) + 0.5*(f_1[i] + f_2[i]); // lax friedrichs flux
        FORCE_flux[i] = 0.5 * (LF_flux[i] + RI_flux[i]); //force flux
    }

    return FORCE_flux;
    
} //FORCE in y

double ComputeTimeStep(const std::vector<std::vector<std::array <double , 4>>>& u , double C , double dx , double dy ,double gamma){
    double dspace = dx;
    if(dy < dx){
        dspace = dy;
    }

    double maxSpeed = 0.0;

    for(int i=2 ; i<=u.size()-2 ; ++i){
        for(int j=2; j<=u[0].size()-2; ++j){
            double rho = u[i][j][0];
            double v_x = u[i][j][1] / rho;
            double v_y = u[i][j][2] / rho;
            double speed = std::sqrt(v_x*v_x + v_y*v_y);
            double KE = 0.5 * rho* speed * speed;
            double pressure = (gamma -1.0) * (u[i][j][3] - KE);

            double sound_speed = std::sqrt(gamma * pressure / rho);
            double new_speed = speed + sound_speed;

            if(new_speed > maxSpeed){
                maxSpeed = new_speed;
            }
            
        }
    }  
    
    double dt = C * dspace / maxSpeed;
    return dt;

} 

// Function to apply boundary conditions
void applyBoundaryConditions(std::vector<std::vector<std::array<double, 4>>>& u, int nxCells, int nyCells) {
    // ------ TRANSMISSIVE --------

    // Left and Right boundaries
    for (int j = 0; j < nyCells + 4; ++j) {
        u[0][j] = u[2][j];      // Left boundary: copy from first interior cell
        u[1][j] = u[2][j];      // Left ghost cell
        u[nxCells + 2][j] = u[nxCells + 1][j];  // Right ghost cell
        u[nxCells + 3][j] = u[nxCells + 1][j];  // Right boundary
    }

    // Bottom and Top boundaries
    for (int i = 0; i < nxCells + 4; ++i) {
        u[i][0] = u[i][2];      // Bottom boundary: copy from first interior cell
        u[i][1] = u[i][2];      // Bottom ghost cell
        u[i][nyCells + 2] = u[i][nyCells + 1];  // Top ghost cell
        u[i][nyCells + 3] = u[i][nyCells + 1];  // Top boundary
    }


    // -------- PERIODIC --------

    //left and right
    // for (int j = 0; j < nyCells + 2; ++j) {
    //     u[0][j] = u[nxCells][j]; 
    //     u[1][j] = u[nxCells+1][j];      // Left boundary
    //     u[nxCells + 2][j] = u[2][j]; 
    //     u[nxCells + 3][j] = u[3][j];  // Right boundary
    // }

    //bottom and top
    // for (int i = 0; i < nxCells + 2; ++i) {
    //     u[i][0] = u[i][nyCells];      // Bottom boundary
    //     u[i][1] = u[i][nyCells+1];    
    //     u[i][nyCells + 2] = u[i][2];  // Top boundary
    //     u[i][nyCells + 3] = u[i][3];  
    // }


    // ------ REFLECTIVE -------

    //left and right
    // for (int j = 0; j < nyCells + 2; ++j) {//reflect in u_y 
    //     u[0][j] = u[2][j];
    //     u[1][j] = u[3][j];        // Bottom boundary
    //     u[nyCells + 2][j] = u[nyCells][j];
    //     u[nyCells + 3][j] = u[nyCells + 1][j];  // Top boundary
    //     u[0][j][2] = -u[2][j][2]; 
    //     u[1][j][2] = -u[3][j][2]; 
    //     u[nyCells + 2][j][2] = -u[nyCells][j][2]; 
    //     u[nyCells + 3][j][2] = -u[nyCells + 1][j][2]; 
    // }

    // Bottom and Top boundaries
    // for (int i = 0; i < nxCells + 2; ++i) {//reflect in u_y
    //     u[i][0] = u[i][2];
    //     u[i][1] = u[i][3];        // Bottom boundary
    //     u[i][nyCells + 2] = u[i][nyCells];
    //     u[i][nyCells + 3] = u[i][nyCells + 1];  // Top boundary
    //     u[i][0][2] = -u[i][2][2]; 
    //     u[i][1][2] = -u[i][3][2]; 
    //     u[i][nyCells + 2][2] = -u[i][nyCells][2]; 
    //     u[i][nyCells + 3][2] = -u[i][nyCells + 1][2]; 
    // }
}



std::vector<std::vector<std::array<double, 4> > > XthenY(std::vector<std::vector<std::array<double, 4> > > u , double dx , double dy , double dt , int nxCells , int nyCells , double x0 ,double x1 , double y0 ,double y1 ,double tStart , double tStop , double C,double gamma,double omega ){
    applyBoundaryConditions(u, nxCells , nyCells);
    std::vector<std::vector<std::array<double, 4> > > uPlus1;
    uPlus1.resize(nxCells+4, std::vector<std::array<double, 4> >(nyCells + 4)); //set up uPlus1

    std::vector<std::vector<std::array<double, 4> > > fluxX;
    fluxX.resize(nxCells+4, std::vector<std::array<double, 4> >(nyCells + 4)); //fluxX

    std::vector<std::vector<std::array<double, 4> > > fluxY;
    fluxY.resize(nxCells+4, std::vector<std::array<double, 4> >(nyCells + 4)); //fluxY

    std::vector<std::vector<std::array<double, 4> > > uPlus1Bar;
    uPlus1Bar.resize(nxCells+4, std::vector<std::array<double, 4> >(nyCells + 4)); //uPlus1Bar

    std::vector<std::vector<std::array<double, 4> > > uBarL;
    uBarL.resize(nxCells+4, std::vector<std::array<double, 4> >(nyCells + 4));

    std::vector<std::vector<std::array<double, 4> > > uBarR;
    uBarR.resize(nxCells+4, std::vector<std::array<double, 4> >(nyCells + 4));

    std::vector<std::vector<std::array<double, 4> > > uBarHalfL;
    uBarHalfL.resize(nxCells+4, std::vector<std::array<double, 4> >(nyCells + 4));

    std::vector<std::vector<std::array<double, 4> > > uBarHalfR;
    uBarHalfR.resize(nxCells+4, std::vector<std::array<double, 4> >(nyCells + 4));

    

        //calculate the flux in the x-direction
        

        //start with uBar
        for(int j = 2; j < nyCells+2; j++) {
            for(int i = 2; i < nxCells+2; i++) {
                for(int k=0 ; k<4 ; ++k){
                    double DeltaPlus = u[i+1][j][k] - u[i][j][k];
                    double DeltaMinus = u[i][j][k] - u[i-1][j][k];
                    double r = DeltaMinus / (DeltaPlus  + 1e-8);
                    
                    double xi_L = 2.0*r/(1+r);
                    double xi_R = 2.0/(1+r);
                    double xi;
                    double Delta = 0.5*(1+omega)*DeltaMinus + 0.5*(1-omega)*DeltaPlus;
                    
                    if(r<=0){ xi=0;}
                    else if(r>0 && r<=1){ xi=r;}
                    else{ 
                        xi=std::fmin(1, xi_R);
                        
                    }

                    uBarL[i][j][k] = u[i][j][k] - 0.5 * xi * Delta;
                    uBarR[i][j][k] = u[i][j][k] + 0.5 * xi * Delta;
                    
                }

            }
        }
        //then uBarHalf
        for(int j = 2; j < nyCells+2; j++) {
            for(int i = 2; i < nxCells+2; i++) {
                for(int k=0 ; k<4 ; ++k){
                    uBarHalfL[i][j][k] = uBarL[i][j][k] - 0.5*(dt/dx)*(fluxX_def(uBarR[i][j] , gamma)[k]-fluxX_def(uBarL[i][j] , gamma)[k]);
                    uBarHalfR[i][j][k] = uBarR[i][j][k] - 0.5*(dt/dx)*(fluxX_def(uBarR[i][j] , gamma)[k]-fluxX_def(uBarL[i][j] , gamma)[k]);
                }
            }
        }
        applyBoundaryConditions(uBarHalfL, nxCells , nyCells);
        applyBoundaryConditions(uBarHalfR, nxCells , nyCells);
        //now we have flux in x direction
        for (int j = 1; j < nyCells + 3; j++) {
            for (int i = 1; i < nxCells + 3; i++) { 
                fluxX[i][j] = getFluxX(uBarHalfR[i][j], uBarHalfL[i+1][j], dx, dt, gamma);
            }
        }
        //from this we get our intermediate ubar
        for(int j = 2; j < nyCells+2; j++) {
            for(int i = 2; i < nxCells+2; i++) {
                for(int v = 0; v <= 3; v++) {
                    uPlus1Bar[i][j][v] = u[i][j][v] - (dt/dx)*(fluxX[i][j][v]-fluxX[i-1][j][v]);
                }
            }
        }
        // apply boundary conditions to uPlus1Bar
        applyBoundaryConditions(uPlus1Bar , nxCells , nyCells);


        //now move in the y direction


        //start with uBar
        for(int j = 2; j < nyCells+2; j++) {
            for(int i = 2; i < nxCells+2; i++) {
                for(int k=0 ; k<4 ; ++k){
                    double DeltaPlus = uPlus1Bar[i][j+1][k] - uPlus1Bar[i][j][k];
                    double DeltaMinus = uPlus1Bar[i][j][k] - uPlus1Bar[i][j-1][k];
                    double r = DeltaMinus / (DeltaPlus  + 1e-8);
                    
                    double xi_L = 2.0*r/(1+r);
                    double xi_R = 2.0/(1+r);
                    double xi;
                    double Delta = 0.5*(1+omega)*DeltaMinus + 0.5*(1-omega)*DeltaPlus;
                    
                    if(r<=0){ xi=0;}
                    else if(r>0 && r<=1){ xi=r;}
                    else{ 
                        xi=std::fmin(1, xi_R);
                        
                    }

                    uBarL[i][j][k] = uPlus1Bar[i][j][k] - 0.5 * xi * Delta;  
                    uBarR[i][j][k] = uPlus1Bar[i][j][k] + 0.5 * xi * Delta;  
                    
                }

            }
        }
        //then uBarHalf
        for(int j = 2; j < nyCells+2; j++) {
            for(int i = 2; i < nxCells+2; i++) {
                for(int k=0 ; k<4 ; ++k){
                    uBarHalfL[i][j][k] = uBarL[i][j][k] - 0.5*(dt/dy)*(fluxY_def(uBarR[i][j] , gamma)[k]-fluxY_def(uBarL[i][j] , gamma)[k]);
                    uBarHalfR[i][j][k] = uBarR[i][j][k] - 0.5*(dt/dy)*(fluxY_def(uBarR[i][j] , gamma)[k]-fluxY_def(uBarL[i][j] , gamma)[k]);
                }
            }
        }
        applyBoundaryConditions(uBarHalfL, nxCells , nyCells);
        applyBoundaryConditions(uBarHalfR, nxCells , nyCells);
        //now we have flux in y direction
        for (int j = 1; j < nyCells + 3; j++) {
            for (int i = 1; i < nxCells + 3; i++) { 
                fluxY[i][j] = getFluxY(uBarHalfR[i][j], uBarHalfL[i][j+1], dy, dt, gamma);
            }
        }


        //update the results!
        for(int j = 2; j < nyCells+2; j++) {
            for(int i = 2; i < nxCells+2; i++) {
                for(int v = 0; v <= 3; v++) {
                    uPlus1[i][j][v] = uPlus1Bar[i][j][v] - (dt/dy)*(fluxY[i][j][v]-fluxY[i][j-1][v]);
                }
            }
        }
    
    return uPlus1;
}


std::vector<std::vector<std::array<double, 4> > > YthenX(std::vector<std::vector<std::array<double, 4> > > u , double dx , double dy , double dt , int nxCells , int nyCells , double x0 ,double x1 , double y0 ,double y1 ,double tStart , double tStop , double C,double gamma,double omega ){
    applyBoundaryConditions(u, nxCells , nyCells);
    std::vector<std::vector<std::array<double, 4> > > uPlus1;
    uPlus1.resize(nxCells+4, std::vector<std::array<double, 4> >(nyCells + 4)); //set up uPlus1

    std::vector<std::vector<std::array<double, 4> > > fluxX;
    fluxX.resize(nxCells+4, std::vector<std::array<double, 4> >(nyCells + 4)); //fluxX

    std::vector<std::vector<std::array<double, 4> > > fluxY;
    fluxY.resize(nxCells+4, std::vector<std::array<double, 4> >(nyCells + 4)); //fluxY

    std::vector<std::vector<std::array<double, 4> > > uPlus1Bar;
    uPlus1Bar.resize(nxCells+4, std::vector<std::array<double, 4> >(nyCells + 4)); //uPlus1Bar

    std::vector<std::vector<std::array<double, 4> > > uBarL;
    uBarL.resize(nxCells+4, std::vector<std::array<double, 4> >(nyCells + 4));

    std::vector<std::vector<std::array<double, 4> > > uBarR;
    uBarR.resize(nxCells+4, std::vector<std::array<double, 4> >(nyCells + 4));

    std::vector<std::vector<std::array<double, 4> > > uBarHalfL;
    uBarHalfL.resize(nxCells+4, std::vector<std::array<double, 4> >(nyCells + 4));

    std::vector<std::vector<std::array<double, 4> > > uBarHalfR;
    uBarHalfR.resize(nxCells+4, std::vector<std::array<double, 4> >(nyCells + 4));

    

        //calculate the flux in the x-direction
        

        //start with uBar
        for(int j = 2; j < nyCells+2; j++) {
            for(int i = 2; i < nxCells+2; i++) {
                for(int k=0 ; k<4 ; ++k){
                    double DeltaPlus = u[i][j+1][k] - u[i][j][k];
                    double DeltaMinus = u[i][j][k] - u[i][j-1][k];
                    double r = DeltaMinus / (DeltaPlus  + 1e-8);
                    
                    double xi_L = 2.0*r/(1+r);
                    double xi_R = 2.0/(1+r);
                    double xi;
                    double Delta = 0.5*(1+omega)*DeltaMinus + 0.5*(1-omega)*DeltaPlus;
                    
                    if(r<=0){ xi=0;}
                    else if(r>0 && r<=1){ xi=r;}
                    else{ 
                        xi=std::fmin(1, xi_R);
                        
                    }

                    uBarL[i][j][k] = u[i][j][k] - 0.5 * xi * Delta;
                    uBarR[i][j][k] = u[i][j][k] + 0.5 * xi * Delta;
                    
                }

            }
        }
        //then uBarHalf
        for(int j = 2; j < nyCells+2; j++) {
            for(int i = 2; i < nxCells+2; i++) {
                for(int k=0 ; k<4 ; ++k){
                    uBarHalfL[i][j][k] = uBarL[i][j][k] - 0.5*(dt/dy)*(fluxY_def(uBarR[i][j] , gamma)[k]-fluxY_def(uBarL[i][j] , gamma)[k]);
                    uBarHalfR[i][j][k] = uBarR[i][j][k] - 0.5*(dt/dy)*(fluxY_def(uBarR[i][j] , gamma)[k]-fluxY_def(uBarL[i][j] , gamma)[k]);
                }
            }
        }
        applyBoundaryConditions(uBarHalfL, nxCells , nyCells);
        applyBoundaryConditions(uBarHalfR, nxCells , nyCells);
        //now we have flux in y direction
        for (int j = 1; j < nyCells + 3; j++) {
            for (int i = 1; i < nxCells + 3; i++) { 
                fluxY[i][j] = getFluxY(uBarHalfR[i][j], uBarHalfL[i][j+1], dy, dt, gamma);
            }
        }
        //from this we get our intermediate ubar
        for(int j = 2; j < nyCells+2; j++) {
            for(int i = 2; i < nxCells+2; i++) {
                for(int v = 0; v <= 3; v++) {
                    uPlus1Bar[i][j][v] = u[i][j][v] - (dt/dy)*(fluxY[i][j][v]-fluxY[i][j-1][v]);
                }
            }
        }
        // apply boundary conditions to uPlus1Bar
        applyBoundaryConditions(uPlus1Bar , nxCells , nyCells);


        //now move in the x direction


        //start with uBar
        for(int j = 2; j < nyCells+2; j++) {
            for(int i = 2; i < nxCells+2; i++) {
                for(int k=0 ; k<4 ; ++k){
                    double DeltaPlus = uPlus1Bar[i+1][j][k] - uPlus1Bar[i][j][k];
                    double DeltaMinus = uPlus1Bar[i][j][k] - uPlus1Bar[i-1][j][k];
                    double r = DeltaMinus / (DeltaPlus  + 1e-8);
                    
                    double xi_L = 2.0*r/(1+r);
                    double xi_R = 2.0/(1+r);
                    double xi;
                    double Delta = 0.5*(1+omega)*DeltaMinus + 0.5*(1-omega)*DeltaPlus;
                    
                    if(r<=0){ xi=0;}
                    else if(r>0 && r<=1){ xi=r;}
                    else{ 
                        xi=std::fmin(1, xi_R);
                        
                    }

                    uBarL[i][j][k] = uPlus1Bar[i][j][k] - 0.5 * xi * Delta;  
                    uBarR[i][j][k] = uPlus1Bar[i][j][k] + 0.5 * xi * Delta;  
                    
                }

            }
        }
        //then uBarHalf
        for(int j = 2; j < nyCells+2; j++) {
            for(int i = 2; i < nxCells+2; i++) {
                for(int k=0 ; k<4 ; ++k){
                    uBarHalfL[i][j][k] = uBarL[i][j][k] - 0.5*(dt/dx)*(fluxX_def(uBarR[i][j] , gamma)[k]-fluxX_def(uBarL[i][j] , gamma)[k]);
                    uBarHalfR[i][j][k] = uBarR[i][j][k] - 0.5*(dt/dx)*(fluxX_def(uBarR[i][j] , gamma)[k]-fluxX_def(uBarL[i][j] , gamma)[k]);
                }
            }
        }
        applyBoundaryConditions(uBarHalfL, nxCells , nyCells);
        applyBoundaryConditions(uBarHalfR, nxCells , nyCells);
        //now we have flux in y direction
        for (int j = 1; j < nyCells + 3; j++) {
            for (int i = 1; i < nxCells + 3; i++) { 
                fluxY[i][j] = getFluxX(uBarHalfR[i][j], uBarHalfL[i+1][j], dx, dt, gamma);
            }
        }


        //update the results!
        for(int j = 2; j < nyCells+2; j++) {
            for(int i = 2; i < nxCells+2; i++) {
                for(int v = 0; v <= 3; v++) {
                    uPlus1[i][j][v] = uPlus1Bar[i][j][v] - (dt/dx)*(fluxX[i][j][v]-fluxX[i-1][j][v]);
                }
            }
        }
    
    return uPlus1;
}


int main(){

    Mutation::Mixture mix("air_5");

    int nxCells = 100;
    int nyCells = 100;
    double x0 = 0.0;
    double x1 = 1.0;
    double y0 = 0.0;
    double y1 = 1.0;
    double tStart = 0.0;
    double tStop = 0.5/u_ref;
    double C = 0.8;
    double gamma = 1.4;
    double omega =0;

    std::vector<std::vector<std::array<double, 4> > > u;
    u.resize(nxCells+4, std::vector<std::array<double, 4> >(nyCells + 4)); //set up u
    

    double dx = (x1 - x0) / nxCells;
    double dy = (y1 - y0) / nyCells;
    double dt = ComputeTimeStep(u,C,dx,dy,gamma);

    //intial conditions!
    for(int i = 0; i < u.size(); i++) { 
        for(int j = 0; j < u[0].size(); j++) {
            std::array<double,4> prim;
            // coordinates
            double x = x0 + (i - 1.5) * dx;
            double y = y0 + (j - 1.5) * dy;

            
            if ( x>0.5 && y>0.5) {
                prim[0] = 1.5;   // density (rho)
                prim[1] = 0;   // x-velocity (u)
                prim[2] = 0;   // y-velocity (v)
                prim[3] = 1.5;   // pressure (p)
            } 
            else if(x<=0.5 && y > 0.5) {
                prim[0] = 0.5323;
                prim[1] = 1.206;
                prim[2] = 0;
                prim[3] = 0.3;
            }
            else if(x<=0.5 && y<=0.5) {
                prim[0] = 0.138;
                prim[1] = 1.206;
                prim[2] = 1.206;
                prim[3] = 0.029;
            }
            else{
                prim[0] = 0.5323;
                prim[1] = 0;
                prim[2] = 1.206;
                prim[3] = 0.3;
            }

            //dimensionalise
            prim[0] *= rho_ref;     // kg/m³
            prim[1] *= u_ref;       // m/s
            prim[2] *= u_ref;       // m/s
            prim[3] *= p_ref;       // Pa

            u[i][j] = PrimativeToConservative2(prim, mix);
        }
    }

    applyBoundaryConditions(u , nxCells , nyCells);

    double t = tStart;

    do{
        dt = ComputeTimeStep(u,C,dx,dy,gamma);
        t +=dt;

        std::cout << "t = "<< t<<" dt = "<< dt<< std::endl; 
        applyBoundaryConditions(u , nxCells , nyCells);

        std::vector<std::vector<std::array<double, 4> > > uPlus1x = XthenY(u ,  dx ,  dy ,  dt ,  nxCells ,  nyCells ,  x0 , x1 ,  y0 , y1 , tStart ,  tStop ,  C, gamma, omega);
        std::vector<std::vector<std::array<double, 4> > > uPlus1y = YthenX(u ,  dx ,  dy ,  dt ,  nxCells ,  nyCells ,  x0 , x1 ,  y0 , y1 , tStart ,  tStop ,  C, gamma, omega);

        for(int j = 2; j < nyCells+2; j++) {
            for(int i = 2; i < nxCells+2; i++) {
                for(int k=0 ; k<4 ; ++k){
                    u[i][j][k] = 0.5*(uPlus1x[i][j][k] + uPlus1y[i][j][k]);
                }
            }
        }
        applyBoundaryConditions(u , nxCells , nyCells);
    }while(t<tStop);

    //now convert back to primitive
    std::vector<std::vector<std::array<double, 4> > > results;
    results.resize(nxCells+2, std::vector<std::array<double, 4> >(nyCells + 2)); //results
 
    for(int j = 1; j < nyCells+2; j++) { 
        for(int i = 1; i < nxCells+2; i++) {
            results[i][j] = ConservativeToPrimative(u[i][j], mix);
        }
    }


    //output the results
    std::ofstream output("SLIC2d.dat");
    for(int j = 1; j < nyCells+2; j++) { 
        for(int i = 1; i < nxCells+2; i++) {
            double x = x0 + (i - 1)*dx;
            double y = y0 + (j - 1)*dy;
            output << x << " " << y << " " << results[i][j][0] << " " << results[i][j][1] << " " << results[i][j][2] << " " << results[i][j][3] << std::endl;
        }
        output<<std::endl;
    }

}