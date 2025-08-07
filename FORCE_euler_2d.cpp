#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <array>
#include <filesystem>

//start by converting between primitive and conservative variables

std::array<double , 4> PrimativeToConservative(std::array<double , 4> prim , double gamma){
    std::array<double , 4> consv;
    consv[0] = prim[0];
    consv[1] = prim[0] * prim[1];
    consv[2] = prim[0] * prim[2];
    consv[3] = prim[3] / (gamma -1) + ( 0.5 * prim[0] * (prim[1]*prim[1] + prim[2]*prim[2]));
    return consv;
}

std::array<double , 4> ConservativeToPrimative(std::array<double , 4> consv , double gamma){
    std::array < double , 4> prim;
    prim[0] = consv[0];
    prim[1] = consv[1] / consv[0];
    prim[2] = consv[2] / consv[0];
    prim[3] = (gamma-1)* (consv[3] - 0.5*prim[0]*(prim[1]*prim[1] + prim[2]*prim[2]));
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

//then need to find the FORCE flux in the x-direction by using the same method as in 1D

//here, we're inputting x as our u_{i} and y as u_{i+1}

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


//next we need to calculate dt but taking the max over i AND j so doing two for loops

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
}






int main(){
    int nxCells = 100;
    int nyCells = 100;
    double x0 = 0.0;
    double x1 = 1.0;
    double y0 = 0.0;
    double y1 = 1.0;
    double tStart = 0.0;
    double tStop = 0.5;
    double C = 0.3;
    double gamma = 1.4;

    std::vector<std::vector<std::array<double, 4> > > u;
    u.resize(nxCells+4, std::vector<std::array<double, 4> >(nyCells + 4)); //set up u

    std::vector<std::vector<std::array<double, 4> > > uPlus1;
    uPlus1.resize(nxCells+4, std::vector<std::array<double, 4> >(nyCells + 4)); //set up uPlus1

    std::vector<std::vector<std::array<double, 4> > > fluxX;
    fluxX.resize(nxCells+4, std::vector<std::array<double, 4> >(nyCells + 4)); //fluxX

    std::vector<std::vector<std::array<double, 4> > > fluxY;
    fluxY.resize(nxCells+4, std::vector<std::array<double, 4> >(nyCells + 4)); //fluxY

    std::vector<std::vector<std::array<double, 4> > > uPlus1Bar;
    uPlus1Bar.resize(nxCells+4, std::vector<std::array<double, 4> >(nyCells + 4)); //uPlus1Bar

    double dx = (x1-x0) / nxCells;
    double dy = (y1-y0) / nyCells; //space steps
    double time;

    //intial conditions!
    for(int i = 0; i < u.size(); i++) { 
        for(int j = 0; j < u[0].size(); j++) {
            std::array<double,4> prim;
            // coordinates
            double x = x0 + (i - 1.5) * dx;
            double y = y0 + (j - 1.5) * dy;

            
            if ( x<=0.5 && y<=0.5) {
                prim[0] = 1.0;   // density (rho)
                prim[1] = 0;   // x-velocity (u)
                prim[2] = 0;   // y-velocity (v)
                prim[3] = 1.0;   // pressure (p)
            } 
            else if(x>0.5 && y<= 0.5) {
                prim[0] = 2.0;
                prim[1] = 0;
                prim[2] = 0;
                prim[3] = 1.0;
            }
            else if(x<=0.5 && y>0.5) {
                prim[0] = 1.0;
                prim[1] = 0;
                prim[2] = 0;
                prim[3] = 2.0;
            }
            else{
                prim[0] = 3.0;
                prim[1] = 0;
                prim[2] = 0;
                prim[3] = 1.0;
            }

            u[i][j] = PrimativeToConservative(prim, gamma);
        }
    }


    double dt = ComputeTimeStep(u,C,dx,dy,gamma);

    //apply transmissive boundary conditions
    applyBoundaryConditions(u , nxCells , nyCells);


    //time loop is the same


    double t = tStart;
    do{

        //compute the time step
        dt = ComputeTimeStep(u,C,dx,dy,gamma);
        t = t+dt;
        std::cout << "t = " << t << ", dt = " << dt << std::endl;

        //transmissive boundary conditions along the edges of the mesh
        //waves are moving through the bounderies

        applyBoundaryConditions(u, nxCells , nyCells);

        //calculate the flux in the x-direction
        for (int j = 1; j < nyCells + 3; j++) {
            for (int i = 1; i < nxCells + 3; i++) { 
                fluxX[i][j] = getFluxX(u[i][j], u[i+1][j], dx, dt, gamma);
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

        // Apply boundary conditions to uPlus1Bar
        applyBoundaryConditions(uPlus1Bar , nxCells , nyCells);

        // //now we can move in the y-direction
        // //calculate the flux in the y-direction
        for (int j = 1; j < nyCells + 3; j++) { 
            for (int i = 2; i < nxCells + 2; i++) {
                fluxY[i][j] = getFluxY(uPlus1Bar[i][j], uPlus1Bar[i][j+1], dy, dt, gamma);
            }
        }


        //now we can update the data
        for(int j = 2; j < nyCells+2; j++) {
            for(int i = 2; i < nxCells+2; i++) {
                for(int v = 0; v <= 3; v++) {
                    uPlus1[i][j][v] = uPlus1Bar[i][j][v] - (dt/dy)*(fluxY[i][j][v]-fluxY[i][j-1][v]);
                }
            }
        }

        


        //replace with updated data
        u = uPlus1;
        applyBoundaryConditions(u, nxCells, nyCells);

    }while(t<tStop);


    //now convert back to primitive

    std::vector<std::vector<std::array<double, 4> > > results;
    results.resize(nxCells+2, std::vector<std::array<double, 4> >(nyCells + 2)); //results
 
    for(int j = 1; j < nyCells+2; j++) { 
        for(int i = 1; i < nxCells+2; i++) {
            results[i][j] = ConservativeToPrimative(u[i][j], gamma);
        }
    }

    //output the results

    std::ofstream output("euler_2d.dat");
    for(int j = 1; j < nyCells+2; j++) { 
        for(int i = 1; i < nxCells+2; i++) {
            double x = x0 + (i - 1)*dx;
            double y = y0 + (j - 1)*dy;
            output << x << " " << y << " " << results[i][j][0] << " " << results[i][j][1] << " " << results[i][j][2] << " " << results[i][j][3] << std::endl;
            //std::cout << x << " " << y << " " << results[i][j][0] << " " << results[i][j][1] << " " << results[i][j][2] << " " << results[i][j][3] << std::endl;
        }
    }

}