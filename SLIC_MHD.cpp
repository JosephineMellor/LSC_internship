#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <array>
#include <tuple>
#include <filesystem>

//Primitive to Conservative
std::array<double , 9> PrimitiveToConservative(const std::array<double , 9>& u , double gamma){
    std::array<double ,9> v;
    v[0] = u[0]; //density
    v[1] = u[0] * u[1]; //vx
    v[2] = u[0] * u[2]; //vy
    v[3] = u[0] * u[3]; //vz
    v[5] = u[5]; //Bx
    v[6] = u[6]; //By
    v[7] = u[7]; //Bz
    v[4] = u[4] / (gamma-1) + 0.5*u[0]*(u[1]*u[1] + u[2]*u[2] + u[3]*u[3]) + 0.5*(u[5]*u[5] + u[6]*u[6] + u[7]*u[7]); //energy
    v[8] = u[8];
    return v;
}
std::array<double, 9> ConservativeToPrimitive(const std::array<double , 9>& u , double gamma){
    std::array<double , 9> v;
    v[0] = u[0];
    v[1] = u[1] / u[0];
    v[2] = u[2] / u[0];
    v[3] = u[3] / u[0];
    v[5] = u[5];
    v[6] = u[6];
    v[7] = u[7];
    v[4] = (u[4] - 0.5*v[0]*(v[1]*v[1] + v[2]*v[2] + v[3]*v[3]) - 0.5*(u[5]*u[5] + u[6]*u[6] + u[7]*u[7]))*(gamma -1);//pressure
    v[8] = u[8];
    return v;
}
std::array<double , 9> FluxDefX(const std::array<double , 9>& u , double gamma ){//conservative variables go into this function
    std::array<double , 9> v = ConservativeToPrimitive(u , gamma);
    std::array<double , 9> flux;
    flux[0] = u[1]; //rho v_x
    flux[1] = u[0]*v[1]*v[1] + v[4] + 0.5*( u[5]*u[5] + u[6]*u[6] + u[7]*u[7] ) - u[5]*u[5]; 
    flux[2] = u[0]*v[1]*v[2] - u[5]*u[6];
    flux[3] = u[0]*v[1]*v[3] - u[5]*u[7];
    flux[4] = (u[4] + v[4] + 0.5*( u[5]*u[5] + u[6]*u[6] + u[7]*u[7] ))*v[1] - (v[1]*v[5]+v[2]*v[6]+v[3]*v[7])*v[5];
    flux[5] = 0;
    flux[6] = v[6]*v[1]-v[5]*v[2];
    flux[7] = v[7]*v[1]-v[5]*v[3];
    flux[8] = 0;
    return flux;
}

std::array<double , 9> FluxDefY(const std::array<double , 9>& u , double gamma ){//conservative variables go into this function
    std::array<double , 9> v = ConservativeToPrimitive(u , gamma);
    std::array<double , 9> flux;
    flux[0] = u[2]; //rho v_y
    flux[1] = u[0]*v[1]*v[2] - v[5]*v[6];
    flux[2] = u[0]*v[2]*v[2]+ v[4] + 0.5*( u[5]*u[5] + u[6]*u[6] + u[7]*u[7] ) - u[6]*u[6];
    flux[3] = u[0]*v[2]*v[3] - u[6]*u[7];
    flux[4] = (u[4] + v[4] + 0.5*( u[5]*u[5] + u[6]*u[6] + u[7]*u[7] ))*v[2] - (v[1]*v[5]+v[2]*v[6]+v[3]*v[7])*v[6];
    flux[5] = v[5]*v[2] - v[6]*v[1];
    flux[6] = 0;
    flux[7] = v[7]*v[2]-v[6]*v[3];
    flux[8] = 0;
    return flux;
}
double psiUpdate(double psi, double ch, double dt){
    double NewPsi = psi*std::exp(-1.0*dt*ch/0.18);
    return NewPsi;
}

std::tuple <double, double > ComputeTimeStep(const std::vector<std::vector<std::array <double , 9>>>& u  , double C, double dx, double dy, double gamma) {
    double maxSpeed = 0.0;
    for(int i=1; i<u.size() -2; ++i){
        for (int j=1; j<u[0].size()-2; ++j) {
            double rho = u[i][j][0];
            double mom_x = u[i][j][1];
            double mom_y = u[i][j][2];
            double mom_z = u[i][j][3];
            double E = u[i][j][4];
            double Bx = u[i][j][5];
            double By = u[i][j][6];
            double Bz = u[i][j][7];

            double u_x = mom_x / rho;
            double u_y = mom_y / rho;
            double u_z = mom_z / rho;
            double intermediate = 0.5 * rho *( u_x * u_x + u_y*u_y + u_z*u_z) + 0.5*(Bx*Bx + By*By + Bz*Bz);
            double BmagSquared = Bx*Bx + By*By + Bz*Bz;
            double pressure = (gamma - 1.0) * (E - intermediate);

            double sound_speed = std::sqrt(gamma * pressure / rho);
            double alfven_speed = std::abs(Bx) / std::sqrt(rho);
            double slow_ma_speed = std::sqrt( 0.5*(sound_speed*sound_speed + (BmagSquared/rho) - std::sqrt((sound_speed*sound_speed + BmagSquared/rho)*(sound_speed*sound_speed + BmagSquared/rho) - 4.0*sound_speed*sound_speed*Bx*Bx / rho)));
            
            // fast speed in x-direction
            double fast_ma_speed_x = sqrt(0.5 * (sound_speed*sound_speed + BmagSquared/rho + sqrt(pow(sound_speed*sound_speed + BmagSquared/rho, 2) - 4 * sound_speed*sound_speed * Bx*Bx / rho)));

            // fast speed in y-direction
            double fast_ma_speed_y = sqrt(0.5 * (sound_speed*sound_speed + BmagSquared/rho + sqrt(pow(sound_speed*sound_speed + BmagSquared/rho, 2) - 4 * sound_speed*sound_speed * By*By / rho)));
            
            // fast speed in z-direction
            double fast_ma_speed_z = sqrt(0.5 * (sound_speed*sound_speed + BmagSquared/rho + sqrt(pow(sound_speed*sound_speed + BmagSquared/rho, 2) - 4 * sound_speed*sound_speed * Bz*Bz / rho)));
            
            double speedx = std::abs(u_x) + fast_ma_speed_x;
            double speedy = std::abs(u_y) + fast_ma_speed_y;
            double speedz = std::abs(u_z) + fast_ma_speed_z;

            double speed = std::max(speedx,speedy);
            speed = std::max(speed , speedz);

            if (speed > maxSpeed) {
                maxSpeed = speed;
            }
        }
    }

    double dspace = std::min(dx, dy);
    double dt = C * dspace / maxSpeed;
     return {dt, maxSpeed};
}
//define a function to get the flux x is ui and y is ui+1 they are both arrays of length 3
//it will spit out an array of size 3 as well give each variable their flux
std::tuple<double,double> psiAndBx(const std::array<double ,9>&x,std::array<double ,9>y, double gamma, double ch){
    double psiL = x[8];
    double psiR = y[8];
    double BxL = x[5];
    double BxR = y[5];
    double BxTilde = 0.5*(BxL + BxR) - (1.0/(2.0*ch))*(psiR - psiL);
    double psiTilde = 0.5*(psiL + psiR) - (ch/2.0)*(BxR - BxL);
    return{BxTilde, psiTilde};
}
std::tuple<double,double> psiAndBy(std::array<double ,9>x,std::array<double ,9>y, double gamma, double ch){
    double psiL = x[8];
    double psiR = y[8];
    double ByL = x[6];
    double ByR = y[6];
    double ByTilde = 0.5*(ByL + ByR) - (1.0/(2.0*ch))*(psiR - psiL);
    double psiTilde = 0.5*(psiL + psiR) - (ch/2.0)*(ByR - ByL);
    return{ByTilde, psiTilde};
}
std::array<double, 9>  getFluxX(std::array<double, 9>  x , std::array<double, 9>  y, double dx , double dt,double gamma, double ch){
    //impliment general flux function 
    auto [ BxTilde, psiTilde] = psiAndBx(x,y,gamma, ch);
    x[5] = BxTilde;
    x[8] = psiTilde;
    y[5] = BxTilde;
    y[8] = psiTilde; 
    std::array<double, 9>  f_1 = FluxDefX(x, gamma); //f(ui)
    std::array<double, 9>  f_2 = FluxDefX(y, gamma); //f(i+1)
    std::array<double, 9>  uPlusHalf; //u_{i+1/2}

    for(int i=0; i<=8; ++i){
        uPlusHalf[i] = 0.5*(x[i] + y[i]) - 0.5*(dt/dx)*(f_2[i] - f_1[i]);
    }

    std::array<double, 9> RI_flux = FluxDefX(uPlusHalf, gamma); //richtmyer flux
    std::array<double, 9> LF_flux; //set up 3 length array for LF flux
    std::array<double, 9> FORCE_flux; //set up 3 length array for FORCE

    for (int i=0; i<=8; ++i) {
        LF_flux[i] = (0.5 * (dx/dt) * (x[i] - y[i])) + 0.5 * (f_1[i] + f_2[i]);
        FORCE_flux[i] = 0.5 * (LF_flux[i] + RI_flux[i]);
    }

    return FORCE_flux;
}
std::array<double, 9>  getFluxY(std::array<double, 9>  x , std::array<double, 9>  y, double dy , double dt,double gamma, double ch){
    //impliment general flux function 
    auto [ByTilde, psiTilde] = psiAndBy(x,y,gamma,ch);
    x[6] = ByTilde;
    x[8] = psiTilde;
    y[6] = ByTilde;
    y[8] = psiTilde;
    std::array<double, 9>  f_1 = FluxDefY(x, gamma); //f(ui)
    std::array<double, 9>  f_2 = FluxDefY(y, gamma); //f(i+1)
    std::array<double, 9>  uPlusHalf; //u_{i+1/2}

    for(int i=0; i<=8; ++i){
        uPlusHalf[i] = 0.5*(x[i] + y[i]) - 0.5*(dt/dy)*(f_2[i] - f_1[i]);
    }

    std::array<double, 9> RI_flux = FluxDefY(uPlusHalf, gamma); //richtmyer flux
    std::array<double, 9> LF_flux; //set up 3 length array for LF flux
    std::array<double, 9> FORCE_flux; //set up 3 length array for FORCE

    for (int i=0; i<=8; ++i) {
        LF_flux[i] = (0.5 * (dy/dt) * (x[i] - y[i])) + 0.5 * (f_1[i] + f_2[i]);
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
    r =0.0;
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

std::array<double, 9> getxiDeltas(const std::array<double, 9>& u_left, const std::array<double, 9>& u_center,  const std::array<double, 9>& u_right) {
    std::array<double, 9> slopes;
    
    for (int k = 0; k < 9; k++) {
        double backward_diff = u_center[k] - u_left[k];
        double forward_diff = u_right[k] - u_center[k];
        slopes[k] = minbee(backward_diff, forward_diff);
    }
    
    return slopes;
}

// this function takes in u and spits out uBarL and uBarR





int main() { 
    int nxCells = 256; 
    int nyCells = 256;
    double x0 = 0.0;
    double x1 = 800.0;
    double y0 = 0.0;
    double y1 = 800.0;
    double tStart = 0.0; //set the start and finish time steps the same
    double tStop = 80;
    double C = 0.75;
    double gamma = 5.0 / 3.0;
    double dx = (x1 - x0) / nxCells; 
    double dy = (y1 - y0) / nyCells;
    const double pi = 3.14159265358979323846;
    double omega =0.0;

    // Allocate matrices with 2 extra points for transmissive BCs
    std::vector<std::vector<std::array<double, 9> > > u;
    u.resize(nxCells+2, std::vector<std::array<double, 9> >(nyCells + 2)); //set up u

    std::vector<std::vector<std::array<double, 9> > > uPlus1;
    uPlus1.resize(nxCells+2, std::vector<std::array<double, 9> >(nyCells + 2)); //set up uPlus1

    std::vector<std::vector<std::array<double, 9> > > fluxX;
    fluxX.resize(nxCells+2, std::vector<std::array<double, 9> >(nyCells + 2)); //fluxX

    std::vector<std::vector<std::array<double, 9> > > fluxY;
    fluxY.resize(nxCells+2, std::vector<std::array<double, 9> >(nyCells + 2)); //fluxY

    std::vector<std::vector<std::array<double, 9> > > uBar;
    uBar.resize(nxCells+2, std::vector<std::array<double, 9> >(nyCells + 2)); //intermediate u

    std::vector<std::vector<std::array<double, 9> > > uBarLx;
    uBarLx.resize(nxCells+2, std::vector<std::array<double, 9> >(nyCells + 2)); //intermediate u

    std::vector<std::vector<std::array<double, 9> > > uBarLy;
    uBarLy.resize(nxCells+2, std::vector<std::array<double, 9> >(nyCells + 2)); //intermediate u

    std::vector<std::vector<std::array<double, 9> > > uBarRx;
    uBarRx.resize(nxCells+2, std::vector<std::array<double, 9> >(nyCells + 2)); //intermediate u

    std::vector<std::vector<std::array<double, 9> > > uBarRy;
    uBarRy.resize(nxCells+2, std::vector<std::array<double, 9> >(nyCells + 2)); //intermediate u

    std::vector<std::vector<std::array<double, 9> > > uBarHalfL;
    uBarHalfL.resize(nxCells+2, std::vector<std::array<double, 9> >(nyCells + 2)); //intermediate u

    std::vector<std::vector<std::array<double, 9> > > uBarHalfR;
    uBarHalfR.resize(nxCells+2, std::vector<std::array<double, 9> >(nyCells + 2)); //intermediate u

  
    // Initial conditions!

    for(int i = 0; i < u.size(); i++) {
        for(int j =0 ; j<u[0].size(); j++){
            // x 0 is at point i=1/2
            double x = x0 + (i-0.5) * dx;
            double y = y0 + (j-0.5) * dy;
            std::array<double, 9> prim;

            //Brio-Wu
            if(y <= 400 && x<= 400) {
                prim[0] = 1; // Density
                prim[1] = 0; // Velocity
                prim[2] = 0;
                prim[3] = 0;
                prim[4] = 1; // pressure
                prim[5] = 1; // magnetic field
                prim[6] = 0.75;
                prim[7] = 0; 
                prim[8] =0;
                } 
            else if(y <= 400 && x> 400) {
                prim[0] = 1; // Density
                prim[1] = 0; // Velocity
                prim[2] = 0;
                prim[3] = 0;
                prim[4] = 1; // pressure
                prim[5] = 1; // magnetic field
                prim[6] = 0.75;
                prim[7] = 0; 
                prim[8] =0;
            }
            else if(y > 400 && x<= 400) {
                prim[0] = 0.125; // Density
                prim[1] = 0; // Velocity
                prim[2] = 0;
                prim[3] = 0;
                prim[4] = 0.1; // pressure
                prim[5] = -1; // magnetic field
                prim[6] = 0.75;
                prim[7] = 0;  
                prim[0] =0;
                } 
            else {
                prim[0] = 0.125; // Density
                prim[1] = 0; // Velocity
                prim[2] = 0;
                prim[3] = 0;
                prim[4] = 0.1; // pressure
                prim[5] = -1; // magnetic field
                prim[6] = 0.75;
                prim[7] = 0; 
                prim[0] =0;
            }

            //Kelvin-Helmhotz
            // prim[0] = 1.0; // Density
            // prim[1] = 1.0/(2.0*std::tanh(20*y)); // Velocity
            // prim[2] = 0.1*std::sin(2.0*pi*x)*std::exp(-y*y/(0.1*0.1));
            // prim[3] = 0.0;
            // prim[4] = 1.0 / gamma; // pressure
            // prim[5] = 0.1*std::sqrt(1.0)*std::cos(pi/3.0); // magnetic field
            // prim[6] = 0.0;
            // prim[7] = 0.1*std::sqrt(1.0)*std::sin(pi/3.0); 
            // prim[8] = 0;

            //Orszag-Tang
            // prim[0] = gamma*gamma; // Density
            // prim[1] = -std::sin(2.0*pi*y); // Velocity
            // prim[2] = std::sin(2.0*pi*x);
            // prim[3] = 0.0;
            // prim[4] = gamma; // pressure
            // prim[5] = -std::sin(2.0*pi*y); // magnetic field
            // prim[6] = std::sin(4.0*pi*x);
            // prim[7] = 0.0; 
            // prim[8] = 0;


            u[i][j] = PrimitiveToConservative(prim, gamma);

        }
    }  

    

    auto [dt,ch] = ComputeTimeStep(u , C , dx, dy,  gamma); //the time steps

    double t = tStart;
    do {

        for(int j = 0 ; j<nyCells+2; j++){
            for(int i=0 ; i<nxCells+2; i++){
                u[i][j][8] = psiUpdate(u[i][j][8] , ch , 0.5*dt);
            }
        }
        auto [dt,ch] = ComputeTimeStep(u , C , dx, dy,gamma); 
        t = t + dt;
        std::cout << "t= " << t << " dt= " << dt << std::endl;

        // Transmissive boundary conditions
        for (int j = 0; j <= nyCells + 1; ++j) {
            u[0][j] = u[1][j];
            u[nxCells + 1][j] = u[nxCells][j];
        }
        for (int i = 0; i <= nxCells + 1; ++i) {
            u[i][0] = u[i][1];
            u[i][nyCells + 1] = u[i][nyCells];
        }

        // Compute uBar in x-direction
        for (int i = 1; i <= nxCells; ++i) {
            for (int j = 1; j <= nyCells; ++j) {
                for (int k = 0; k < 9; ++k) {
                    double DeltaPlus = u[i+1][j][k] - u[i][j][k];
                    double DeltaMinus = u[i][j][k] - u[i-1][j][k];
                    double r = DeltaMinus / (DeltaPlus + 1e-8);
                    r=0;
                    double xi = (r <= 0) ? 0 : ((r <= 1) ? r : std::fmin(1, 2.0 / (1 + r)));
                    double Delta = 0.5 * (1 + omega) * DeltaMinus + 0.5 * (1 - omega) * DeltaPlus;

                    uBarLx[i][j][k] = u[i][j][k] - 0.5 * xi * Delta;
                    uBarRx[i][j][k] = u[i][j][k] + 0.5 * xi * Delta;
                }
            }
        }

        

        // Predictor step (half time-step flux correction in both directions)
        for (int i = 1; i <= nxCells; ++i) {
            for (int j = 1; j <= nyCells; ++j) {
                for (int k = 0; k < 9; ++k) {
                    auto FxL = FluxDefX(uBarLx[i][j], gamma);
                    auto FxR = FluxDefX(uBarRx[i][j], gamma);

                    uBarHalfL[i][j][k] = uBarLx[i][j][k] - 0.5 * (dt / dx) * (FxR[k] - FxL[k]);

                    uBarHalfR[i][j][k] = uBarRx[i][j][k] - 0.5 * (dt / dx) * (FxR[k] - FxL[k]);
                }
            }
        }

        // Boundary conditions for uBarHalfL/R
        // Repeat similar to u

        // Compute fluxes in x-direction
        for (int i = 0; i <= nxCells; ++i) {
            for (int j = 1; j <= nyCells; ++j) {
                fluxX[i][j] = getFluxX(uBarHalfR[i][j], uBarHalfL[i+1][j], dx, dt, gamma, ch);
            }
        }

        for (int i = 1; i <= nxCells; ++i) {
            for (int j = 1; j <= nyCells; ++j) {
                for (int k = 0; k < 9; ++k) {
                    uBar[i][j][k] = u[i][j][k]- (dt / dx) * (fluxX[i][j][k] - fluxX[i-1][j][k]);
                }
            }
        }

        // Compute uBar in y-direction
        for (int i = 1; i <= nxCells; ++i) {
            for (int j = 1; j <= nyCells; ++j) {
                for (int k = 0; k < 9; ++k) {
                    double DeltaPlus = u[i][j+1][k] - u[i][j][k];
                    double DeltaMinus = u[i][j][k] - u[i][j-1][k];
                    double r = DeltaMinus / (DeltaPlus + 1e-8);
                    r=0;
                    double xi = (r <= 0) ? 0 : ((r <= 1) ? r : std::fmin(1, 2.0 / (1 + r)));
                    double Delta = 0.5 * (1 + omega) * DeltaMinus + 0.5 * (1 - omega) * DeltaPlus;

                    uBarLy[i][j][k] = uBar[i][j][k] - 0.5 * xi * Delta;
                    uBarRy[i][j][k] = uBar[i][j][k] + 0.5 * xi * Delta;
                }
            }
        }
        // Predictor step (half time-step flux correction in both directions)
        for (int i = 1; i <= nxCells; ++i) {
            for (int j = 1; j <= nyCells; ++j) {
                for (int k = 0; k < 9; ++k) {
                    auto FyL = FluxDefY(uBarLy[i][j], gamma);
                    auto FyR = FluxDefY(uBarRy[i][j], gamma);

                    uBarHalfL[i][j][k] = uBarLy[i][j][k]  - 0.5 * (dt / dy) * (FyR[k] - FyL[k]);

                    uBarHalfR[i][j][k] = uBarRy[i][j][k] - 0.5 * (dt / dy) * (FyR[k] - FyL[k]);
                }
            }
        }

        // Compute fluxes in y-direction
        for (int i = 1; i <= nxCells; ++i) {
            for (int j = 0; j <= nyCells; ++j) {
                fluxY[i][j] = getFluxY(uBarHalfR[i][j], uBarHalfL[i][j+1], dy, dt, gamma, ch);
            }
        }

        // Final update
        for (int i = 1; i <= nxCells; ++i) {
            for (int j = 1; j <= nyCells; ++j) {
                for (int k = 0; k < 9; ++k) {
                    uPlus1[i][j][k] = uBar[i][j][k] - (dt / dy) * (fluxY[i][j][k] - fluxY[i][j-1][k]);
                }
            }
        }
        for(int j = 0 ; j<nyCells+2; j++){
            for(int i=0 ; i<nxCells+2; i++){
                uPlus1[i][j][8] = psiUpdate(uPlus1[i][j][8] , ch , 0.5*dt);
            }
        }

        // Prepare for next step
        u = uPlus1;
    } while (t < tStop);

    //still need to convert it back to primitive
    //define final results


    std::vector<std::vector<std::array<double,9>> >results;
    results.resize(nxCells+2, std::vector<std::array<double, 9> >(nyCells + 2)); //set up results

    for(int j = 1; j < nyCells+2; j++) { 
        for(int i = 1; i < nxCells+2; i++) {
            results[i][j] = ConservativeToPrimitive(u[i][j], gamma);
        }
    }




    // Output the results
    std::ofstream output("MHD.dat");
    for(int j = 1; j < nyCells+2; j++) { 
        for(int i = 1; i < nxCells+2; i++) {
            double x = x0 + (i - 1)*dx;
            double y = y0 + (j - 1)*dy;
            output << x << " " << y << " " << results[i][j][0] << " " << results[i][j][1] << " " << results[i][j][2] << " " << results[i][j][3] <<  " " << results[i][j][4]<< " " << std::sqrt(results[i][j][5]*results[i][j][5] + results[i][j][6]*results[i][j][6]) / results[i][j][7] <<" " << results[i][j][6] <<  " " << results[i][j][7] << " " << results[i][j][8]<<std::endl;
        }
        output<<std::endl;
    }
}