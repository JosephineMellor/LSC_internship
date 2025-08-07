#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <array>
#include <filesystem>

std::array<double , 8> PrimitiveToConservative(const std::array<double , 8>& u , double gamma){
    std::array<double ,8> v;
    v[0] = u[0]; //density
    v[1] = u[0] * u[1]; //vx
    v[2] = u[0] * u[2]; //vy
    v[3] = u[0] * u[3]; //vz
    v[5] = u[5]; //Bx
    v[6] = u[6]; //By
    v[7] = u[7]; //Bz
    v[4] = u[4] / (gamma-1) + 0.5*u[0]*(u[1]*u[1] + u[2]*u[2] + u[3]*u[3]) + 0.5*(u[5]*u[5] + u[6]*u[6] + u[7]*u[7]); //energy
    return v;
}
std::array<double, 8> ConservativeToPrimitive(const std::array<double , 8>& u , double gamma){
    std::array<double , 8> v;
    v[0] = u[0];
    v[1] = u[1] / u[0];
    v[2] = u[2] / u[0];
    v[3] = u[3] / u[0];
    v[5] = u[5];
    v[6] = u[6];
    v[7] = u[7];
    v[4] = (u[4] - 0.5*v[0]*(v[1]*v[1] + v[2]*v[2] + v[3]*v[3]) - 0.5*(u[5]*u[5] + u[6]*u[6] + u[7]*u[7]))*(gamma -1);//pressure
    return v;
}
std::array<double , 8> FluxDef(const std::array<double , 8>& u , double gamma ){//conservative variables go into this function
    std::array<double , 8> v = ConservativeToPrimitive(u , gamma);
    std::array<double , 8> flux;
    flux[0] = u[1]; //rho v_x
    flux[1] = u[0]*v[1]*v[1] + v[4] + 0.5*( u[5]*u[5] + u[6]*u[6] + u[7]*u[7] ) - u[5]*u[5]; 
    flux[2] = u[0]*v[1]*v[2] - u[5]*u[6];
    flux[3] = u[0]*v[1]*v[3] - u[5]*u[7];
    flux[4] = (u[4] + v[4] + 0.5*( u[5]*u[5] + u[6]*u[6] + u[7]*u[7] ))*v[1] - (v[1]*u[5]+v[2]*u[6]+v[3]*u[7])*u[5];    flux[5] = 0;
    flux[6] = u[6]*v[1] - u[5]*v[2];
    flux[7] = u[7]*v[1] - u[5]*v[3];
    return flux;
}
double ComputeTimeStep(const std::vector<std::array<double,8>>& u , double C, double dx, double gamma) {
    double maxSpeed = 0.0;

    for (const auto& state : u) {
        double rho = state[0];
        double mom_x = state[1];
        double mom_y = state[2];
        double mom_z = state[3];
        double E = state[4];
        double Bx = state[5];
        double By = state[6];
        double Bz = state[7];

        double u_x = mom_x / rho;
        double u_y = mom_y / rho;
        double u_z = mom_z / rho;
        double intermediate = 0.5 * rho *( u_x * u_x + u_y*u_y + u_z*u_z) + 0.5*(Bx*Bx + By*By + Bz*Bz);
        double BmagSquared = Bx*Bx + By*By + Bz*Bz;
        double pressure = (gamma - 1.0) * (E - intermediate);

        double sound_speed = std::sqrt(gamma * pressure / rho);
        double alfven_speed = std::abs(Bx) / std::sqrt(rho);
        double slow_ma_speed = std::sqrt( 0.5*(sound_speed*sound_speed + (BmagSquared/rho) - std::sqrt((sound_speed*sound_speed + BmagSquared/rho)*(sound_speed*sound_speed + BmagSquared/rho) - 4.0*sound_speed*sound_speed*Bx*Bx / rho)));
        
        double fast_ma_speed = std::sqrt( 0.5*(sound_speed*sound_speed + (BmagSquared/rho) + std::sqrt((sound_speed*sound_speed + BmagSquared/rho)*(sound_speed*sound_speed + BmagSquared/rho) - 4.0*sound_speed*sound_speed*Bx*Bx / rho)));
        
        double speed = std::abs(u_x) + fast_ma_speed;

        if (speed > maxSpeed) {
            maxSpeed = speed;
        }
    }

    double dt = C * dx / maxSpeed;
     return dt;
}

std::array<double , 8> uHLL ( std::array<double , 8>x , std::array<double , 8> y , double gamma){
        //define all the variables
    double rhoL = x[0];
    double mom_xL = x[1];
    double mom_yL = x[2];
    double mom_zL = x[3];
    double EL = x[4];
    double BxL = x[5];
    double ByL = x[6];
    double BzL = x[7];

    double u_xL = mom_xL / rhoL;
    double u_yL = mom_yL / rhoL;
    double u_zL = mom_zL / rhoL;
    double intermediateL = 0.5 * rhoL *( u_xL * u_xL + u_yL*u_yL + u_zL*u_zL) + 0.5*(BxL*BxL + ByL*ByL + BzL*BzL);
    double BmagSquaredL = BxL*BxL + ByL*ByL + BzL*BzL;
    double pressureL = (gamma - 1.0) * (EL - intermediateL);
    double pTL = pressureL + 0.5*BmagSquaredL;

    double sound_speedL = std::sqrt(gamma * pressureL / rhoL);
    double alfven_speedL = std::abs(BxL) / std::sqrt(rhoL);
    double slow_ma_speedL = std::sqrt( 0.5*(sound_speedL*sound_speedL + (BmagSquaredL/rhoL) - std::sqrt((sound_speedL*sound_speedL + BmagSquaredL/rhoL)*(sound_speedL*sound_speedL + BmagSquaredL/rhoL) - 4.0*sound_speedL*sound_speedL*BxL*BxL / rhoL)));
    double fast_ma_speedL = std::sqrt( 0.5*(sound_speedL*sound_speedL + (BmagSquaredL/rhoL) + std::sqrt((sound_speedL*sound_speedL + BmagSquaredL/rhoL)*(sound_speedL*sound_speedL + BmagSquaredL/rhoL) - 4.0*sound_speedL*sound_speedL*BxL*BxL / rhoL)));
    
    double rhoR = y[0];
    double mom_xR = y[1];
    double mom_yR = y[2];
    double mom_zR = y[3];
    double ER = y[4];
    double BxR = y[5];
    double ByR = y[6];
    double BzR = y[7];

    double u_xR = mom_xR / rhoR;
    double u_yR = mom_yR / rhoR;
    double u_zR = mom_zR / rhoR;
    double intermediateR = 0.5 * rhoR * (u_xR * u_xR + u_yR * u_yR + u_zR * u_zR) + 0.5 * (BxR * BxR + ByR * ByR + BzR * BzR);
    double BmagSquaredR = BxR * BxR + ByR * ByR + BzR * BzR;
    double pressureR = (gamma - 1.0) * (ER - intermediateR);
    double pTR = pressureR + 0.5*BmagSquaredR;

    double sound_speedR = std::sqrt(gamma * pressureR / rhoR);
    double alfven_speedR = std::abs(BxR) / std::sqrt(rhoR);
    double slow_ma_speedR = std::sqrt(0.5 * (sound_speedR * sound_speedR + (BmagSquaredR / rhoR) - std::sqrt((sound_speedR * sound_speedR + BmagSquaredR / rhoR) * (sound_speedR * sound_speedR + BmagSquaredR / rhoR) - 4.0 * sound_speedR * sound_speedR * BxR * BxR / rhoR)));
    double fast_ma_speedR = std::sqrt(0.5 * (sound_speedR * sound_speedR + (BmagSquaredR / rhoR) + std::sqrt((sound_speedR * sound_speedR + BmagSquaredR / rhoR) * (sound_speedR * sound_speedR + BmagSquaredR / rhoR) - 4.0 * sound_speedR * sound_speedR * BxR * BxR / rhoR)));

    double sL = std::min(u_xL , u_xR) - std::max(fast_ma_speedL , fast_ma_speedR);
    double sR = std::max(u_xL , u_xR) + std::max(fast_ma_speedL , fast_ma_speedR);

    std::array<double , 8> uStar;
    auto FluxX = FluxDef(x, gamma);
    auto FluxY = FluxDef(y, gamma);
    for (int i = 0; i < 8; ++i) {
        uStar[i] = (sR * y[i] - sL * x[i] - FluxY[i] + FluxX[i]) / (sR - sL);
    }

    std::array<double , 8> uHLL;
    if (sL >= 0) {
        uHLL = x;
    } else if (sL < 0 && sR > 0) {
        uHLL = uStar;
    } else { // sR <= 0
        uHLL = y;
    }

    return uHLL;
}

std::array<double , 8> FluxHLL(std::array<double , 8>x , std::array<double , 8> y , double gamma){
    double rhoL = x[0];
    double mom_xL = x[1];
    double mom_yL = x[2];
    double mom_zL = x[3];
    double EL = x[4];
    double BxL = x[5];
    double ByL = x[6];
    double BzL = x[7];

    double u_xL = mom_xL / rhoL;
    double u_yL = mom_yL / rhoL;
    double u_zL = mom_zL / rhoL;
    double intermediateL = 0.5 * rhoL *( u_xL * u_xL + u_yL*u_yL + u_zL*u_zL) + 0.5*(BxL*BxL + ByL*ByL + BzL*BzL);
    double BmagSquaredL = BxL*BxL + ByL*ByL + BzL*BzL;
    double pressureL = (gamma - 1.0) * (EL - intermediateL);
    double pTL = pressureL + 0.5*BmagSquaredL;

    double sound_speedL = std::sqrt(gamma * pressureL / rhoL);
    double alfven_speedL = std::abs(BxL) / std::sqrt(rhoL);
    double slow_ma_speedL = std::sqrt( 0.5*(sound_speedL*sound_speedL + (BmagSquaredL/rhoL) - std::sqrt((sound_speedL*sound_speedL + BmagSquaredL/rhoL)*(sound_speedL*sound_speedL + BmagSquaredL/rhoL) - 4.0*sound_speedL*sound_speedL*BxL*BxL / rhoL)));
    double fast_ma_speedL = std::sqrt( 0.5*(sound_speedL*sound_speedL + (BmagSquaredL/rhoL) + std::sqrt((sound_speedL*sound_speedL + BmagSquaredL/rhoL)*(sound_speedL*sound_speedL + BmagSquaredL/rhoL) - 4.0*sound_speedL*sound_speedL*BxL*BxL / rhoL)));
    
    double rhoR = y[0];
    double mom_xR = y[1];
    double mom_yR = y[2];
    double mom_zR = y[3];
    double ER = y[4];
    double BxR = y[5];
    double ByR = y[6];
    double BzR = y[7];

    double u_xR = mom_xR / rhoR;
    double u_yR = mom_yR / rhoR;
    double u_zR = mom_zR / rhoR;
    double intermediateR = 0.5 * rhoR * (u_xR * u_xR + u_yR * u_yR + u_zR * u_zR) + 0.5 * (BxR * BxR + ByR * ByR + BzR * BzR);
    double BmagSquaredR = BxR * BxR + ByR * ByR + BzR * BzR;
    double pressureR = (gamma - 1.0) * (ER - intermediateR);
    double pTR = pressureR + 0.5*BmagSquaredR;

    double sound_speedR = std::sqrt(gamma * pressureR / rhoR);
    double alfven_speedR = std::abs(BxR) / std::sqrt(rhoR);
    double slow_ma_speedR = std::sqrt(0.5 * (sound_speedR * sound_speedR + (BmagSquaredR / rhoR) - std::sqrt((sound_speedR * sound_speedR + BmagSquaredR / rhoR) * (sound_speedR * sound_speedR + BmagSquaredR / rhoR) - 4.0 * sound_speedR * sound_speedR * BxR * BxR / rhoR)));
    double fast_ma_speedR = std::sqrt(0.5 * (sound_speedR * sound_speedR + (BmagSquaredR / rhoR) + std::sqrt((sound_speedR * sound_speedR + BmagSquaredR / rhoR) * (sound_speedR * sound_speedR + BmagSquaredR / rhoR) - 4.0 * sound_speedR * sound_speedR * BxR * BxR / rhoR)));

    double sL = std::min(u_xL , u_xR) - std::max(fast_ma_speedL , fast_ma_speedR);
    double sR = std::max(u_xL , u_xR) + std::max(fast_ma_speedL , fast_ma_speedR);

    std::array<double , 8> Flux;
    if(sL >= 0){
        Flux = FluxDef(x,gamma);
    }
    else if(sL < 0 && sR > 0){
        auto FluxX = FluxDef(x,gamma);
        auto FluxY = FluxDef(y,gamma);
        for(int i=0 ; i<8; ++i){
            Flux[i] = (sR*FluxX[i] - sL*FluxY[i] + sL*sR*(y[i] - x[i]))/(sR - sL);
        }
    }
    else { // sR <= 0
        Flux = FluxDef(y,gamma);
    }
    return Flux;
}

std::array<double , 8> FluxHLLC ( std::array<double , 8>x , std::array<double , 8> y , double gamma){
    //define all the variables
    double rhoL = x[0];
    double mom_xL = x[1];
    double mom_yL = x[2];
    double mom_zL = x[3];
    double EL = x[4];
    double BxL = x[5];
    double ByL = x[6];
    double BzL = x[7];

    double u_xL = mom_xL / rhoL;
    double u_yL = mom_yL / rhoL;
    double u_zL = mom_zL / rhoL;
    double intermediateL = 0.5 * rhoL *( u_xL * u_xL + u_yL*u_yL + u_zL*u_zL) + 0.5*(BxL*BxL + ByL*ByL + BzL*BzL);
    double BmagSquaredL = BxL*BxL + ByL*ByL + BzL*BzL;
    double pressureL = (gamma - 1.0) * (EL - intermediateL);
    double pTL = pressureL + 0.5*BmagSquaredL;

    double sound_speedL = std::sqrt(gamma * pressureL / rhoL);
    double alfven_speedL = std::abs(BxL) / std::sqrt(rhoL);
    double slow_ma_speedL = std::sqrt( 0.5*(sound_speedL*sound_speedL + (BmagSquaredL/rhoL) - std::sqrt((sound_speedL*sound_speedL + BmagSquaredL/rhoL)*(sound_speedL*sound_speedL + BmagSquaredL/rhoL) - 4.0*sound_speedL*sound_speedL*BxL*BxL / rhoL)));
    double fast_ma_speedL = std::sqrt( 0.5*(sound_speedL*sound_speedL + (BmagSquaredL/rhoL) + std::sqrt((sound_speedL*sound_speedL + BmagSquaredL/rhoL)*(sound_speedL*sound_speedL + BmagSquaredL/rhoL) - 4.0*sound_speedL*sound_speedL*BxL*BxL / rhoL)));
    
    double rhoR = y[0];
    double mom_xR = y[1];
    double mom_yR = y[2];
    double mom_zR = y[3];
    double ER = y[4];
    double BxR = y[5];
    double ByR = y[6];
    double BzR = y[7];

    double u_xR = mom_xR / rhoR;
    double u_yR = mom_yR / rhoR;
    double u_zR = mom_zR / rhoR;
    double intermediateR = 0.5 * rhoR * (u_xR * u_xR + u_yR * u_yR + u_zR * u_zR) + 0.5 * (BxR * BxR + ByR * ByR + BzR * BzR);
    double BmagSquaredR = BxR * BxR + ByR * ByR + BzR * BzR;
    double pressureR = (gamma - 1.0) * (ER - intermediateR);
    double pTR = pressureR + 0.5*BmagSquaredR;

    double sound_speedR = std::sqrt(gamma * pressureR / rhoR);
    double alfven_speedR = std::abs(BxR) / std::sqrt(rhoR);
    double slow_ma_speedR = std::sqrt(0.5 * (sound_speedR * sound_speedR + (BmagSquaredR / rhoR) - std::sqrt((sound_speedR * sound_speedR + BmagSquaredR / rhoR) * (sound_speedR * sound_speedR + BmagSquaredR / rhoR) - 4.0 * sound_speedR * sound_speedR * BxR * BxR / rhoR)));
    double fast_ma_speedR = std::sqrt(0.5 * (sound_speedR * sound_speedR + (BmagSquaredR / rhoR) + std::sqrt((sound_speedR * sound_speedR + BmagSquaredR / rhoR) * (sound_speedR * sound_speedR + BmagSquaredR / rhoR) - 4.0 * sound_speedR * sound_speedR * BxR * BxR / rhoR)));

    double sL = std::min(u_xL , u_xR) - std::max(fast_ma_speedL , fast_ma_speedR);
    double sR = std::max(u_xL , u_xR) + std::max(fast_ma_speedL , fast_ma_speedR);

    //double qStar = (rhoR*u_xR*(sR - u_xR) - rhoL*u_xL*(sL - u_xL) + pTL - pTR - BxL*BxL + BxR*BxR) / ( rhoR*(sR - u_xR) - rhoL*(sL - u_xL));
    double qStar = (pTR - pTL + rhoL * u_xL * (sL - u_xL) - rhoR * u_xR * (sR - u_xR)) / (rhoL * (sL - u_xL) - rhoR * (sR - u_xR));

    std::array<double , 8> u_HLL = uHLL(x,y,gamma);
    double BxHLL = u_HLL[5];
    double ByHLL = u_HLL[6];
    double BzHLL = u_HLL[7];

    double rhoHLL = u_HLL[0];
    double u_xHLL = u_HLL[1] / rhoHLL;  
    double u_yHLL = u_HLL[2] / rhoHLL;  
    double u_zHLL = u_HLL[3] / rhoHLL;  

    double BxStarL = BxL;
    double ByStarL = ByHLL;
    double BzStarL = BzHLL;

    double pStarL = rhoL*(sL - u_xL)*(qStar - u_xL) + pTL - BxL*BxL + BxStarL*BxStarL;
    //double pStarL = rhoL*(sL - u_xL)*(qStar - u_xL) + pressureL - 0.5*BxL*BxL + 0.5*BxStarL*BxStarL;

    double rhoStarL = rhoL*(sL - u_xL)/(sL - qStar);
    double mom_xStarL = rhoStarL * qStar;
    double mom_yStarL = mom_yL*(sL-u_xL)/(sL - qStar) - (BxStarL*ByStarL - BxL*ByL)/(sL - qStar);
    double mom_zStarL = mom_zL*(sL-u_xL)/(sL - qStar) - (BxStarL*BzStarL - BxL*BzL)/(sL - qStar);
    double u_xStarL = mom_xStarL / rhoStarL;
    double u_yStarL = mom_yStarL / rhoStarL;
    double u_zStarL = mom_zStarL / rhoStarL;

    double vDotB_StarL = BxStarL * u_xStarL + ByStarL * u_yStarL + BzStarL * u_zStarL;
    double vDotB_L = BxL * u_xL + ByL * u_yL + BzL * u_zL;
    double energyStarL = (EL * (sL - u_xL) - pTL * u_xL + pStarL * qStar + BxL * (vDotB_L - vDotB_StarL)) / (sL - qStar);


    std::array<double ,8> uStarL;
    uStarL[0] = rhoStarL;
    uStarL[1] = mom_xStarL;
    uStarL[2] = mom_yStarL;
    uStarL[3] = mom_zStarL;
    uStarL[4] = energyStarL;
    uStarL[5] = BxStarL;
    uStarL[6] = ByStarL;
    uStarL[7] = BzStarL;


    double BxStarR = BxR;
    double ByStarR = ByHLL;
    double BzStarR = BzHLL;

    //double pStarR = rhoR*(sR - u_xR)*(qStar - u_xR) + pressureR - BxR*BxR + BxStarR*BxStarR;
    double pStarR = rhoR*(sR - u_xR)*(qStar - u_xR) + pTR - 0.5*BxR*BxR + 0.5*BxStarR*BxStarR;

    double rhoStarR = rhoR*(sR - u_xR)/(sR - qStar);
    double mom_xStarR = rhoStarR * qStar;
    double mom_yStarR = mom_yR*(sR-u_xR)/(sR - qStar) - (BxStarR*ByStarR - BxR*ByR)/(sR - qStar);
    double mom_zStarR = mom_zR*(sR-u_xR)/(sR - qStar) - (BxStarR*BzStarR - BxR*BzR)/(sR - qStar);
    double u_xStarR = mom_xStarR / rhoStarR;
    double u_yStarR = mom_yStarR / rhoStarR;
    double u_zStarR = mom_zStarR / rhoStarR;

    double vDotB_StarR = BxStarR * u_xStarR + ByStarR * u_yStarR + BzStarR * u_zStarR;
    double vDotB_R = BxR * u_xR + ByR * u_yR + BzR * u_zR;
    double energyStarR = (ER * (sR - u_xR) - pTR * u_xR + pStarR * qStar + BxR * (vDotB_R - vDotB_StarR)) / (sR - qStar);


    std::array<double ,8> uStarR;
    uStarR[0] = rhoStarR;
    uStarR[1] = mom_xStarR;
    uStarR[2] = mom_yStarR;
    uStarR[3] = mom_zStarR;
    uStarR[4] = energyStarR;
    uStarR[5] = BxStarR;
    uStarR[6] = ByStarR;
    uStarR[7] = BzStarR;

    std::array<double , 8> FluxL = FluxDef(x,gamma);
    std::array<double , 8> FluxR = FluxDef(y,gamma);
    std::array<double , 8> FluxStarL = FluxDef(uStarL,gamma);
    std::array<double , 8> FluxStarR = FluxDef(uStarR , gamma);
    std::array<double , 8> FluxHLLC;

    if (sL >= 0) {
        FluxHLLC = FluxL;
    } else if (sL <= 0 && qStar >= 0) {
        for (int i = 0; i < 8; ++i)
            FluxHLLC[i] = FluxL[i] + sL * (uStarL[i] - x[i]);
    } else if (qStar <= 0 && sR >= 0) {
        for (int i = 0; i < 8; ++i)
            FluxHLLC[i] = FluxR[i] + sR * (uStarR[i] - y[i]);
    } else {
        FluxHLLC = FluxR;
    }
    return FluxHLLC;
}

int main() { 
    int nCells = 1000; //the distance between points is 0.01
    double x0 = 0.0;
    double x1 = 800;
    double tStart = 0.0; //set the start and finish time steps the same
    double tStop = 80;
    double C = 1.0;
    double gamma = 2.0;
    double omega = 0;
    double dx = (x1 - x0) / nCells; //the space steps 

    // Allocate matrices with 2 extra points for transmissive BCs
    std::vector<std::array<double,8>> u(nCells+2);
    std::vector<std::array<double,8>> uPlus1(nCells+2);
    std::vector<std::array<double,8>> flux(nCells+2);

    // Initial conditions!

    for(int i = 0; i < u.size(); i++) {
        // x 0 is at point i=1/2
        double x = x0 + (i-0.5) * dx;
        std::array<double, 8> prim;
        if(x <= 400) {
            prim[0] = 1; // Density
            prim[1] = 0; // Velocity
            prim[2] = 0;
            prim[3] = 0;
            prim[4] = 1; // pressure
            prim[5] = 0.75; // magnetic field
            prim[6] = 1;
            prim[7] = 0; 
            } else {
            prim[0] = 0.125; // Density
            prim[1] = 0; // Velocity
            prim[2] = 0;
            prim[3] = 0;
            prim[4] = 0.1; // pressure
            prim[5] = 0.75; // magnetic field
            prim[6] = -1;
            prim[7] = 0; 
        }

        u[i] = PrimitiveToConservative(prim, gamma);
    }

    double dt = ComputeTimeStep(u , C , dx, gamma); //the time steps


    double t = tStart;
    do {

        dt = ComputeTimeStep(u , C , dx, gamma); 
        t = t + dt;
        std::cout << "t= " << t << " dt= "<< dt << std::endl;

        // Trasmissive boundary conditions
        u[0] = u[1];
        u[nCells + 1] = u[nCells];


        for(int i = 0; i < nCells+1; i++) { //Define the fluxes
            // flux[i] corresponds to cell i+1/2 
            flux[i] = FluxHLLC( u[i], u[i+1] , gamma);
        }

        for(int i = 1; i < nCells+1; i++) { //Update the data
            for(int j=0; j<8; ++j){
                uPlus1[i][j] = u[i][j] - (dt/dx) * (flux[i][j] - flux[i-1][j]);
            }
        }

    
        // Now replace u with the updated data for the next time step

        u = uPlus1;
    } while (t < tStop);

    //convert it back to primitive and define final results

    std::vector<std::array<double,8>> results(nCells+2);

    for(int i=0; i<= results.size() -1; ++i){
        results[i] = ConservativeToPrimitive(u[i], gamma);
    }


    // Output the results
    std::ofstream output("MHD.dat");
    for (int i = 1; i <= nCells; ++i) {
        double x = x0 + (i - 1) * dx;
        output << x << " " << results[i][0] <<  " " << results[i][1] <<  " " << results[i][2] << " " << results[i][3] <<" " << results[i][4] <<" " << results[i][5] <<" " << results[i][6] <<" " << results[i][7] << std::endl;
        // std::cout << x << " " << u[i][0] <<  " " << u[i][1] <<  " " << u[i][2] << std::endl;
    }
}
