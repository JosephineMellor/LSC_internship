#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <array>
#include <filesystem>
#include <tuple>

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
            double fast_ma_speed_x = std::sqrt(0.5 * (sound_speed*sound_speed + BmagSquared/rho + std::sqrt(pow(sound_speed*sound_speed + BmagSquared/rho, 2) - 4 * sound_speed*sound_speed * Bx*Bx / rho)));

            // fast speed in y-direction
            double fast_ma_speed_y = std::sqrt(0.5 * (sound_speed*sound_speed + BmagSquared/rho + std::sqrt(pow(sound_speed*sound_speed + BmagSquared/rho, 2) - 4 * sound_speed*sound_speed * By*By / rho)));
            
            // fast speed in z-direction
            double fast_ma_speed_z = std::sqrt(0.5 * (sound_speed*sound_speed + BmagSquared/rho + std::sqrt(pow(sound_speed*sound_speed + BmagSquared/rho, 2) - 4 * sound_speed*sound_speed * Bz*Bz / rho)));
            
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

std::array<double , 9> uHLLx ( std::array<double , 9>x , std::array<double , 9> y , double gamma,double ch){
    auto [BxTilde , psiTilde] = psiAndBx(x,y,gamma,ch);
        //define all the variables
    double rhoL = x[0];
    double mom_xL = x[1];
    double mom_yL = x[2];
    double mom_zL = x[3];
    double EL = x[4];
    double BxL = BxTilde;
    double ByL = x[6];
    double BzL = x[7];

    double u_xL = mom_xL / rhoL;
    double u_yL = mom_yL / rhoL;
    double u_zL = mom_zL / rhoL;
    double intermediateL = 0.5 * rhoL *( u_xL * u_xL + u_yL*u_yL + u_zL*u_zL) + 0.5*(BxL*BxL + ByL*ByL + BzL*BzL);
    double BmagSquaredL = BxL*BxL + ByL*ByL + BzL*BzL;
    double pressureL = (gamma - 1.0) * (EL - intermediateL);

    double sound_speedL = std::sqrt(gamma * pressureL / rhoL);
    double fast_ma_speedL = std::sqrt( 0.5*(sound_speedL*sound_speedL + (BmagSquaredL/rhoL) + std::sqrt((sound_speedL*sound_speedL + BmagSquaredL/rhoL)*(sound_speedL*sound_speedL + BmagSquaredL/rhoL) - 4.0*sound_speedL*sound_speedL*BxL*BxL / rhoL)));
    
    double rhoR = y[0];
    double mom_xR = y[1];
    double mom_yR = y[2];
    double mom_zR = y[3];
    double ER = y[4];
    double BxR = BxTilde;
    double ByR = y[6];
    double BzR = y[7];

    double u_xR = mom_xR / rhoR;
    double u_yR = mom_yR / rhoR;
    double u_zR = mom_zR / rhoR;
    double intermediateR = 0.5 * rhoR * (u_xR * u_xR + u_yR * u_yR + u_zR * u_zR) + 0.5 * (BxR * BxR + ByR * ByR + BzR * BzR);
    double BmagSquaredR = BxR * BxR + ByR * ByR + BzR * BzR;
    double pressureR = (gamma - 1.0) * (ER - intermediateR);

    double sound_speedR = std::sqrt(gamma * pressureR / rhoR);
    double fast_ma_speedR = std::sqrt(0.5 * (sound_speedR * sound_speedR + (BmagSquaredR / rhoR) + std::sqrt((sound_speedR * sound_speedR + BmagSquaredR / rhoR) * (sound_speedR * sound_speedR + BmagSquaredR / rhoR) - 4.0 * sound_speedR * sound_speedR * BxR * BxR / rhoR)));

    double sL = std::min(u_xL , u_xR) - std::max(fast_ma_speedL , fast_ma_speedR);
    double sR = std::max(u_xL , u_xR) + std::max(fast_ma_speedL , fast_ma_speedR);

    std::array<double , 9> uStar;
    auto FluxX = FluxDefX(x, gamma);
    auto FluxY = FluxDefX(y, gamma);
    for (int i = 0; i < 9; ++i) {
        uStar[i] = (sR * y[i] - sL * x[i] - FluxY[i] + FluxX[i]) / (sR - sL);
    }

    std::array<double , 9> uHLL;
    if (sL >= 0) {
        uHLL = x;
    } else if (sL < 0 && sR > 0) {
        uHLL = uStar;
    } else { // sR <= 0
        uHLL = y;
    }

    return uHLL;
}

std::array<double , 9> uHLLy ( std::array<double , 9>x , std::array<double , 9> y , double gamma, double ch){
    auto [ByTilde , psiTilde] = psiAndBy(x,y,gamma,ch);
        //define all the variables
    double rhoL = x[0];
    double mom_xL = x[1];
    double mom_yL = x[2];
    double mom_zL = x[3];
    double EL = x[4];
    double BxL = x[5];
    double ByL = ByTilde;
    double BzL = x[7];

    double u_xL = mom_xL / rhoL;
    double u_yL = mom_yL / rhoL;
    double u_zL = mom_zL / rhoL;
    double intermediateL = 0.5 * rhoL *( u_xL * u_xL + u_yL*u_yL + u_zL*u_zL) + 0.5*(BxL*BxL + ByL*ByL + BzL*BzL);
    double BmagSquaredL = BxL*BxL + ByL*ByL + BzL*BzL;
    double pressureL = (gamma - 1.0) * (EL - intermediateL);
    double pTL = pressureL + 0.5*BmagSquaredL;

    double sound_speedL = std::sqrt(gamma * pressureL / rhoL);
    double fast_ma_speedL = std::sqrt( 0.5*(sound_speedL*sound_speedL + (BmagSquaredL/rhoL) + std::sqrt((sound_speedL*sound_speedL + BmagSquaredL/rhoL)*(sound_speedL*sound_speedL + BmagSquaredL/rhoL) - 4.0*sound_speedL*sound_speedL*ByL*ByL / rhoL)));
    
    double rhoR = y[0];
    double mom_xR = y[1];
    double mom_yR = y[2];
    double mom_zR = y[3];
    double ER = y[4];
    double BxR = y[5];
    double ByR = ByTilde;
    double BzR = y[7];

    double u_xR = mom_xR / rhoR;
    double u_yR = mom_yR / rhoR;
    double u_zR = mom_zR / rhoR;
    double intermediateR = 0.5 * rhoR * (u_xR * u_xR + u_yR * u_yR + u_zR * u_zR) + 0.5 * (BxR * BxR + ByR * ByR + BzR * BzR);
    double BmagSquaredR = BxR * BxR + ByR * ByR + BzR * BzR;
    double pressureR = (gamma - 1.0) * (ER - intermediateR);
    double pTR = pressureR + 0.5*BmagSquaredR;

    double sound_speedR = std::sqrt(gamma * pressureR / rhoR);
    double fast_ma_speedR = std::sqrt(0.5 * (sound_speedR * sound_speedR + (BmagSquaredR / rhoR) + std::sqrt((sound_speedR * sound_speedR + BmagSquaredR / rhoR) * (sound_speedR * sound_speedR + BmagSquaredR / rhoR) - 4.0 * sound_speedR * sound_speedR * ByR * ByR / rhoR)));

    double sL = std::min(u_yL , u_yR) - std::max(fast_ma_speedL , fast_ma_speedR);
    double sR = std::max(u_yL , u_yR) + std::max(fast_ma_speedL , fast_ma_speedR);

    std::array<double , 9> uStar;
    auto FluxX = FluxDefY(x, gamma);
    auto FluxY = FluxDefY(y, gamma);
    for (int i = 0; i < 9; ++i) {
        uStar[i] = (sR * y[i] - sL * x[i] - FluxY[i] + FluxX[i]) / (sR - sL);
    }

    std::array<double , 9> uHLL;
    if (sL >= 0) {
        uHLL = x;
    } else if (sL < 0 && sR > 0) {
        uHLL = uStar;
    } else { // sR <= 0
        uHLL = y;
    }

    return uHLL;
}

std::array<double , 9> FluxHLLCX ( std::array<double , 9>x , std::array<double , 9> y , double gamma , double ch){
    //define all the variables
    auto [BxTilde , psiTilde] = psiAndBx(x,y,gamma,ch);
    // x[5] = y[5] = BxTilde;
    // x[8] = y[8] = psiTilde;
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
    double BmagSquaredL = BxTilde*BxTilde + ByL*ByL + BzL*BzL;
    double pressureL = (gamma - 1.0) * (EL - intermediateL);
    double pTL = pressureL + 0.5*BmagSquaredL - BxTilde*BxTilde;

    double sound_speedL = std::sqrt(gamma * pressureL / rhoL);
    double fast_ma_speedL = std::sqrt( 0.5*(sound_speedL*sound_speedL + (BmagSquaredL/rhoL) + std::sqrt((sound_speedL*sound_speedL + BmagSquaredL/rhoL)*(sound_speedL*sound_speedL + BmagSquaredL/rhoL) - 4.0*sound_speedL*sound_speedL*BxTilde*BxTilde / rhoL)));
    
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
    double BmagSquaredR = BxTilde * BxTilde + ByR * ByR + BzR * BzR;
    double pressureR = (gamma - 1.0) * (ER - intermediateR);
    double pTR = pressureR + 0.5*BmagSquaredR - BxTilde*BxTilde;

    double sound_speedR = std::sqrt(gamma * pressureR / rhoR);
    double fast_ma_speedR = std::sqrt(0.5 * (sound_speedR * sound_speedR + (BmagSquaredR / rhoR) + std::sqrt((sound_speedR * sound_speedR + BmagSquaredR / rhoR) * (sound_speedR * sound_speedR + BmagSquaredR / rhoR) - 4.0 * sound_speedR * sound_speedR * BxTilde * BxTilde / rhoR)));

    double sL = std::min(u_xL , u_xR) - std::max(fast_ma_speedL , fast_ma_speedR);
    double sR = std::max(u_xL , u_xR) + std::max(fast_ma_speedL , fast_ma_speedR);

    //double qStar = (rhoR*u_xR*(sR - u_xR) - rhoL*u_xL*(sL - u_xL) + pTL - pTR - BxL*BxL + BxR*BxR) / ( rhoR*(sR - u_xR) - rhoL*(sL - u_xL));
    double qStar = (pTR - pTL + rhoL * u_xL * (sL - u_xL) - rhoR * u_xR * (sR - u_xR)) / (rhoL * (sL - u_xL) - rhoR * (sR - u_xR));


    std::array<double , 9> u_HLL = uHLLx(x,y,gamma, ch);
    double BxHLL = u_HLL[5];
    double ByHLL = u_HLL[6];
    double BzHLL = u_HLL[7];

    double rhoHLL = u_HLL[0];
    double u_xHLL = u_HLL[1] / rhoHLL;  
    double u_yHLL = u_HLL[2] / rhoHLL;  
    double u_zHLL = u_HLL[3] / rhoHLL;  

    double BxStarL = BxTilde;
    double ByStarL = ByHLL;
    double BzStarL = BzHLL;
    pressureL += 0.5*(BxTilde*BxTilde + ByL*ByL + BzL*BzL );

    //double pStarL = rhoL*(sL - u_xL)*(qStar - u_xL) + pTL - BxL*BxL + BxStarL*BxStarL;
    double pStarL = rhoL*(sL - u_xL)*(qStar - u_xL) + pressureL - 0.5*BxL*BxL + 0.5*BxStarL*BxStarL;

    double rhoStarL = rhoL*(sL - u_xL)/(sL - qStar);
    double mom_xStarL = rhoStarL * qStar;
    double mom_yStarL = mom_yL*(sL-u_xL)/(sL - qStar) - (BxStarL*ByStarL - BxL*ByL)/(sL - qStar);
    double mom_zStarL = mom_zL*(sL-u_xL)/(sL - qStar) - (BxStarL*BzStarL - BxL*BzL)/(sL - qStar);
    double u_xStarL = mom_xStarL / rhoStarL;
    double u_yStarL = mom_yStarL / rhoStarL;
    double u_zStarL = mom_zStarL / rhoStarL;

    double vDotB_StarL = BxStarL * u_xHLL + ByStarL * u_yHLL + BzStarL * u_zHLL;
    double vDotB_L = BxL * u_xL + ByL * u_yL + BzL * u_zL;
    double energyStarL = (EL * (sL - u_xL) - pressureL * u_xL + pStarL * qStar + BxL * vDotB_L - BxStarL* vDotB_StarL) / (sL - qStar);


    std::array<double ,9> uStarL;
    uStarL[0] = rhoStarL;
    uStarL[1] = mom_xStarL;
    uStarL[2] = mom_yStarL;
    uStarL[3] = mom_zStarL;
    uStarL[4] = energyStarL;
    uStarL[5] = BxStarL;
    uStarL[6] = ByStarL;
    uStarL[7] = BzStarL;
    uStarL[8] = psiTilde;


    double BxStarR = BxTilde;
    double ByStarR = ByHLL;
    double BzStarR = BzHLL;
    pressureR += 0.5*(BxTilde*BxTilde + ByR*ByR + BzR*BzR );

    //double pStarR = rhoR*(sR - u_xR)*(qStar - u_xR) + pressureR - BxR*BxR + BxStarR*BxStarR;
    double pStarR = rhoR*(sR - u_xR)*(qStar - u_xR) + pTR - 0.5*BxR*BxR + 0.5*BxStarR*BxStarR;

    double rhoStarR = rhoR*(sR - u_xR)/(sR - qStar);
    double mom_xStarR = rhoStarR * qStar;
    double mom_yStarR = mom_yR*(sR-u_xR)/(sR - qStar) - (BxStarR*ByStarR - BxR*ByR)/(sR - qStar);
    double mom_zStarR = mom_zR*(sR-u_xR)/(sR - qStar) - (BxStarR*BzStarR - BxR*BzR)/(sR - qStar);
    double u_xStarR = mom_xStarR / rhoStarR;
    double u_yStarR = mom_yStarR / rhoStarR;
    double u_zStarR = mom_zStarR / rhoStarR;

    double vDotB_StarR = BxStarR * u_xHLL + ByStarR * u_yHLL + BzStarR * u_zHLL;
    double vDotB_R = BxR * u_xR + ByR * u_yR + BzR * u_zR;
    double energyStarR = (ER * (sR - u_xR) - pressureR * u_xR + pStarR * qStar + BxR * vDotB_R - BxStarL*vDotB_StarR) / (sR - qStar);


    std::array<double ,9> uStarR;
    uStarR[0] = rhoStarR;
    uStarR[1] = mom_xStarR;
    uStarR[2] = mom_yStarR;
    uStarR[3] = mom_zStarR;
    uStarR[4] = energyStarR;
    uStarR[5] = BxStarR;
    uStarR[6] = ByStarR;
    uStarR[7] = BzStarR;
    uStarR[8] = psiTilde;

    std::array<double , 9> FluxL = FluxDefX(x,gamma);
    std::array<double , 9> FluxR = FluxDefX(y,gamma);
    std::array<double , 9> FluxStarL = FluxDefX(uStarL,gamma);
    std::array<double , 9> FluxStarR = FluxDefX(uStarR , gamma);
    std::array<double , 9> FluxHLLC;

    if (sL >= 0) {
        FluxHLLC = FluxL;
    } else if (sL <= 0 && qStar >= 0) {
        for (int i = 0; i < 9; ++i)
            FluxHLLC[i] = FluxL[i] + sL * (uStarL[i] - x[i]);
    } else if (qStar <= 0 && sR >= 0) {
        for (int i = 0; i < 9; ++i)
            FluxHLLC[i] = FluxR[i] + sR * (uStarR[i] - y[i]);
    } else {
        FluxHLLC = FluxR;
    }
    FluxHLLC[5] = psiTilde;
    FluxHLLC[8] = ch*ch*BxTilde;

    return FluxHLLC;
}

std::array<double , 9> FluxHLLCY ( std::array<double , 9>x , std::array<double , 9> y , double gamma, double ch){
    //define all the variables

    auto [ByTilde , psiTilde] = psiAndBy(x,y,gamma,ch);

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
    double BmagSquaredL = BxL*BxL + ByTilde*ByTilde + BzL*BzL;
    double pressureL = (gamma - 1.0) * (EL - intermediateL);
    double pTL = pressureL + 0.5*BmagSquaredL - ByTilde*ByTilde;

    double sound_speedL = std::sqrt(gamma * pressureL / rhoL);
    double fast_ma_speedL = std::sqrt( 0.5*(sound_speedL*sound_speedL + (BmagSquaredL/rhoL) + std::sqrt((sound_speedL*sound_speedL + BmagSquaredL/rhoL)*(sound_speedL*sound_speedL + BmagSquaredL/rhoL) - 4.0*sound_speedL*sound_speedL*ByTilde*ByTilde / rhoL)));
    
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
    double BmagSquaredR = BxR * BxR + ByTilde * ByTilde + BzR * BzR;
    double pressureR = (gamma - 1.0) * (ER - intermediateR);
    double pTR = pressureR + 0.5*BmagSquaredR- ByTilde*ByTilde;

    double sound_speedR = std::sqrt(gamma * pressureR / rhoR);
    double fast_ma_speedR = std::sqrt(0.5 * (sound_speedR * sound_speedR + (BmagSquaredR / rhoR) + std::sqrt((sound_speedR * sound_speedR + BmagSquaredR / rhoR) * (sound_speedR * sound_speedR + BmagSquaredR / rhoR) - 4.0 * sound_speedR * sound_speedR * ByTilde * ByTilde / rhoR)));

    double sL = std::min(u_yL , u_yR) - std::max(fast_ma_speedL , fast_ma_speedR);
    double sR = std::max(u_yL , u_yR) + std::max(fast_ma_speedL , fast_ma_speedR);

    double qStar = (pTR - pTL + rhoL * u_yL * (sL - u_yL) - rhoR * u_yR * (sR - u_yR) ) / (rhoL * (sL - u_yL) - rhoR * (sR - u_yR));



    std::array<double , 9> u_HLL = uHLLy(x,y,gamma, ch);
    double BxHLL = u_HLL[5];
    double ByHLL = u_HLL[6];
    double BzHLL = u_HLL[7];

    double rhoHLL = u_HLL[0];
    double u_xHLL = u_HLL[1] / rhoHLL;  
    double u_yHLL = u_HLL[2] / rhoHLL;  
    double u_zHLL = u_HLL[3] / rhoHLL;  

    double BxStarL = BxHLL;
    double ByStarL = ByTilde;
    double BzStarL = BzHLL;
    pressureL += 0.5*(BxL*BxL + ByTilde*ByTilde + BzL*BzL );

    //double pStarL = rhoL*(sL - u_yL)*(qStar - u_yL) + pTL - ByL*ByL + ByStarL*ByStarL;
    double pStarL = rhoL*(sL - u_yL)*(qStar - u_yL) + pressureL - 0.5*ByL*ByL + 0.5*ByStarL*ByStarL;

    double rhoStarL = rhoL*(sL - u_yL)/(sL - qStar);
    double mom_xStarL = mom_xL*(sL-u_yL)/(sL - qStar) - (BxStarL*ByStarL - BxL*ByL)/(sL - qStar);
    double mom_yStarL = rhoStarL*qStar;
    double mom_zStarL = mom_zL*(sL-u_yL)/(sL - qStar) - (ByStarL*BzStarL - ByL*BzL)/(sL - qStar);
    double u_xStarL = mom_xStarL / rhoStarL;
    double u_yStarL = mom_yStarL / rhoStarL;
    double u_zStarL = mom_zStarL / rhoStarL;

    double vDotB_StarL = BxStarL * u_xHLL + ByStarL * u_yHLL + BzStarL * u_zHLL;
    double vDotB_L = BxL * u_xL + ByL * u_yL + BzL * u_zL;
    double energyStarL = (EL * (sL - u_yL) - pressureL * u_yL + pStarL * qStar + ByL * vDotB_L - ByStarL* vDotB_StarL) / (sL - qStar);


    std::array<double ,9> uStarL;
    uStarL[0] = rhoStarL;
    uStarL[1] = mom_xStarL;
    uStarL[2] = mom_yStarL;
    uStarL[3] = mom_zStarL;
    uStarL[4] = energyStarL;
    uStarL[5] = BxStarL;
    uStarL[6] = ByStarL;
    uStarL[7] = BzStarL;
    uStarL[8] = psiTilde;


    double BxStarR = BxHLL;
    double ByStarR = ByTilde;
    double BzStarR = BzHLL;
    pressureR += 0.5*(BxR*BxR + ByTilde*ByTilde + BzR*BzR );

    double pStarR = rhoR*(sR - u_yR)*(qStar - u_yR) + pTR - 0.5*ByR*ByR + 0.5*ByStarR*ByStarR;
    //double pStarL = rhoL*(sL - u_xL)*(qStar - u_xL) + pressureL - 0.5*BxL*BxL + 0.5*BxStarL*BxStarL;

    double rhoStarR = rhoR*(sR - u_yR)/(sR - qStar);
    double mom_xStarR = mom_xR*(sR-u_yR)/(sR - qStar) - (BxStarR*ByStarR - BxR*ByR)/(sR - qStar);
    double mom_yStarR = rhoStarR*qStar;
    double mom_zStarR = mom_zR*(sR-u_yR)/(sR - qStar) - (ByStarR*BzStarR - ByR*BzR)/(sR - qStar);
    double u_xStarR = mom_xStarR / rhoStarR;
    double u_yStarR = mom_yStarR / rhoStarR;
    double u_zStarR = mom_zStarR / rhoStarR;

    double vDotB_StarR = BxStarR * u_xHLL + ByStarR * u_yHLL + BzStarR * u_zHLL;
    double vDotB_R = BxR * u_xR + ByR * u_yR + BzR * u_zR;
    double energyStarR = (ER * (sR - u_yR) - pressureR * u_yR + pStarR * qStar + ByR * vDotB_R - ByStarR*vDotB_StarR) / (sR - qStar);

    std::array<double ,9> uStarR;
    uStarR[0] = rhoStarR;
    uStarR[1] = mom_xStarR;
    uStarR[2] = mom_yStarR;
    uStarR[3] = mom_zStarR;
    uStarR[4] = energyStarR;
    uStarR[5] = BxStarR;
    uStarR[6] = ByStarR;
    uStarR[7] = BzStarR;
    uStarR[8] = psiTilde;

    std::array<double , 9> FluxL = FluxDefY(x,gamma);
    std::array<double , 9> FluxR = FluxDefY(y,gamma);
    std::array<double , 9> FluxStarL = FluxDefY(uStarL,gamma);
    std::array<double , 9> FluxStarR = FluxDefY(uStarR , gamma);
    std::array<double , 9> FluxHLLC;

    if (sL >= 0) {
        FluxHLLC = FluxL;
    } else if (sL <= 0 && qStar >= 0) {
        for (int i = 0; i < 9; ++i)
            FluxHLLC[i] = FluxL[i] + sL * (uStarL[i] - x[i]);
    } else if (qStar <= 0 && sR >= 0) {
        for (int i = 0; i < 9; ++i)
            FluxHLLC[i] = FluxR[i] + sR * (uStarR[i] - y[i]);
    } else {
        FluxHLLC = FluxR;
    }
    FluxHLLC[6] = psiTilde;
    FluxHLLC[8] = ch*ch*ByTilde;
    return FluxHLLC;
}

std::array<double , 9> FluxHLLX(std::array<double , 9>x , std::array<double , 9> y , double gamma,double ch){
    auto [BxTilde , psiTilde] = psiAndBx(x,y,gamma,ch);
    // x[5] = y[5] = BxTilde;
    // x[8] = y[8] = psiTilde;
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
    double pTL = pressureL + 0.5*BmagSquaredL - ByL*ByL;

    double sound_speedL = std::sqrt(gamma * pressureL / rhoL);
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
    double pTR = pressureR + 0.5*BmagSquaredR- ByR*ByR;

    double sound_speedR = std::sqrt(gamma * pressureR / rhoR);
    double fast_ma_speedR = std::sqrt(0.5 * (sound_speedR * sound_speedR + (BmagSquaredR / rhoR) + std::sqrt((sound_speedR * sound_speedR + BmagSquaredR / rhoR) * (sound_speedR * sound_speedR + BmagSquaredR / rhoR) - 4.0 * sound_speedR * sound_speedR * BxR * BxR / rhoR)));

    double sL = std::min(u_xL , u_xR) - std::max(fast_ma_speedL , fast_ma_speedR);
    double sR = std::max(u_xL , u_xR) + std::max(fast_ma_speedL , fast_ma_speedR);

    std::array<double , 9> Flux;
    if(sL >= 0){
        Flux = FluxDefX(x,gamma);
    }
    else if(sL < 0 && sR > 0){
        auto FluxX = FluxDefX(x,gamma);
        auto FluxY = FluxDefX(y,gamma);
        for(int i=0 ; i<9; ++i){
            Flux[i] = (sR*FluxX[i] - sL*FluxY[i] + sL*sR*(y[i] - x[i]))/(sR - sL);
        }
    }
    else { // sR <= 0
        Flux = FluxDefX(y,gamma);
    }
    Flux[5] = psiTilde;
    Flux[8] = BxTilde*ch*ch;
    return Flux;
}

std::array<double , 9> FluxHLLY(std::array<double , 9>x , std::array<double , 9> y , double gamma,double ch){
    auto [ByTilde , psiTilde] = psiAndBy(x,y,gamma,ch);
    // x[6] = y[6] = ByTilde;
    // x[8] = y[8] = psiTilde;
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
    double slow_ma_speedL = std::sqrt( 0.5*(sound_speedL*sound_speedL + (BmagSquaredL/rhoL) - std::sqrt((sound_speedL*sound_speedL + BmagSquaredL/rhoL)*(sound_speedL*sound_speedL + BmagSquaredL/rhoL) - 4.0*sound_speedL*sound_speedL*ByL*ByL / rhoL)));
    double fast_ma_speedL = std::sqrt( 0.5*(sound_speedL*sound_speedL + (BmagSquaredL/rhoL) + std::sqrt((sound_speedL*sound_speedL + BmagSquaredL/rhoL)*(sound_speedL*sound_speedL + BmagSquaredL/rhoL) - 4.0*sound_speedL*sound_speedL*ByL*ByL / rhoL)));
    
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
    double slow_ma_speedR = std::sqrt(0.5 * (sound_speedR * sound_speedR + (BmagSquaredR / rhoR) - std::sqrt((sound_speedR * sound_speedR + BmagSquaredR / rhoR) * (sound_speedR * sound_speedR + BmagSquaredR / rhoR) - 4.0 * sound_speedR * sound_speedR * ByR * ByR / rhoR)));
    double fast_ma_speedR = std::sqrt(0.5 * (sound_speedR * sound_speedR + (BmagSquaredR / rhoR) + std::sqrt((sound_speedR * sound_speedR + BmagSquaredR / rhoR) * (sound_speedR * sound_speedR + BmagSquaredR / rhoR) - 4.0 * sound_speedR * sound_speedR * ByR * ByR / rhoR)));

    double sL = std::min(u_yL , u_yR) - std::max(fast_ma_speedL , fast_ma_speedR);
    double sR = std::max(u_yL , u_yR) + std::max(fast_ma_speedL , fast_ma_speedR);

    std::array<double , 9> Flux;
    if(sL >= 0){
        Flux = FluxDefY(x,gamma);
    }
    else if(sL < 0 && sR > 0){
        auto FluxX = FluxDefY(x,gamma);
        auto FluxY = FluxDefY(y,gamma);
        for(int i=0 ; i<9; ++i){
            Flux[i] = (sR*FluxX[i] - sL*FluxY[i] + sL*sR*(y[i] - x[i]))/(sR - sL);
        }
    }
    else { // sR <= 0
        Flux = FluxDefY(y,gamma);
    }
    Flux[6] = psiTilde;
    Flux[8] = ByTilde*ch*ch;
    return Flux;
}

void applyBoundaryConditions(std::vector<std::vector<std::array<double, 9>>>& u, int nxCells, int nyCells) {
    // Left and Right boundaries periodic
    for (int j = 0; j < nyCells + 2; ++j) {
        u[0][j] = u[nxCells][j]; 
        u[1][j] = u[nxCells+1][j];      // Left boundary
        u[nxCells + 2][j] = u[2][j]; 
        u[nxCells + 3][j] = u[3][j];  // Right boundary
    }

    // Left and Right boundaries transmissive
    // for (int j = 0; j < nyCells + 2; ++j) {
    //     u[0][j] = u[2][j]; 
    //     u[1][j] = u[3][j];      // Left boundary
    //     u[nxCells + 2][j] = u[nxCells][j]; 
    //     u[nxCells + 3][j] = u[nxCells + 1][j];  // Right boundary
    // }

    // Bottom and Top boundaries, reflective
    // for (int i = 0; i < nxCells + 2; ++i) {//reflect in u_y and B_y
    //     u[i][0] = u[i][1];        // Bottom boundary
    //     u[i][nyCells + 1] = u[i][nyCells];  // Top boundary
    //     u[i][0][2] = -u[i][1][2]; 
    //     u[i][0][6] = -u[i][1][6]; 
    //     u[i][nyCells + 1][2] = -u[i][nyCells][2]; 
    //     u[i][nyCells + 1][6] = -u[i][nyCells][6]; 
    // }

    //bottom and top periodic
    for (int i = 0; i < nxCells + 2; ++i) {
        u[i][0] = u[i][nyCells];      // Bottom boundary
        u[i][1] = u[i][nyCells+1];    
        u[i][nyCells + 2] = u[i][2];  // Top boundary
        u[i][nyCells + 3] = u[i][3];  
    }

    //bottom and top transmissive
    // for (int i = 0; i < nxCells + 2; ++i) {
    //     u[i][0] = u[i][2];      // Bottom boundary
    //     u[i][1] = u[i][3];    
    //     u[i][nyCells + 2] = u[i][nyCells];  // Top boundary
    //     u[i][nyCells + 3] = u[i][nyCells + 1];  
    // }
}



int main() { 
    int nxCells = 800; 
    int nyCells = 800;
    double x0 = 0.0;
    double x1 = 1.0;
    double y0 = 0.0;
    double y1 = 1.0;
    double tStart = 0.0; //set the start and finish time steps the same
    double tStop = 0.5;
    double C = 0.8;
    double gamma = 5.0/3.0;
    double dx = (x1 - x0) / nxCells; 
    double dy = (y1 - y0) / nyCells;
    const double pi = 3.14159265358979323846;

    std::vector<std::vector<std::array<double, 9> > > u;
    u.resize(nxCells+4, std::vector<std::array<double, 9> >(nyCells + 4)); //set up u

    std::vector<std::vector<std::array<double, 9> > > v;
    v.resize(nxCells+4, std::vector<std::array<double, 9> >(nyCells + 4)); //set up v

    std::vector<std::vector<std::array<double, 9> > > uPlus1;
    uPlus1.resize(nxCells+4, std::vector<std::array<double, 9> >(nyCells + 4)); //set up uPlus1

    std::vector<std::vector<std::array<double, 9> > > fluxX;
    fluxX.resize(nxCells+4, std::vector<std::array<double, 9> >(nyCells + 4)); //fluxX

    std::vector<std::vector<std::array<double, 9> > > fluxY;
    fluxY.resize(nxCells+4, std::vector<std::array<double, 9> >(nyCells + 4)); //fluxY

    std::vector<std::vector<std::array<double, 9> > > uBar;
    uBar.resize(nxCells+4, std::vector<std::array<double, 9> >(nyCells + 4)); //intermediate u

    std::vector<std::vector<std::array<double, 9> > > uBarR_x;
    uBarR_x.resize(nxCells+4, std::vector<std::array<double, 9> >(nyCells + 4)); //intermediate u

    std::vector<std::vector<std::array<double, 9> > > uBarL_x;
    uBarL_x.resize(nxCells+4, std::vector<std::array<double, 9> >(nyCells + 4)); //intermediate u
    
    std::vector<std::vector<std::array<double, 9> > > uBarB_y;
    uBarB_y.resize(nxCells+4, std::vector<std::array<double, 9> >(nyCells + 4)); //intermediate u

    std::vector<std::vector<std::array<double, 9> > > uBarT_y;
    uBarT_y.resize(nxCells+4, std::vector<std::array<double, 9> >(nyCells + 4)); //intermediate u

    std::vector<std::vector<std::array<double, 9> > > uBarHalf_yB;
    uBarHalf_yB.resize(nxCells+4, std::vector<std::array<double, 9> >(nyCells + 4)); //intermediate u

    std::vector<std::vector<std::array<double, 9> > > uBarHalf_yT;
    uBarHalf_yT.resize(nxCells+4, std::vector<std::array<double, 9> >(nyCells + 4)); //intermediate u

    std::vector<std::vector<std::array<double, 9> > > uBarHalf_xL;
    uBarHalf_xL.resize(nxCells+4, std::vector<std::array<double, 9> >(nyCells + 4)); //intermediate u

    std::vector<std::vector<std::array<double, 9> > > uBarHalf_xR;
    uBarHalf_xR.resize(nxCells+4, std::vector<std::array<double, 9> >(nyCells + 4)); //intermediate u

    

    // Initial conditions!

    for(int i = 0; i < u.size(); i++) {
        for(int j =0 ; j<u[0].size(); j++){
            // x 0 is at point i=1/2
            double x = x0 + (i-0.5) * dx;
            double y = y0 + (j-0.5) * dy;
            std::array<double, 9> prim;

            //Brio-Wu 
            // if(y <= 400 && x<= 400) {
            //     prim[0] = 1; // Density
            //     prim[1] = 0; // Velocity
            //     prim[2] = 0;
            //     prim[3] = 0;
            //     prim[4] = 1; // pressure
            //     prim[5] = 1.0;// magnetic field
            //     prim[6] = 0.75;
            //     prim[7] = 0; 
            // }
            // else if(y <= 400 && x> 400) {
            //     prim[0] = 1; // Density
            //     prim[1] = 0; // Velocity
            //     prim[2] = 0;
            //     prim[3] = 0;
            //     prim[4] = 1; // pressure
            //     prim[5] = 1.0; // magnetic field
            //     prim[6] = 0.75;
            //     prim[7] = 0; 
            // }
            // else if(y > 400 && x<= 400) {
            //     prim[0] = 0.125; // Density
            //     prim[1] = 0; // Velocity
            //     prim[2] = 0;
            //     prim[3] = 0;
            //     prim[4] = 0.1; // pressure
            //     prim[5] = -1.0; // magnetic field
            //     prim[6] = 0.75;
            //     prim[7] = 0;  
            //     }
            // else{
            //     prim[0] = 0.125; // Density
            //     prim[1] = 0; // Velocity
            //     prim[2] = 0;
            //     prim[3] = 0;
            //     prim[4] = 0.1; // pressure
            //     prim[5] = -1.0; // magnetic field
            //     prim[6] = 0.75;
            //     prim[7] = 0;  
            //     } 
            
                
            // double cos45 = 0.70710678118;
            // double sin45 = 0.70710678118;

            // if (y <= x) {
            //     prim[0] = 1;      // Density
            //     prim[1] = 0;      // Velocity x
            //     prim[2] = 0;      // Velocity y
            //     prim[3] = 0;      // Velocity z
            //     prim[4] = 1;      // Pressure
                
            //     double Bx_old = 1.0;
            //     double By_old = 0.75;
            //     prim[5] = Bx_old * cos45 - By_old * sin45;  // rotated Bx
            //     prim[6] = Bx_old * sin45 + By_old * cos45;  // rotated By
            //     prim[7] = 0;      // Bz
            // } else {
            //     prim[0] = 0.125;
            //     prim[1] = 0;
            //     prim[2] = 0;
            //     prim[3] = 0;
            //     prim[4] = 0.1;
                
            //     double Bx_old = -1.0;
            //     double By_old = 0.75;
            //     prim[5] = Bx_old * cos45 - By_old * sin45;
            //     prim[6] = Bx_old * sin45 + By_old * cos45;
            //     prim[7] = 0;
            // }


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
            prim[0] = gamma*gamma; // Density
            prim[1] = -std::sin(2.0*pi*y); // Velocity
            prim[2] = std::sin(2.0*pi*x);
            prim[3] = 0.0;
            prim[4] = gamma; // pressure
            prim[5] = -std::sin(2.0*pi*y); // magnetic field
            prim[6] = std::sin(4.0*pi*x);
            prim[7] = 0.0; 
            prim[8] = 0;


            u[i][j] = PrimitiveToConservative(prim, gamma);

        }
    }

    auto [dt, maxSpeed] = ComputeTimeStep(u , C , dx, dy, gamma); //the time steps

    double t = tStart;
    do {

        for(int j = 0 ; j<nyCells+4; j++){
            for(int i=0 ; i<nxCells+4; i++){
                u[i][j][8] = psiUpdate(u[i][j][8] , maxSpeed , 0.5*dt);
            }
        }

        auto [dt,maxSpeed] = ComputeTimeStep(u, C, dx, dy, gamma); 
        t += dt;
        std::cout << "t= " << t << " dt= " << dt << " ch= "<< maxSpeed<<std::endl;

        applyBoundaryConditions( u , nxCells , nyCells);

        // Reconstruct left/right states in x-direction
        for (int i = 1; i < nxCells+3; ++i) {
            for (int j = 0; j < nyCells+4; ++j) {
                v[i][j] = ConservativeToPrimitive(u[i][j],gamma);
                for (int k = 0; k < 9; ++k) {
                    double dPlus = v[i+1][j][k] - v[i][j][k];
                    double dMinus = v[i][j][k] - v[i-1][j][k];//calculate limiting with primitive
                    double r = dMinus / (dPlus + 1e-8);
                    dPlus = u[i+1][j][k] - u[i][j][k];
                    dMinus = u[i][j][k] - u[i-1][j][k];
                    r=0;
                    double xi ;
                    if(r<=0){ xi=0;}
                    else{ 
                        xi=std::fmin(2.0*r/(1+r), 2.0/(1+r));  
                    }
                    
                    double Delta = 0.5*dMinus + 0.5*dPlus;
                    uBarL_x[i][j][k] = u[i][j][k] - 0.5 * xi * Delta;
                    uBarR_x[i][j][k] = u[i][j][k] + 0.5 * xi * Delta;
                }
            }
        }

        // Predictor Step: update uBar to half timestep
        for (int i = 1; i < nxCells+3; ++i) {
            for (int j = 1; j < nyCells+3; ++j) {
                std::array<double,9> fluxLx = FluxDefX(uBarL_x[i][j], gamma);
                std::array<double,9> fluxRx = FluxDefX(uBarR_x[i][j], gamma);
                for (int k = 0; k < 9; ++k) {
                    uBarHalf_xL[i][j][k] = uBarL_x[i][j][k] - 0.5 * (dt / dx) * (fluxRx[k] - fluxLx[k]);
                    uBarHalf_xR[i][j][k] = uBarR_x[i][j][k] - 0.5 * (dt / dx) * (fluxRx[k] - fluxLx[k]);
                }
            }
        }

        // Apply boundary conditions to predicted states (optional depending on implementation)

        applyBoundaryConditions(uBarHalf_xL , nxCells , nyCells);
        applyBoundaryConditions(uBarHalf_xR , nxCells , nyCells);

        // Compute fluxes in x directions
        for (int i = 1; i < nxCells+3; ++i) {
            for (int j = 1; j < nyCells+3; ++j) {
                fluxX[i][j] = FluxHLLCX(uBarHalf_xR[i][j], uBarHalf_xL[i+1][j], gamma, maxSpeed);
            }
        }

        for (int i = 2; i < nxCells+2; ++i) {
            for (int j = 2; j < nyCells+2; ++j) {
                for (int k = 0; k < 9; ++k) {
                    uBar[i][j][k] = u[i][j][k] - (dt / dx) * (fluxX[i][j][k] - fluxX[i-1][j][k]);
                }
            }
        }

        applyBoundaryConditions(uBar , nxCells , nyCells);

        // Reconstruct bottom/top states in y-direction
        for (int i = 0; i < nxCells+4; ++i) {
            for (int j = 1; j < nyCells+3; ++j) {
                v[i][j] = ConservativeToPrimitive(u[i][j],gamma);
                for (int k = 0; k < 9; ++k) {
                    double dPlus = v[i][j+1][k] - v[i][j][k];
                    double dMinus = v[i][j][k] - v[i][j-1][k];//calculate r with primitive
                    double r = dMinus / (dPlus + 1e-8);
                    dPlus = uBar[i][j+1][k] - uBar[i][j][k];
                    dMinus = uBar[i][j][k] - uBar[i][j-1][k];
                    r=0;
                    double xi ;
                    if(r<=0){ xi=0;}
                    else{ 
                        xi=std::fmin(2.0*r/(1+r), 2.0/(1+r));  
                    }

                    double Delta = 0.5*dMinus + 0.5*dPlus;
                    uBarB_y[i][j][k] = uBar[i][j][k] - 0.5 * xi * Delta;
                    uBarT_y[i][j][k] = uBar[i][j][k] + 0.5 * xi * Delta;
                }
            }
        }

        // Predictor Step: update uBar to half timestep
         for (int i = 1; i < nxCells+3; ++i) {
            for (int j = 1; j < nyCells+3; ++j) {
                std::array<double,9> fluxBy = FluxDefY(uBarB_y[i][j], gamma);
                std::array<double,9> fluxTy = FluxDefY(uBarT_y[i][j], gamma);
                for (int k = 0; k < 9; ++k) {
                    uBarHalf_yB[i][j][k] = uBarB_y[i][j][k] - 0.5 * (dt / dy) * (fluxTy[k] - fluxBy[k]);
                    uBarHalf_yT[i][j][k] = uBarT_y[i][j][k] - 0.5 * (dt / dy) * (fluxTy[k] - fluxBy[k]);
                }
            }
        }

        applyBoundaryConditions(uBarHalf_yT , nxCells , nyCells);
        applyBoundaryConditions(uBarHalf_yB , nxCells , nyCells);

        for (int i = 1; i < nxCells+3; ++i) {
            for (int j = 1; j < nyCells+3; ++j) {
                fluxY[i][j] = FluxHLLCY(uBarHalf_yT[i][j], uBarHalf_yB[i][j+1], gamma, maxSpeed);
            }
        }

        // Final update
        for (int i = 2; i < nxCells+2; ++i) {
            for (int j = 2; j < nyCells+2; ++j) {
                for (int k = 0; k < 9; ++k) {
                    uPlus1[i][j][k] = uBar[i][j][k] - (dt / dy) * (fluxY[i][j][k] - fluxY[i][j-1][k]);
                }
            }
        }

        applyBoundaryConditions(uPlus1 , nxCells , nyCells);

        //update for source
        for(int j = 0 ; j<nyCells+4; j++){
            for(int i=0 ; i<nxCells+4; i++){
                uPlus1[i][j][8] = psiUpdate(uPlus1[i][j][8] , maxSpeed , 0.5*dt);
            }
        }


        u=uPlus1;

    } while (t < tStop);


    applyBoundaryConditions(u , nxCells , nyCells);

    //convert it back to primitive and define final results

    std::vector<std::vector<std::array<double,9>> >results;
    results.resize(nxCells+4, std::vector<std::array<double, 9> >(nyCells + 4)); //set up results

    for(int j = 0; j < nyCells+4; j++) { 
        for(int i = 0; i < nxCells+4; i++) {
            results[i][j] = ConservativeToPrimitive(u[i][j], gamma);
        }
    }




    // Output the results
    std::ofstream output("MHD2.dat");
    for(int j = 1; j < nyCells+4; j++) { 
        for(int i = 1; i < nxCells+4; i++) {
            double x = x0 + (i - 1)*dx;
            double y = y0 + (j - 1)*dy;
            output << x << " " << y << " " << results[i][j][0] << " " << results[i][j][1] << " " << results[i][j][2] << " " << results[i][j][3] <<  " " << results[i][j][4]<< " " << std::sqrt(results[i][j][5]*results[i][j][5] + results[i][j][6]*results[i][j][6]) / results[i][j][7] <<" " << results[i][j][6] <<  " " << results[i][j][7] << " " << results[i][j][8]<<std::endl;
        }
        output<<std::endl;
    }
}
