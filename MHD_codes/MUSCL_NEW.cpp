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

std::array<double , 9> FluxDefX(const std::array<double , 9>& u , double gamma , double ch){//conservative variables go into this function
    std::array<double , 9> v = ConservativeToPrimitive(u , gamma);
    std::array<double , 9> flux;
    flux[0] = u[1]; //rho v_x
    flux[1] = u[0]*v[1]*v[1] + v[4] + 0.5*( u[5]*u[5] + u[6]*u[6] + u[7]*u[7] ) - u[5]*u[5]; 
    flux[2] = u[0]*v[1]*v[2] - u[5]*u[6];
    flux[3] = u[0]*v[1]*v[3] - u[5]*u[7];
    flux[4] = (u[4] + v[4] + 0.5*( u[5]*u[5] + u[6]*u[6] + u[7]*u[7] ))*v[1] - (v[1]*v[5]+v[2]*v[6]+v[3]*v[7])*v[5];
    flux[5] = u[8];
    flux[6] = v[6]*v[1]-v[5]*v[2];
    flux[7] = v[7]*v[1]-v[5]*v[3];
    flux[8] = ch*ch*u[5];
    return flux;
}

std::array<double , 9> FluxDefY(const std::array<double , 9>& u , double gamma, double ch){//conservative variables go into this function
    std::array<double , 9> v = ConservativeToPrimitive(u , gamma);
    std::array<double , 9> flux;
    flux[0] = u[2]; //rho v_y
    flux[1] = u[0]*v[1]*v[2] - v[5]*v[6];
    flux[2] = u[0]*v[2]*v[2]+ v[4] + 0.5*( u[5]*u[5] + u[6]*u[6] + u[7]*u[7] ) - u[6]*u[6];
    flux[3] = u[0]*v[2]*v[3] - u[6]*u[7];
    flux[4] = (u[4] + v[4] + 0.5*( u[5]*u[5] + u[6]*u[6] + u[7]*u[7] ))*v[2] - (v[1]*v[5]+v[2]*v[6]+v[3]*v[7])*v[6];
    flux[5] = v[5]*v[2] - v[6]*v[1];
    flux[6] = u[8];
    flux[7] = v[7]*v[2]-v[6]*v[3];
    flux[8] = ch*ch*u[6];
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
std::tuple<double,double> psiAndBx(const std::array<double ,9>&x,const std::array<double ,9>y, double gamma, double ch){
    double psiL = x[8];
    double psiR = y[8];
    double BxL = x[5];
    double BxR = y[5];
    double BxTilde = 0.5*(BxL + BxR) - (1.0/(2.0*ch))*(psiR - psiL);
    double psiTilde = 0.5*(psiL + psiR) - (ch/2.0)*(BxR - BxL);
    return{BxTilde, psiTilde};
}
std::tuple<double,double> psiAndBy(const std::array<double ,9>&x,const std::array<double ,9>&y, double gamma, double ch){
    double psiL = x[8];
    double psiR = y[8];
    double ByL = x[6];
    double ByR = y[6];
    double ByTilde = 0.5*(ByL + ByR) - (1.0/(2.0*ch))*(psiR - psiL);
    double psiTilde = 0.5*(psiL + psiR) - (ch/2.0)*(ByR - ByL);
    return{ByTilde, psiTilde};
}

std::tuple<double,double> speedsx( const std::array<double , 9>&x , double gamma, double BxTilde){
    
    // Collect conservative and primitive variables
    auto primX = ConservativeToPrimitive(x,gamma);
    double rho = x[0];
    double u_x = primX[1];
    double u_y = primX[2];
    double u_z = primX[3];
    double E = x[4];
    double p = primX[4];
    double Bx = x[5];
    double By = x[6];
    double Bz = x[7];
    
    // Calculate magnetic magnitude and total pressure
    double B2 = Bx*Bx + By*By + Bz*Bz;
    double pT = p + 0.5*B2;

    //find the sound speed
    double cs2 = gamma*p / rho;

    //find the fast speed
    double cfa2 = 0.5*(cs2 + B2/rho + std::sqrt((cs2 + B2/rho)*(cs2 + B2/rho) - 4.0*cs2*Bx*Bx/rho));

    return {cs2,cfa2};
}

std::tuple<double,double,double> wavespeedx(const std::array<double,9>&x,const std::array<double,9>&y,double gamma, double ch, double BxTilde){
    
    //collect our wave speeds squared from the previous function
    auto[cs2L,cfa2L] = speedsx(x,gamma,BxTilde);
    auto[cs2R,cfa2R] = speedsx(y,gamma,BxTilde);

    //cprimitive variables for other calculations
    auto primL = ConservativeToPrimitive(x,gamma);
    auto primR = ConservativeToPrimitive(y,gamma);

    //left and right wave speeds
    double sL = std::min(primL[1], primR[1]) - std::max(std::sqrt(cfa2L),std::sqrt(cfa2R));
    double sR = std::max(primL[1], primR[1]) + std::max(std::sqrt(cfa2L),std::sqrt(cfa2R));
    
    //magnetic magnitude, total pressure, modified pressure (subtractive normal magnetic component squared)
    double B2L = x[5]*x[5]+x[6]*x[6]+x[7]*x[7];
    double pTL = primL[4] + 0.5*B2L;
    double pModTL = pTL - x[5]*x[5];

    //right
    double B2R = y[5]*y[5]+y[6]*y[6]+y[7]*y[7];
    double pTR = primR[4] + 0.5*B2R;
    double pModTR = pTR - y[5]*y[5];

    //intermediate wavespeed qStar
    double qStar = (y[0]*primR[1]*(sR - primR[1]) - x[0]*primL[1]*(sL - primL[1]) + pModTL - pModTR) / (y[0]*(sR - primR[1]) - x[0]*(sL - primL[1]));

    return {sL, qStar, sR};
}

std::tuple<double,double> speedsy( const std::array<double , 9>&x , double gamma , double ByTilde){
    
    // collect primitive and conservative variables
    auto primX = ConservativeToPrimitive(x,gamma);
    double rho = x[0];
    double u_x = primX[1];
    double u_y = primX[2];
    double u_z = primX[3];
    double E = x[4];
    double p = primX[4];
    double By = x[6];
    double Bx = x[5];
    double Bz = x[7];
    
    //calculate magnetic magnitude and total pressure
    double B2 = Bx*Bx + By*By + Bz*Bz;
    double pT = p + 0.5*B2;

    //calculate sound speed
    double cs2 = gamma*p / rho;

    //find the fast speed
    double cfa2 = 0.5*(cs2 + B2/rho + std::sqrt((cs2 + B2/rho)*(cs2 + B2/rho) - 4.0*cs2*By*By/rho));

    return {cs2,cfa2};
}

std::tuple<double,double,double> wavespeedy(const std::array<double,9>x,const std::array<double,9>y,double gamma, double ch, double ByTilde){
   
    //collect the wave speeds squared from the previous function
    auto[cs2L,cfa2L] = speedsy(x,gamma, ByTilde);
    auto[cs2R,cfa2R] = speedsy(y,gamma, ByTilde);

    //get primitive variables for other calculations
    auto primL = ConservativeToPrimitive(x,gamma);
    auto primR = ConservativeToPrimitive(y,gamma);

    //left and right wave speeds
    double sL = std::min(primL[2], primR[2]) - std::max(std::sqrt(cfa2L),std::sqrt(cfa2R));
    double sR = std::max(primL[2], primR[2]) + std::max(std::sqrt(cfa2L),std::sqrt(cfa2R));
    
    //magnetic magnitudes and total pressure, modified pressure ( subtracting the normal magnetic field component) 
    //left
    double B2L = x[5]*x[5]+x[6]*x[6]+x[7]*x[7];
    double pTL = primL[4] + 0.5*B2L;
    double pModTL = pTL - x[6]*x[6];

    //right
    double B2R = y[5]*y[5]+y[6]*y[6]+y[7]*y[7];
    double pTR = primR[4] + 0.5*B2R;
    double pModTR = pTR - y[6]*y[6];

    //intermedaite wave speed qStar
    double qStar = (y[0]*primR[2]*(sR - primR[2]) - x[0]*primL[2]*(sL - primL[2]) + pModTL - pModTR) / (y[0]*(sR - primR[2]) - x[0]*(sL - primL[2]));

    return {sL, qStar, sR};
}

std::array<double , 9> uHLLx ( const std::array<double , 9>&x , const std::array<double , 9>& y , double gamma,double ch, double BxTilde){
    
    //get wave speeds for left riht and intermediate
    auto [sL, qStar, sR] = wavespeedx(x, y, gamma, ch, BxTilde);

    std::array<double , 9> uStar;
    auto FluxX = FluxDefX(x, gamma, ch);
    auto FluxY = FluxDefX(y, gamma,ch);

    //calculate the intemediate value U*
    for (int i = 0; i < 9; ++i) {
        uStar[i] = (sR * y[i] - sL * x[i] - FluxY[i] + FluxX[i]) / (sR - sL);
    }

    //use the wavespeeds to decide the intermediate value to take
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

std::array<double , 9> uHLLy ( const std::array<double , 9>&x , const std::array<double , 9>& y , double gamma, double ch, double ByTilde){
   
    //get wave speeds for left right and intermediate
    auto [sL, qStar, sR] = wavespeedy(x, y, gamma, ch, ByTilde);
    
    std::array<double , 9> uStar;
    auto FluxX = FluxDefY(x, gamma,ch);
    auto FluxY = FluxDefY(y, gamma,ch);

    //calculate the intermediate value U*
    for (int i = 0; i < 9; ++i) {
        uStar[i] = (sR * y[i] - sL * x[i] - FluxY[i] + FluxX[i]) / (sR - sL);
    }

    //decide which intermediate value to take as per the wave speeds
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

std::tuple<std::array<double , 9>, std::array<double , 9>> uStarX( const std::array<double , 9>&x , const std::array<double , 9>& y , double gamma , double ch, double BxTilde, double psiTilde){

    auto [sL, qStar, sR] = wavespeedx(x, y, gamma, ch, BxTilde);
    auto primL = ConservativeToPrimitive(x, gamma);
    auto primR = ConservativeToPrimitive(y, gamma);
    auto uHLL = uHLLx(x, y, gamma, ch, BxTilde);

    // Left state primitives and conserved variables
    double rhoL = x[0];
    double u_xL = primL[1], u_yL = primL[2], u_zL = primL[3];
    double EL = x[4];
    double pL = primL[4];
    double BxL = BxTilde, ByL = x[6], BzL = x[7];

    // Right state primitives and conserved variables
    double rhoR = y[0];
    double u_xR = primR[1], u_yR = primR[2], u_zR = primR[3];
    double ER = y[4];
    double pR = primR[4];
    double BxR = BxTilde, ByR = y[6], BzR = y[7];

    // HLL intermediate state 
    double u_xHLL = uHLL[1], u_yHLL = uHLL[2], u_zHLL = uHLL[3];
    double BxHLL = uHLL[5], ByHLL = uHLL[6], BzHLL = uHLL[7];

    // Magnetic field magnitudes and total pressure in left and right states
    double B2L = BxL*BxL + ByL*ByL + BzL*BzL;
    double B2R = BxR*BxR + ByR*ByR + BzR*BzR;
    double pTL = pL + 0.5 * B2L;
    double pTR = pR + 0.5 * B2R;

    // Compute intermediate normal magnetic field BxStar 
    double BxStar = (sR * BxR - sL * BxL) / (sR - sL);
    
    // Compute intermediate total pressures pStarL and pStarR 
    double pStarL = rhoL * (sL - u_xL) * (qStar - u_xL) + pTL - BxL*BxL + BxStar*BxStar;
    double pStarR = rhoR * (sR - u_xR) * (qStar - u_xR) + pTR - BxR*BxR + BxStar*BxStar;

    // Compute intermediate star density 
    double rhoStarL = rhoL * (sL - u_xL) / (sL - qStar);
    double rhoStarR = rhoR * (sR - u_xR) / (sR - qStar);

    //Compute intemediate momentum values
    double momxStarL = rhoStarL * qStar;                          
    double momyStarL = u_yL * rhoStarL - (BxStar * ByHLL - BxL * ByL) / (sL - qStar);  
    double momzStarL = u_zL * rhoStarL - (BxStar * BzHLL - BxL * BzL) / (sL - qStar); 
    double momxStarR = rhoStarR * qStar;                          
    double momyStarR = u_yR * rhoStarR - (BxStar * ByHLL - BxR * ByR) / (sR - qStar);  
    double momzStarR = u_zR * rhoStarR - (BxStar * BzHLL - BxR * BzR) / (sR - qStar); 

    // Dot products B.v for left, right, and HLL states (needed for energy and momentum)
    double vDotBL = BxL * u_xL + ByL * u_yL + BzL * u_zL;
    double vDotBR = BxR * u_xR + ByR * u_yR + BzR * u_zR;
    double vDotBHLLL = (BxStar * momxStarL + ByHLL * momyStarL + BzHLL * momzStarL) / rhoStarL;
    double vDotBHLLR = (BxStar * momxStarR + ByHLL * momyStarR + BzHLL * momzStarR) / rhoStarR;


    // Compute intermediate star energies 
    double energyStarL = EL * (rhoStarL / rhoL) + ((pStarL * qStar - pTL * u_xL) - (BxStar * vDotBHLLL - BxL * vDotBL)) / (sL - qStar);
    double energyStarR = ER * (rhoStarR / rhoR) + ((pStarR * qStar - pTR * u_xR) - (BxStar * vDotBHLLR - BxR * vDotBR)) / (sR - qStar);

    // Compose intermediate left star conserved state vector U*_L
    std::array<double, 9> uStarL;
    uStarL[0] = rhoStarL;
    uStarL[1] = momxStarL;                          
    uStarL[2] = momyStarL;  
    uStarL[3] = momzStarL;  
    uStarL[4] = energyStarL;
    uStarL[5] = BxStar;
    uStarL[6] = ByHLL;
    uStarL[7] = BzHLL;
    uStarL[8] = psiTilde;   

    // Compose intermediate right star conserved state vector U*_R
    std::array<double, 9> uStarR;
    uStarR[0] = rhoStarR;
    uStarR[1] = momxStarR;
    uStarR[2] = momyStarR;
    uStarR[3] = momzStarR;
    uStarR[4] = energyStarR;
    uStarR[5] = BxStar;
    uStarR[6] = ByHLL;
    uStarR[7] = BzHLL;
    uStarR[8] = psiTilde;   

    return {uStarL, uStarR};
}

std::tuple<std::array<double , 9>, std::array<double , 9>> uStarY( const std::array<double , 9>&x , const std::array<double , 9>& y , double gamma , double ch, double ByTilde, double psiTilde){
   
    auto [sL,qStar,sR] = wavespeedy(x,y,gamma,ch, ByTilde);
    auto primL = ConservativeToPrimitive(x,gamma);
    auto primR = ConservativeToPrimitive(y,gamma);
    auto uHLL = uHLLy(x,y,gamma,ch, ByTilde);

    // Left state primitive and conservative variables
    double rhoL = x[0];
    double u_xL = primL[1],u_yL = primL[2],u_zL = primL[3];
    double EL = x[4], pL = primL[4];
    double BxL = x[5], ByL=ByTilde, BzL=x[7];
    
    // Right state primitive and conservative variables
    double rhoR = y[0];
    double u_xR = primR[1],u_yR = primR[2],u_zR = primR[3];
    double ER = y[4], pR = primR[4];
    double BxR = y[5], ByR=ByTilde, BzR=y[7];
    
    // hLL intermediate states
    double u_xHLL = uHLL[1], u_yHLL = uHLL[2], u_zHLL = uHLL[3], BxHLL=uHLL[5], ByHLL = uHLL[6], BzHLL = uHLL[7];
    
    // Magnetic feild magnitudes and total pressure in left and right states
    double B2L = BxL*BxL + ByL*ByL + BzL*BzL;
    double B2R = BxR*BxR + ByR*ByR + BzR*BzR;
    double pTL = pL + 0.5*B2L;
    double pTR = pR + 0.5*B2R;

    // Intermediate normal magnetic field ByStar
    double ByStar = (sR*ByR - sL*ByL)/(sR-sL);

    //Intermediate total pressures pStarL and pStarR
    double pStarL = rhoL*(sL - u_yL)*(qStar - u_yL) + pTL - ByL*ByL + ByStar*ByStar;
    double pStarR = rhoR*(sR - u_yR)*(qStar - u_yR) + pTR - ByR*ByR + ByStar*ByStar;

    // Compute intermediate star densities
    double rhoStarL = rhoL*(sL - u_yL)/(sL - qStar);
    double rhoStarR = rhoR*(sR - u_yR)/(sR - qStar);

    //Compute intemediate momentum values
    double momxStarL = rhoStarL*u_xL - (BxHLL*ByStar - BxL*ByL)/(sL - qStar);                          
    double momyStarL = rhoStarL*qStar;  
    double momzStarL = rhoStarL*u_zL - (ByStar*BzHLL - ByL*BzL)/(sL - qStar); 
    double momxStarR = rhoStarR*u_xR - (BxHLL*ByStar - BxR*ByR)/(sR - qStar);                          
    double momyStarR = rhoStarR*qStar;  
    double momzStarR = rhoStarL*u_zL - (ByStar*BzHLL - ByL*BzL)/(sL - qStar); 

    //Dot products B.v for left and right, and HLL states (for the energy calculation)
    double vDotBL = BxL * u_xL + ByL * u_yL + BzL * u_zL;
    double vDotBR = BxR * u_xR + ByR * u_yR + BzR * u_zR;
    double vDotBHLLL = (BxHLL * momxStarL + ByStar * momyStarL + BzHLL * momzStarL) / rhoStarL;
    double vDotBHLLR = (BxHLL * momxStarR + ByStar * momyStarR + BzHLL * momzStarR) / rhoStarR;


    // Compute intermediate star energies
    double energyStarL = EL*(rhoStarL/rhoL) + ((pStarL*qStar - pTL*u_yL) - (ByStar*vDotBHLLL - ByL*vDotBL))/(sL - qStar);
    double energyStarR = ER*(rhoStarR/rhoR) + ((pStarR*qStar - pTR*u_yR) - (ByStar*vDotBHLLR - ByR*vDotBR))/(sR - qStar);

    // Compose intermediate left star conserved state U*_L
    std::array<double ,9> uStarL;
    uStarL[0] = rhoStarL;
    uStarL[1] = momxStarL;                          
    uStarL[2] = momyStarL;  
    uStarL[3] = momzStarL; 
    uStarL[4] = energyStarL;
    uStarL[5] = BxHLL;
    uStarL[6] = ByStar;
    uStarL[7] = BzHLL;
    uStarL[8] = psiTilde;

    // Compute intermediate right star conserved state U*_R
    std::array<double ,9> uStarR;
    uStarR[0] = rhoStarR;
    uStarR[1] = momxStarR;
    uStarR[2] = momyStarR;
    uStarR[3] = momzStarR;
    uStarR[4] = energyStarR;
    uStarR[5] = BxHLL;
    uStarR[6] = ByStar;
    uStarR[7] = BzHLL;
    uStarR[8] = psiTilde;   

    return {uStarL, uStarR};

}

std::array<double,9> FluxHLLCX(const std::array<double,9>& x, const std::array<double,9>& y, double gamma, double ch){
    auto [BxTilde, psiTilde] = psiAndBx(x,y,gamma,ch);
    
    auto [uStarL, uStarR] = uStarX(x,y,gamma,ch,BxTilde, psiTilde);
    auto [sL, qStar, sR] = wavespeedx(x, y, gamma, ch, BxTilde);

    //construct a new vector with cleaning to use in the flux calculations
    std::array<double,9> x1,y1;
    x1=x;
    y1=y;
    x1[5] = y1[5] = BxTilde;
    x1[8] = y1[8] = psiTilde;

    //calculate the flux with these values
    std::array<double , 9> FluxL = FluxDefX(x1,gamma,ch);
    std::array<double , 9> FluxR = FluxDefX(y1,gamma,ch);
    std::array<double , 9> FluxHLLC;

    //use the conditions from the wave speeds to decide which intermediate to take
    if (sL >= 0) {
        FluxHLLC = FluxL;
    } else if (sL <= 0 && qStar >= 0) {
        for (int i = 0; i < 9; ++i)
            FluxHLLC[i] = FluxL[i] + sL * (uStarL[i] - x1[i]);
    } else if (qStar <= 0 && sR >= 0) {
        for (int i = 0; i < 9; ++i)
            FluxHLLC[i] = FluxR[i] + sR * (uStarR[i] - y1[i]);
    } else {
        FluxHLLC = FluxR;
    }

    return FluxHLLC;
}

std::array<double,9> FluxHLLCY(const std::array<double,9> &x, const std::array<double,9> &y, double gamma, double ch){
    auto [ByTilde, psiTilde] = psiAndBy(x,y,gamma,ch); //start function by applying divergence cleaning
    
    auto [uStarL, uStarR] = uStarY(x,y,gamma,ch, ByTilde, psiTilde);
    auto [sL, qStar, sR] = wavespeedy(x, y, gamma, ch, ByTilde);

    //construct a new vector with the cleaning for divergence to use in the flux functions
    std::array<double,9> x1,y1;
    x1=x;
    y1=y;
    x1[6] = y1[6] = ByTilde;
    x1[8] = y1[8] = psiTilde;

    //then calulate the flux with the cleaned vectors
    std::array<double , 9> FluxL = FluxDefY(x1,gamma,ch);
    std::array<double , 9> FluxR = FluxDefY(y1,gamma,ch);
    std::array<double , 9> FluxHLLC;

    //use the conditions from the wave speeds to decide which intermediate values to take
    if (sL >= 0) {
        FluxHLLC = FluxL;
    } else if (sL <= 0 && qStar >= 0) {
        for (int i = 0; i < 9; ++i)
            FluxHLLC[i] = FluxL[i] + sL * (uStarL[i] - x1[i]);
    } else if (qStar <= 0 && sR >= 0) {
        for (int i = 0; i < 9; ++i)
            FluxHLLC[i] = FluxR[i] + sR * (uStarR[i] - y1[i]);
    } else {
        FluxHLLC = FluxR;
    }

    return FluxHLLC;
}

std::array<double , 9> FluxHLLX(const std::array<double , 9>&x , const std::array<double , 9> &y , double gamma,double ch){
    auto [BxTilde, psiTilde] = psiAndBx(x,y,gamma,ch);
    std::array<double,9> x1,y1;
    x1=x;
    y1=y;
    // x1[5] = y1[5] = BxTilde;
    // x1[8] = y1[8] = psiTilde;
    auto [uStarL, uStarR] = uStarX(x1,y1,gamma,ch, BxTilde, psiTilde);
    auto [sL, qStar, sR] = wavespeedx(x1, y1, gamma, ch, BxTilde);
    
    std::array<double , 9> Flux;
    if(sL >= 0){
        Flux = FluxDefX(x1,gamma,ch);
    }
    else if(sL < 0 && sR > 0){
        auto FluxX = FluxDefX(x1,gamma,ch);
        auto FluxY = FluxDefX(y1,gamma,ch);
        for(int i=0 ; i<9; ++i){
            Flux[i] = (sR*FluxX[i] - sL*FluxY[i] + sL*sR*(y1[i] - x1[i]))/(sR - sL);
        }
    }
    else { // sR <= 0
        Flux = FluxDefX(y1,gamma,ch);
    }
    
    return Flux;
}

std::array<double , 9> FluxHLLY(const std::array<double , 9>&x , const std::array<double , 9>& y , double gamma,double ch){
    auto [ByTilde, psiTilde] = psiAndBy(x,y,gamma,ch);
    std::array<double,9> x1,y1;
    x1=x;
    y1=y;
    // x1[6] = y1[6] = ByTilde;
    // x1[8] = y1[8] = psiTilde;
    auto [uStarL, uStarR] = uStarY(x1,y1,gamma,ch, ByTilde, psiTilde);
    auto [sL, qStar, sR] = wavespeedy(x1, y1, gamma, ch, ByTilde);
    
    std::array<double , 9> Flux;
    if(sL >= 0){
        Flux = FluxDefY(x1,gamma,ch);
    }
    else if(sL < 0 && sR > 0){
        auto FluxX = FluxDefY(x1,gamma,ch);
        auto FluxY = FluxDefY(y1,gamma,ch);
        for(int i=0 ; i<9; ++i){
            Flux[i] = (sR*FluxX[i] - sL*FluxY[i] + sL*sR*(y1[i] - x1[i]))/(sR - sL);
        }
    }
    else { // sR <= 0
        Flux = FluxDefY(y1,gamma,ch);
    }
    
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
    for (int i = 0; i < nxCells + 2; ++i) {//reflect in u_y and B_y
        u[i][0] = u[i][2];
        u[i][1] = u[i][3];        // Bottom boundary
        u[i][nyCells + 2] = u[i][nyCells];
        u[i][nyCells + 3] = u[i][nyCells + 1];  // Top boundary
        u[i][0][2] = -u[i][2][2]; 
        u[i][0][6] = -u[i][2][6]; 
        u[i][1][2] = -u[i][3][2]; 
        u[i][1][6] = -u[i][3][6]; 
        u[i][nyCells + 2][2] = -u[i][nyCells][2]; 
        u[i][nyCells + 2][6] = -u[i][nyCells][6]; 
        u[i][nyCells + 3][2] = -u[i][nyCells + 1][2]; 
        u[i][nyCells + 3][6] = -u[i][nyCells + 1][6]; 
    }

    //bottom and top periodic
    // for (int i = 0; i < nxCells + 2; ++i) {
    //     u[i][0] = u[i][nyCells];      // Bottom boundary
    //     u[i][1] = u[i][nyCells+1];    
    //     u[i][nyCells + 2] = u[i][2];  // Top boundary
    //     u[i][nyCells + 3] = u[i][3];  
    // }

    //bottom and top transmissive
    // for (int i = 0; i < nxCells + 2; ++i) {
    //     u[i][0] = u[i][2];      // Bottom boundary
    //     u[i][1] = u[i][3];    
    //     u[i][nyCells + 2] = u[i][nyCells];  // Top boundary
    //     u[i][nyCells + 3] = u[i][nyCells + 1];  
    // }
}


int main() { 
    int nxCells = 128; 
    int nyCells = 256;
    double x0 = 0.0;
    double x1 = 1.0;
    double y0 = -1.0;
    double y1 = 1.0;
    double tStart = 0.0; //set the start and finish time steps the same
    double tStop = 12.0;
    double C = 0.8;
    double gamma = 5.0 / 3.0;
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

    double time;

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
            prim[0] = 1.0; // Density
            prim[1] = 0.5*std::tanh(20*y); // Velocity
            prim[2] = 0.01*std::sin(2.0*pi*x)*std::exp(-y*y/(0.01));
            prim[3] = 0.0;
            prim[4] = 1.0 / gamma; // pressure
            prim[5] = 0.1*std::cos(pi/3.0); // magnetic field
            prim[6] = 0.0;
            prim[7] = 0.1*std::sin(pi/3.0); 
            prim[8] = 0;

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

    auto [dt, maxSpeed] = ComputeTimeStep(u , C , dx, dy, gamma); //the time steps

    double t = tStart;
    int counter =0;
    do {

        auto [dt,maxSpeed] = ComputeTimeStep(u, C, dx, dy, gamma); 
        double ch = maxSpeed;
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
                std::array<double,9> fluxLx = FluxDefX(uBarL_x[i][j], gamma,ch);
                std::array<double,9> fluxRx = FluxDefX(uBarR_x[i][j], gamma,ch);
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

        for(int j = 0; j < nyCells+4; j++){//update psi inbetween
            for(int i = 0; i < nxCells+4; i++){
                uBar[i][j][8] = psiUpdate(uBar[i][j][8], maxSpeed, 0.5*dt);
            }
        }

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
                std::array<double,9> fluxBy = FluxDefY(uBarB_y[i][j], gamma,ch);
                std::array<double,9> fluxTy = FluxDefY(uBarT_y[i][j], gamma,ch);
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

            applyBoundaryConditions(u , nxCells , nyCells);

        //convert it back to primitive and define final results

        if(t > 0.1*counter && t-dt <= 0.1*counter){

            std::vector<std::vector<std::array<double,9>> >results;
            results.resize(nxCells+4, std::vector<std::array<double, 9> >(nyCells + 4)); //set up results

            for(int j = 0; j < nyCells+4; j++) { 
                for(int i = 0; i < nxCells+4; i++) {
                    results[i][j] = ConservativeToPrimitive(u[i][j], gamma);
                }
            }




            // Output the results
            std::string filename = "MHD" + std::to_string(counter) + ".dat";
            std::ofstream output(filename);
            for(int j = 1; j < nyCells+4; j++) { 
                for(int i = 1; i < nxCells+4; i++) {
                    double x = x0 + (i - 1)*dx;
                    double y = y0 + (j - 1)*dy;
                    output << x << " " << y << " " << results[i][j][0] << " " << results[i][j][1] << " " << results[i][j][2] << " " << results[i][j][3] <<  " " << results[i][j][4]<< " " << std::sqrt(results[i][j][5]*results[i][j][5] + results[i][j][6]*results[i][j][6]) / results[i][j][7] <<" " << results[i][j][6] <<  " " << results[i][j][7] << " " << results[i][j][8]<<std::endl;
                }
                output<<std::endl;
            }
            counter +=1;
        }

    } while (t < tStop );


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
    std::string filename = "MHD" + std::to_string(counter) + ".dat";
    std::ofstream output(filename);
    for(int j = 1; j < nyCells+4; j++) { 
        for(int i = 1; i < nxCells+4; i++) {
            double x = x0 + (i - 1)*dx;
            double y = y0 + (j - 1)*dy;
            output << x << " " << y << " " << results[i][j][0] << " " << results[i][j][1] << " " << results[i][j][2] << " " << results[i][j][3] <<  " " << results[i][j][4]<< " " << std::sqrt(results[i][j][5]*results[i][j][5] + results[i][j][6]*results[i][j][6]) / results[i][j][7] <<" " << results[i][j][6] <<  " " << results[i][j][7] << " " << results[i][j][8]<<std::endl;
        }
        output<<std::endl;
    }
    
}