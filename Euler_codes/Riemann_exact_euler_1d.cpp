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

double NewtonRaphson(double rhoL , double vL , double pL,double rhoR , double vR , double pR, double pOld ,double gamma){
    double AL = 2.0 / ((gamma + 1.0)*rhoL);
    double AR = 2.0 / ((gamma + 1.0)*rhoR);

    double BL = ((gamma-1.0) * pL)/ (gamma+1.0);
    double BR = ((gamma-1.0) * pR)/ (gamma+1.0);


    double csL = std::sqrt(gamma * pL / rhoL);
    double csR = std::sqrt(gamma * pR / rhoR);

    double fL , fR , fPrimeL , fPrimeR;
    if(pL < pOld){//nrshock for fl
        fL = (pOld - pL)*std::sqrt(AL / (pOld + BL));
        fPrimeL = std::sqrt(AL / (BL + pOld))*(1-((pOld - pL)/(2*(BL+pOld))));
    }
    else{//nr refraction
        fL = (2.0*csL/(gamma -1)) * ( std::pow((pOld / pL) , ((gamma-1)/(2.0*gamma))) - 1);
        fPrimeL = (1.0 / (csL * rhoL)) * std::pow((pOld / pL), (-(gamma + 1) / (2.0 * gamma)));
    }

    if(pR < pOld){//shock for fr
        fR = (pOld - pR)*std::sqrt(AR / (pOld + BR));
        fPrimeR = std::sqrt(AR / (BR + pOld))*(1-((pOld - pR)/(2*(BR+pOld))));
    }
    else{//nr refraction
        fR = (2.0*csR/(gamma -1)) * ( std::pow((pOld / pR) , ((gamma-1)/(2.0*gamma))) - 1);
        fPrimeR = (1.0 / (csR * rhoR)) * std::pow((pOld / pR), (-(gamma + 1) / (2.0 * gamma)));
    }

    double f = fR + fL + vR - vL;
    double fPrime = fPrimeR + fPrimeL;
    double pNew = pOld - (f / fPrime);

    return pNew;
}

double FindpStar(std::array<double , 3>uL , std::array<double , 3>uR , double gamma , double epsilon=1e-6){
    double rhoL = uL[0];
    double vL = uL[1] ;
    double pL = uL[2];
    double rhoR = uR[0];
    double vR = uR[1];
    double pR = uR[2];

    double AL = 2.0 / ((gamma + 1.0)*rhoL);
    double AR = 2.0 / ((gamma + 1.0)*rhoR);

    double BL = ((gamma-1.0) * pL)/ (gamma+1.0);
    double BR = ((gamma-1.0) * pR)/ (gamma+1.0);

    double csL = std::sqrt(gamma * pL / rhoL);
    double csR = std::sqrt(gamma * pR / rhoR);

    double ppL = csL / std::pow(pL , (gamma - 1)/(2.0*gamma));
    double ppR = csR / std::pow(pR , (gamma - 1)/(2.0*gamma));

    

    //need to start the process with a guess for pressure (from toros book)

    double p_pvrs = 0.5*(pL+pR) - 0.125*(vR -vL)*(rhoL + rhoR)*(csL + csR);
    double pOld;
    double pNew = std::max(1e-10 , p_pvrs);
    //pNew = 0.5*(pL+pR);
    pNew = std::pow(((csL + csR - 0.5*(gamma-1)*(vR-vL))/((ppL) + (ppR))) , 2.0*gamma/(gamma-1));
    //std::cout << pL <<" "<<pR<<" "<< vL << " " << vR << " "<< rhoL<< " "<<rhoR << " "<<csL << " "  <<csR <<" pNew = "<< pNew<< std::endl;

    int i =0;
    do{
        pOld = pNew;
        
        pNew = NewtonRaphson(rhoL , vL , pL , rhoR , vR , pR , pOld , gamma);
        i+=1;
        //std::cout << "pOld = " << pOld << " pNew = "<< pNew << " i= "<< i << " difference = " << std::abs(2.0*(pOld - pNew)/(pOld + pNew)) << std::endl;

    }while (std::abs(2.0*(pOld - pNew)/(pOld + pNew))> epsilon);
    //std::cout << "p* = "<< pNew << std::endl;
    return pNew;
}

std::array<double , 3> FindvStarrhoStars(std::array<double , 3>uL , std::array<double , 3>uR , double gamma , double epsilon=1e-6){
    double pStar = FindpStar(uL , uR , gamma);
    double rhoL = uL[0];
    double vL = uL[1] ;
    double pL = uL[2];
    double rhoR = uR[0];
    double vR = uR[1];
    double pR = uR[2];

    std::array<double , 3> answers;

    double AL = 2.0 / ((gamma + 1.0)*rhoL);
    double AR = 2.0 / ((gamma + 1.0)*rhoR);

    double BL = ((gamma-1.0) * pL)/ (gamma+1.0);
    double BR = ((gamma-1.0) * pR)/ (gamma+1.0);

    double csL = std::sqrt(gamma * pL / rhoL);
    double csR = std::sqrt(gamma * pR / rhoR);

    double rhoStarL, rhoStarR, vStar , fL , fR , sL , sR , vStarL, vStarR;
    int choice;

    if(pL < pStar){//nrshock for fl
        fL = (pStar - pL)*std::sqrt(AL / (pStar + BL));
        rhoStarL = (rhoL*((pStar/pL) + (gamma-1)/(gamma+1)))/(((gamma-1)/(gamma+1))*(pStar/pL)+1);
    }
    else{//rarefaction
        rhoStarL = rhoL*std::pow((pStar/pL) , 1/gamma);
        fL = (2.0*csL/(gamma -1)) * ( std::pow((pStar/ pL) , ((gamma-1)/(2*gamma))) - 1);
    }
    if(pR < pStar){//nrshock for fl
        fR = (pStar - pR)*std::sqrt(AR / (pStar + BR));
        rhoStarR = (rhoR*((pStar/pR) + (gamma-1)/(gamma+1)))/(((gamma-1)/(gamma+1))*(pStar/pR)+1);
    }
    else{//rarefaction
        rhoStarR = rhoR*std::pow((pStar/pR) , 1/gamma);
        fR = (2.0*csR/(gamma -1)) * ( std::pow((pStar/ pR) , ((gamma-1)/(2*gamma))) - 1);
    }

    vStar = 0.5*(vL+vR) + 0.5*(fR - fL);
    answers[0] = vStar;
    answers[1] = rhoStarL;
    answers[2] = rhoStarR;

    return answers;
}

std::vector<std::array<double,3>> reconstructStates(std::vector<std::array<double,3>>v, double gamma , double nCells){
    std::vector<std::array<double,3>> uStar(nCells+2);
    std::vector<std::array<double,3>> u(nCells+2);
    std::vector<std::array<double,3>> answer(nCells+2);
    std::vector<double> pStar(nCells + 2);
    std::vector<double> v_Star(nCells + 2);
    std::vector<double> rho_StarL(nCells + 2);
    std::vector<double> rho_StarR(nCells + 2);
    std::vector<double> rho_Star(nCells + 2);
    
    for(int j=0 ; j<=nCells+1 ; j++){
        u[j] = conservativeToPrimative(v[j], gamma);
    }

    for(int i = 0; i < nCells+1; i++) { 
        pStar[i] = FindpStar(u[i] , u[i+1] ,  gamma );
        std::array<double , 3> vrhoStars = FindvStarrhoStars(u[i] , u[i+1] ,  gamma );
        v_Star[i] = vrhoStars[0];
        rho_StarL[i] = vrhoStars[1];
        rho_StarR[i] = vrhoStars[2];
        double rhoL = u[i][0];
        double vL = u[i][1] ;
        double pL = u[i][2];
        double rhoR = u[i+1][0];
        double vR = u[i+1][1];
        double pR = u[i+1][2];
        double csL = std::sqrt(gamma * pL / rhoL);
        double csR = std::sqrt(gamma * pR / rhoR);
        
        double sHL = vL - csL;
        double sTL = v_Star[i] - csL*std::pow((pStar[i]/pL), ((gamma-1)/2.0*gamma));
        double sL = vL - csL*std::sqrt(((gamma+1)/(2.0*gamma))*(pStar[i]/pL)+(gamma-1)/(2.0*gamma));

        double sHR = vR + csR;
        double sTR = v_Star[i] + csR*std::pow((pStar[i]/pR), ((gamma-1)/2.0*gamma));
        double sR = vR + csR*std::sqrt(((gamma+1)/(2.0*gamma))*(pStar[i]/pR)+(gamma-1)/(2.0*gamma));

        std::array<double , 3> WLfan;
        WLfan[0] = rhoL*std::pow((2.0/(gamma+1) + ((gamma-1)/((gamma+1)*csL))*vL), 2.0/(gamma-1));//rho
        WLfan[1] = (2.0/(gamma+1))*(csL + ((gamma-1)/2.0)*vL);//velocity
        WLfan[2] = pL*std::pow((2.0/(gamma+1) + ((gamma-1)/((gamma+1)*csL))*vL), 2.0/(gamma-1));//pressure

        std::array<double , 3> WRfan;
        WRfan[0] = rhoR*std::pow((2.0/(gamma+1) - ((gamma-1)/((gamma+1)*csR))*vR), 2.0/(gamma-1));//rho
        WRfan[1] = (2.0/(gamma+1))*(-csR + ((gamma-1)/2.0)*vR);//velocity
        WRfan[2] = pR*std::pow((2.0/(gamma+1) - ((gamma-1)/((gamma+1)*csR))*vR), 2.0/(gamma-1));//pressure

        std::array<double , 3> WStarL;
        WStarL[0] = rho_StarL[i];
        WStarL[1] = v_Star[i];
        WStarL[2] = pStar[i];

        std::array<double , 3> WStarR;
        WStarR[0] = rho_StarR[i];
        WStarR[1] = v_Star[i];
        WStarR[2] = pStar[i];

        std::array<double , 3> WL;
        WL[0] = rhoL;
        WL[1] = vL;
        WL[2] = pL;

        std::array<double , 3> WR;
        WR[0] = rhoR;
        WR[1] = vR;
        WR[2] = pR;


        if (v_Star[i] >= 0) {
            // star region coming from the left
            if(pStar[i]>pL){
                if(0<sL){
                    uStar[i] = WL;
                }
                else{
                    uStar[i] = WStarL;
                }
            }
            else{
                if(0<sHL){
                    uStar[i] = WL;
                }
                else{
                    if(0>sTL){
                        uStar[i] = WStarL;
                    }
                    else{
                        uStar[i] = WLfan;
                    }
                }
            }
        } else {
            if(pStar[i] > pR){
                if(sR >= 0){
                    uStar[i] = WStarR;
                }
                else{
                    uStar[i] = WR;
                }
            }
            else{
                if(0<=sTR){
                    uStar[i] = WStarR;
                }
                else{
                    if(0<sHR){
                        uStar[i] = WRfan;
                    }
                    else{
                        uStar[i] = WR;
                    }
                }
            }
        }
    }

    for(int j=0 ; j<=nCells+1 ; j++){
        answer[j] = primitiveToConservative(uStar[j], gamma);
    }

    return answer;
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


int main(){
    int nCells = 100; //the distance between points is 0.01
    double x0 = 0.0;
    double x1 = 1.0;
    double tStart = 0.0; //set the start and finish time steps the same
    double tStop = 0.15;
    double C = 1.0;
    double gamma = 1.4;

    // Allocate matrices with 2 extra points for transmissive BCs
    std::vector<std::array<double,3>> u(nCells+2);
    std::vector<std::array<double,3>> uPlus1(nCells+2);
    std::vector<std::array<double,3>> uPlusHalf(nCells+2);
    std::vector<std::array<double,3>> uStar(nCells+2);
    std::vector<std::array<double,3>> flux(nCells+2);
    std::vector<double> pStar(nCells + 2);
    std::vector<double> v_Star(nCells + 2);
    std::vector<double> rho_StarL(nCells + 2);
    std::vector<double> rho_StarR(nCells + 2);
    std::vector<double> rho_Star(nCells + 2);

    double dx = (x1 - x0) / nCells; //the space steps 
    double dt;

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

    double t = tStart;
    do {

        dt = computeTimeStep(u , C , dx , gamma); 
        t = t + dt;

        // Trasmissive boundary conditions
        u[0] = u[1];
        u[nCells + 1] = u[nCells];

        uPlusHalf = reconstructStates(u , gamma , nCells);
        //Trasmissive boundary conditions
        uPlusHalf[0] = uPlusHalf[1];
        uPlusHalf[nCells + 1] = uPlusHalf[nCells];

        for(int i = 0; i <= nCells+1; i++) { 
            flux[i] = flux_def(uPlusHalf[i] , gamma); //gives f(u_{i+1/2})
        }
        for(int i = 1; i <= nCells+1; i++) { //Update the data
            for(int j =0; j<3; ++j){
                uPlus1[i][j] = u[i][j] - (dt/dx) * (flux[i][j] - flux[i-1][j]);
            }
        }

        // Trasmissive boundary conditions
        uPlus1[0] = uPlus1[1];
        uPlus1[nCells + 1] = uPlus1[nCells];

        u = uPlus1;
    } while (t < tStop);

    //define final results

    std::vector<std::array<double,3>> results(nCells+2);

    for(int i=0; i<= results.size() -1; ++i){
        results[i] = conservativeToPrimative(u[i], gamma);
    }


    // Output the results
    std::ofstream output("exact_solver.dat");
    for (int i = 1; i < nCells+1; ++i) {
        double x = x0 + (i - 1) * dx;
        output << x << " " << results[i][0] <<  " " << results[i][1] <<  " " << results[i][2] << std::endl;
         //std::cout << x << " " << u[i][0] <<  " " << u[i][1] <<  " " << u[i][2] << std::endl;
    }
}
