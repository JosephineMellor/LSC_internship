#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <array>
#include <filesystem>
#include <tuple>

//defining some useful types
typedef std::array<double , 500> data_vec;
typedef std::vector<double> big_vector;
typedef std::vector<std::vector<double>> data_table;
typedef std::array<double,4> array;

//function to read in the data from the tabulated equation of state
std::tuple<data_vec , data_vec , data_table , data_table , data_table , data_table , data_table> Plasma19(){
    
    //load the file
    std::ifstream file("Plasma19_EoS.txt"); 

    //set up the lines, starting row, and data arrays
    std::string line;
    int row = 0;
    data_vec densities;
    data_vec pressures;

    //we want to save sound speeds, energies, temperatures, electrical conductivity and thermal conductivity in tables
    data_table sound_speeds(500, std::vector<double>(500));
    data_table energies(500, std::vector<double>(500));
    data_table temperatures(500, std::vector<double>(500));
    data_table electrical_conductivity(500, std::vector<double>(500));
    data_table thermal_conductivity(500, std::vector<double>(500));
    
    //now loop over the data file row by row, we only want the fifth and sixth row for densities and pressures
    while(std::getline(file , line)){
        row++;

        if( row ==10 || row ==12 ){
            std::istringstream iss(line); //extract the line

            big_vector row_data; 
            double val;
            while (iss >> val) {
                row_data.push_back(val);
            } //put it temporarily in an array

            //then save the densities and pressures used for indexing the tables

            for(int i=0; i<500; ++i){
                if(row ==10){
                    densities[i] = row_data[i];
                }
                if(row ==12){
                    pressures[i] = row_data[i];
                }
            }
        }

        //sound speed
        if( row >= 14 && row< 514){
            std::istringstream iss(line); 

            big_vector row_data; 
            double val;
            while (iss >> val) {
                row_data.push_back(val);
            }

            sound_speeds[row-14] = std::move(row_data);
        }

        //internal energy
        if( row >= 515 && row< 1015){
            std::istringstream iss(line); 

            big_vector row_data; 
            double val;
            while (iss >> val) {
                row_data.push_back(val);
            }
            if(row_data.size() != 500) {
                std::cerr << "Row " << row << " for energy has " << row_data.size() << " values, expected 500.\n";
                exit(1);
            }

            energies[row-515] = std::move(row_data);
        }

        //temperatures
        if( row >= 1016 && row< 1516){
            std::istringstream iss(line);

            big_vector row_data; 
            double val;
            while (iss >> val) {
                row_data.push_back(val);
            }

            temperatures[row-1016] = std::move(row_data);
        }

        //electrical conductivity
        if( row >= 1517 && row< 2017){
            std::istringstream iss(line); 

            big_vector row_data; 
            double val;
            while (iss >> val) {
                row_data.push_back(val);
            }

            electrical_conductivity[row-1517] = std::move(row_data);
        }

        //thermal conductivity
        if( row >= 2018 && row< 2518){
            std::istringstream iss(line); 

            big_vector row_data; 
            double val;
            while (iss >> val) {
                row_data.push_back(val);
            }

            thermal_conductivity[row-2018] = std::move(row_data);
        }
    }


    return{densities , pressures , sound_speeds , energies , temperatures , electrical_conductivity , thermal_conductivity};
}

//import our data tables and arrays as global variables

auto EoS_data = Plasma19(); // the tuple
const auto& densities = std::get<0>(EoS_data);
const auto& pressures = std::get<1>(EoS_data);
const auto& sound_speeds = std::get<2>(EoS_data);
const auto& energies = std::get<3>(EoS_data);
const auto& temperatures = std::get<4>(EoS_data);
const auto& electrical_conductivity = std::get<5>(EoS_data);
const auto& thermal_conductivity = std::get<6>(EoS_data);

//function to perform bilinear interpolation on a simplified grid (0,0) , (0,1) , (1,0) , (1,1)
double BilinearInterpolation(double f_NE, double f_SE, double f_SW, double f_NW, double x, double y){
    double f_x_y1 = (1-x)*f_SW + x*f_SE;
    double f_x_y2 = (1-x)*f_NW + x*f_NE;
    double f_x_y = (1-y)*f_x_y1 + y*f_x_y2;

    return f_x_y;
}

//function to 'reverse' bilinear interpolation on a simplified grid (0,0) , (0,1) , (1,0) , (1,1)
double ReverseBilinearInterpolation(double f_NE, double f_SE, double f_SW, double f_NW, double x, double f_x_y){
    double f_x_y1 = (1-x)*f_SW + x*f_SE;
    double f_x_y2 = (1-x)*f_NW + x*f_NE;
    double y = (f_x_y - f_x_y1) / (f_x_y2 - f_x_y1);
    if(f_x_y2 == f_x_y1){y = 1;}

    return y;
}


//functions to convert between primitive to conseravative variables
array PrimativeToConservative(array prim){
    array consv;

    //start with the variables that dont require the equation of state
    consv[0] = prim[0]; //rho
    consv[1] = prim[0] * prim[1]; //mom_x
    consv[2] = prim[0] * prim[2]; //mom_y

    //find index for density with binary search
    int point1 =0;
    int point2 =499;
    int half;
    int density_index;
    double rho = prim[0];
    while(std::abs(point1-point2)>1){
        half = floor((point1+point2)/2);
        if(rho > densities[half]){
            point1 = half;
        }
        else{
            point2 = half;
        }
    }

    //handle boundary terms
    if(point1 >= 499){point1 = 499;}
    if(point2 >= 499){point2 = 499;}
    if(point1 <=0){point1 = 0;}
    if(point2 <=0){point2 = 0;}

    //set up the linear interpolation
    double density_top = densities[point2];
    int density_top_i = point2;
    double density_bottom = densities[point1];
    int density_bottom_i = point1;
    double density_ratio = (rho - density_bottom) / (density_top - density_bottom);
    if(density_bottom == density_top){density_ratio =1;}
    

    //find index for pressure also with binary search
    point1 =0;
    point2 =499;
    half;
    int pressure_index;
    double p = prim[3];
    while(std::abs(point1-point2)>1){
        half = floor((point1+point2)/2);
        if(p > pressures[half]){
            point1 = half;
        }
        else{
            point2 = half;
        }
    }

    //handle boundary terms
    if(point1 >= 499){point1 = 499;}
    if(point2 >= 499){point2 = 499;}
    if(point1 <=0){point1 = 0;}
    if(point2 <=0){point2 = 0;}

    //set up linear interpolation
    double pressure_top = pressures[point2];
    int pressure_top_i = point2;
    double pressure_bottom = pressures[point1];
    int pressure_bottom_i = point1;
    double pressure_ratio = (p - pressure_bottom) / (pressure_top - pressure_bottom);
    if(pressure_bottom == pressure_top){pressure_ratio =1;}

    //then apply to the energies table
    double energies_NE = energies[pressure_top_i][density_top_i];
    double energies_SE = energies[pressure_bottom_i][density_top_i];
    double energies_SW = energies[pressure_bottom_i][density_bottom_i];
    double energies_NW = energies[pressure_top_i][density_bottom_i];

    //linear interpolate using our variables from before
    double energy = BilinearInterpolation(energies_NE , energies_SE , energies_SW , energies_NW , density_ratio , pressure_ratio);    

    consv[3] = energy*rho + 0.5*rho*(prim[1]*prim[1] + prim[2]*prim[2]);
    return consv;
}
array ConservativeToPrimative(array consv){
    array prim;

    //find the simpler variables 
    prim[0] = consv[0]; //rho
    prim[1] = consv[1] / consv[0]; //u_x
    prim[2] = consv[2] / consv[0]; //u_y

    //find our index for density with binary search
    int point1 =0;
    int point2 =499;
    int half;
    int density_index;
    double rho = prim[0];
    while(std::abs(point1-point2)>1){
        half = floor((point1+point2)/2);
        if(rho > densities[half]){
            point1 = half;
        }
        else{
            point2 = half;
        }
    }

    //handle boundary terms
    if(point1 >= 499){point1 = 499;}
    if(point2 >= 499){point2 = 499;}
    if(point1 <=0){point1 = 0;}
    if(point2 <=0){point2 = 0;}

    //set up the linear interpolation
    double density_top = densities[point2];
    int density_top_i = point2;
    double density_bottom = densities[point1];
    int density_bottom_i = point1;
    double density_ratio = (rho - density_bottom) / (density_top - density_bottom);
    if(density_bottom == density_top){density_ratio = 1;}

    //find our indices for pressure with binary search on the density-based linear interpolation of energies
    point1 =0;
    point2 =499;
    half;
    int pressure_index;
    double e = (consv[3] - 0.5*rho*(prim[1]*prim[1] + prim[2]*prim[2])) / rho; //internal energy

    //find the indices either side of our energy value in density top and bottom, they should be the same
    while(std::abs(point1-point2)>1){
        half = floor((point1+point2)/2);
        if(e > (energies[half][density_bottom_i] + density_ratio*(energies[half][density_top_i] - energies[half][density_bottom_i]))){
            point1 = half;
        }
        else{
            point2 = half;
        }
    }
    int pressure_top_i = point2;
    int pressure_bottom_i = point1;

    //handle boundary terms
    if(point1 >= 499){point1 = 499;}
    if(point2 >= 499){point2 = 499;}
    if(point1 <=0){point1 = 0;}
    if(point2 <=0){point2 = 0;}

    //set up our quadrant to search in
    double energies_NE = energies[pressure_top_i][density_top_i];
    double energies_SE = energies[pressure_bottom_i][density_top_i];
    double energies_SW = energies[pressure_bottom_i][density_bottom_i];
    double energies_NW = energies[pressure_top_i][density_bottom_i];

    //reverse bilinear interpolation
    double pressure_ratio = ReverseBilinearInterpolation(energies_NE , energies_SE , energies_SW , energies_NW , density_ratio , e); 
    
    //calculate the value of pressure
    double p = pressure_ratio * (pressures[pressure_top_i] - pressures[pressure_bottom_i]) + pressures[pressure_bottom_i];

    prim[3] = p;

    return prim;
}

//function to find the sound speed with bilinear interpolation, assumming you know the conservative form
double SoundSpeed(array u){
    array prim = ConservativeToPrimative(u);

    //find our index for density with binary search
    int point1 =0;
    int point2 =499;
    int half;
    int density_index;
    double rho = u[0]; 
    while(std::abs(point1-point2)>1){
        half = floor((point1+point2)/2);
        if(rho > densities[half]){
            point1 = half;
        }
        else{
            point2 = half;
        }
    }

    //handle boundary terms
    if(point1 >= 499){point1 = 499;}
    if(point2 >= 499){point2 = 499;}
    if(point1 <=0){point1 = 0;}
    if(point2 <=0){point2 = 0;}
    
    //set up the linear interpolation
    double density_top = densities[point2];
    int density_top_i = point2;
    double density_bottom = densities[point1];
    int density_bottom_i = point1;
    double density_ratio = (rho - density_bottom) / (density_top - density_bottom);
    if(density_bottom == density_top){density_ratio = 1;}
    

    //find our index for pressure also with linear bisection
    point1 =0;
    point2 =500;
    half;
    int pressure_index;
    double p = prim[3];
    while(std::abs(point1-point2)>1){
        half = floor((point1+point2)/2);
        if(p > pressures[half]){
            point1 = half;
        }
        else{
            point2 = half;
        }
    }

    //handle boundary terms
    if(point1 >= 499){point1 = 499;}
    if(point2 >= 499){point2 = 499;}
    if(point1 <=0){point1 = 0;}
    if(point2 <=0){point2 = 0;}

    //set up linear interpolation
    double pressure_top = pressures[point2];
    int pressure_top_i = point2;
    double pressure_bottom = pressures[point1];
    int pressure_bottom_i = point1;
    double pressure_ratio = (p - pressure_bottom) / (pressure_top - pressure_bottom);
    if(pressure_bottom == pressure_top){pressure_ratio =1;}

    //then apply to the energies table
    double ss_NE = sound_speeds[pressure_top_i][density_top_i];
    double ss_SE = sound_speeds[pressure_bottom_i][density_top_i];
    double ss_SW = sound_speeds[pressure_bottom_i][density_bottom_i];
    double ss_NW = sound_speeds[pressure_top_i][density_bottom_i];

    //linear interpolate using our variables from before
    double sound_speed = BilinearInterpolation(ss_NE , ss_SE , ss_SW , ss_NW , density_ratio , pressure_ratio);  
    
    return sound_speed;
}

//define the flux function but now we split into two functions fluxX and fluxY

array fluxX_def(std::array <double ,4> x){
    array fluxX;
    array prim = ConservativeToPrimative(x);
    double rho = prim[0];
    double v_x = prim[1];
    double v_y = prim[2];
    double p = prim[3];

    fluxX[0] = x[1];
    fluxX[1] = rho*v_x*v_x + p;
    fluxX[2] = rho * v_x * v_y;
    fluxX[3] = (x[3] + p)*v_x;

    return fluxX;
} //flux for the x-direction

array fluxY_def(std::array <double ,4> x){
    array fluxY;
    array prim = ConservativeToPrimative(x);
    double rho = prim[0];
    double v_x = prim[1];
    double v_y = prim[2];
    double p = prim[3];

    fluxY[0] = x[2];
    fluxY[1] = rho*v_x*v_y;
    fluxY[2] = rho * v_y * v_y + p;
    fluxY[3] = (x[3] + p)*v_y;

    return fluxY;
} //flux for the y-direction


array getFluxX(array x , array y , double dx , double dt){
    //find the flux at u_{i} and u_{i+1}

    array f_1 = fluxX_def(x); //f(ui)
    array f_2 = fluxX_def(y); //f(ui+1)

    //find u_{i+1/2}

    array uPlusHalf;
    for(int i=0; i<=3; ++i){
        uPlusHalf[i] = 0.5*(x[i] + y[i]) - 0.5*(dt/dx)*(f_2[i] - f_1[i]);
    }

    array RI_flux = fluxX_def(uPlusHalf); //richtmyer flux
    array LF_flux; 
    array FORCE_flux;
    
    for(int i=0; i<=3; ++i){
        LF_flux[i] = (0.5*dx/dt)*(x[i] - y[i]) + 0.5*(f_1[i] + f_2[i]); // lax friedrichs flux
        FORCE_flux[i] = 0.5 * (LF_flux[i] + RI_flux[i]); //force flux
    }

    return FORCE_flux;

    
} //FORCE in x


array getFluxY(array x , array y , double dy , double dt){
    //find the flux at u_{i} and u_{i+1}

    array f_1 = fluxY_def(x); //f(ui)
    array f_2 = fluxY_def(y); //f(ui+1)

    //find u_{i+1/2}

    array uPlusHalf;
    for(int i=0; i<=3; ++i){
        uPlusHalf[i] = 0.5*(x[i] + y[i]) - 0.5*(dt/dy)*(f_2[i] - f_1[i]);
    }

    array RI_flux = fluxY_def(uPlusHalf); //richtmyer flux
    array LF_flux; 
    array FORCE_flux;
    
    for(int i=0; i<=3; ++i){
        LF_flux[i] = (0.5*dy/dt)*(x[i] - y[i]) + 0.5*(f_1[i] + f_2[i]); // lax friedrichs flux
        FORCE_flux[i] = 0.5 * (LF_flux[i] + RI_flux[i]); //force flux
    }

    return FORCE_flux;
    
} //FORCE in y

double ComputeTimeStep(const std::vector<std::vector<std::array <double , 4>>>& u , double C , double dx , double dy){
    double dspace = dx;
    if(dy < dx){
        dspace = dy;
    }

    double maxSpeed = 0.0;


    for(int i=2 ; i<=u.size()-2 ; ++i){
        for(int j=2; j<=u[0].size()-2; ++j){
            array prim = ConservativeToPrimative(u[i][j]);
            double rho = prim[0];
            double v_x = prim[1];
            double v_y = prim[2];
            double pressure = prim[3];
            double speed = std::sqrt(v_x*v_x + v_y*v_y);

            double sound_speed = SoundSpeed(u[i][j]);
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



std::vector<std::vector<std::array<double, 4> > > XthenY(std::vector<std::vector<std::array<double, 4> > > u , double dx , double dy , double dt , int nxCells , int nyCells , double x0 ,double x1 , double y0 ,double y1 ,double tStart , double tStop , double C,double omega ){
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
                    uBarHalfL[i][j][k] = uBarL[i][j][k] - 0.5*(dt/dx)*(fluxX_def(uBarR[i][j])[k]-fluxX_def(uBarL[i][j])[k]);
                    uBarHalfR[i][j][k] = uBarR[i][j][k] - 0.5*(dt/dx)*(fluxX_def(uBarR[i][j])[k]-fluxX_def(uBarL[i][j])[k]);
                }
            }
        }
        applyBoundaryConditions(uBarHalfL, nxCells , nyCells);
        applyBoundaryConditions(uBarHalfR, nxCells , nyCells);
        //now we have flux in x direction
        for (int j = 1; j < nyCells + 3; j++) {
            for (int i = 1; i < nxCells + 3; i++) { 
                fluxX[i][j] = getFluxX(uBarHalfR[i][j], uBarHalfL[i+1][j], dx, dt);
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
                    uBarHalfL[i][j][k] = uBarL[i][j][k] - 0.5*(dt/dy)*(fluxY_def(uBarR[i][j])[k]-fluxY_def(uBarL[i][j])[k]);
                    uBarHalfR[i][j][k] = uBarR[i][j][k] - 0.5*(dt/dy)*(fluxY_def(uBarR[i][j])[k]-fluxY_def(uBarL[i][j])[k]);
                }
            }
        }
        applyBoundaryConditions(uBarHalfL, nxCells , nyCells);
        applyBoundaryConditions(uBarHalfR, nxCells , nyCells);
        //now we have flux in y direction
        for (int j = 1; j < nyCells + 3; j++) {
            for (int i = 1; i < nxCells + 3; i++) { 
                fluxY[i][j] = getFluxY(uBarHalfR[i][j], uBarHalfL[i][j+1], dy, dt);
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


std::vector<std::vector<std::array<double, 4> > > YthenX(std::vector<std::vector<std::array<double, 4> > > u , double dx , double dy , double dt , int nxCells , int nyCells , double x0 ,double x1 , double y0 ,double y1 ,double tStart , double tStop , double C,double omega ){
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
                    uBarHalfL[i][j][k] = uBarL[i][j][k] - 0.5*(dt/dy)*(fluxY_def(uBarR[i][j] )[k]-fluxY_def(uBarL[i][j] )[k]);
                    uBarHalfR[i][j][k] = uBarR[i][j][k] - 0.5*(dt/dy)*(fluxY_def(uBarR[i][j] )[k]-fluxY_def(uBarL[i][j] )[k]);
                }
            }
        }
        applyBoundaryConditions(uBarHalfL, nxCells , nyCells);
        applyBoundaryConditions(uBarHalfR, nxCells , nyCells);
        //now we have flux in y direction
        for (int j = 1; j < nyCells + 3; j++) {
            for (int i = 1; i < nxCells + 3; i++) { 
                fluxY[i][j] = getFluxY(uBarHalfR[i][j], uBarHalfL[i][j+1], dy, dt);
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
                    uBarHalfL[i][j][k] = uBarL[i][j][k] - 0.5*(dt/dx)*(fluxX_def(uBarR[i][j] )[k]-fluxX_def(uBarL[i][j] )[k]);
                    uBarHalfR[i][j][k] = uBarR[i][j][k] - 0.5*(dt/dx)*(fluxX_def(uBarR[i][j] )[k]-fluxX_def(uBarL[i][j] )[k]);
                }
            }
        }
        applyBoundaryConditions(uBarHalfL, nxCells , nyCells);
        applyBoundaryConditions(uBarHalfR, nxCells , nyCells);
        //now we have flux in y direction
        for (int j = 1; j < nyCells + 3; j++) {
            for (int i = 1; i < nxCells + 3; i++) { 
                fluxY[i][j] = getFluxX(uBarHalfR[i][j], uBarHalfL[i+1][j], dx, dt);
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
    int nxCells = 100;
    int nyCells = 100;
    double x0 = 0.0;
    double x1 = 1.0;
    double y0 = 0.0;
    double y1 = 1.0;
    double tStart = 0.0;
    double tStop = 0.012/std::pow(10,2.5);
    double C = 0.8;
    double omega =0;

    std::vector<std::vector<std::array<double, 4> > > u;
    u.resize(nxCells+4, std::vector<std::array<double, 4> >(nyCells + 4)); //set up u
    std::vector<std::vector<std::array<double, 4> > > v;
    v.resize(nxCells+4, std::vector<std::array<double, 4> >(nyCells + 4)); //set up v to check prim->conserv

    double dx = (x1 - x0) / nxCells;
    double dy = (y1 - y0) / nyCells;

    //intial conditions!
    for(int i = 0; i < u.size(); i++) { 
        for(int j = 0; j < u[0].size(); j++) {
            array prim;
            // coordinates
            double x = x0 + (i - 1.5) * dx;
            double y = y0 + (j - 1.5) * dy;

            
            if ( x>0.5 && y>0.5) {
                prim[0] = 1.0;   // density (rho)
                prim[1] = 0;   // x-velocity (u)
                prim[2] = 0;   // y-velocity (v)
                prim[3] = 0.01*1e5;   // pressure (p)
            } 
            else if(x<=0.5 && y > 0.5) {
                prim[0] = 1.0;
                prim[1] = 0;
                prim[2] = 0;
                prim[3] = 1000*1e5;
            }
            else if(x<=0.5 && y<=0.5) {
                prim[0] = 1.0;
                prim[1] = 0;
                prim[2] = 0;
                prim[3] = 1000*1e5;
            }
            else{
                prim[0] = 1.0;
                prim[1] = 0;
                prim[2] = 0;
                prim[3] = 0.01*1e5;
            }

            std::cout<<"i= "<<i<<" j= "<<j<<std::endl;
            u[i][j] = PrimativeToConservative(prim);
            v[i][j] = ConservativeToPrimative(u[i][j]);
            std::cout<<u[i][j][0]<<" "<<u[i][j][1]<<" "<<u[i][j][2]<<" "<<u[i][j][3]<<std::endl;
            std::cout<<v[i][j][0]<<" "<<v[i][j][1]<<" "<<v[i][j][2]<<" "<<v[i][j][3]<<std::endl;
            std::cout<<prim[0]<<" "<<prim[1]<<" "<<prim[2]<<" "<<prim[3]<<std::endl;
        }
    }
    std::cout<<"finished loop!"<<std::endl;

    applyBoundaryConditions(u , nxCells , nyCells);

    double dt;
    double t = tStart;

    do{
        dt = ComputeTimeStep(u,C,dx,dy);
        t +=dt;

        std::cout << "t = "<< t<<" dt = "<< dt<< std::endl; 
        applyBoundaryConditions(u , nxCells , nyCells);

        std::vector<std::vector<std::array<double, 4> > > uPlus1x = XthenY(u ,  dx ,  dy ,  dt ,  nxCells ,  nyCells ,  x0 , x1 ,  y0 , y1 , tStart ,  tStop ,  C, omega);
        std::vector<std::vector<std::array<double, 4> > > uPlus1y = YthenX(u ,  dx ,  dy ,  dt ,  nxCells ,  nyCells ,  x0 , x1 ,  y0 , y1 , tStart ,  tStop ,  C, omega);

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
            results[i][j] = ConservativeToPrimative(u[i][j]);
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
