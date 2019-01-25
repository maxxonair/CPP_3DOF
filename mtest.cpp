//                                                  
//                                                        Mass LOOP 
//  Author: Max Braun
//  Date: 29/10/2016 
//
// 
//                         required Input files:
//                                                * interim2.dat (Vtd,htd,thrust, mfuel)
//
//                                                * c_mass_loop.dat : mass loop controll file
//                                                                    (1) Payload mass [kg]
//                                                                    (2) Phi_shield [kg/m²]
//                                                                    (3) Safety Margin [-]
//                                                                    (4) Fuel combination [1] LOX/LCH4 [2] MMH/NTO 
//                                                                    (5) Pressure tank [Pa]
//                                                                    (6) Average drag coefficient [-]
//
//                                                * c_ship.dat : required for BC 
// Tank material: titanium (fixed)
// Turbo pump engine
//
//                      -> Output:    M0 - Entry mass [kg] 
//
//----------------------------------------------------------------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;
//----------------------------------------------------------------------------------------------------------------------------
    char file_interim[50]    = "interim2.dat";
    char file_c_mass[50]     = "c_mass_loop.dat";
    char file_cship[50]      = "INPUT/c_ship.dat";
//----------------------------------------------------------------------------------------------------------------------------
const int column = 4; 
//----------------------------------------------------------------------------------------------------------------------------
int line_counter(char file[50])
{
std::ifstream infile(file);
char cc;
int line_count = 0;
while (infile.get(cc))
{
   if (cc == '\n')
        ++line_count;
}
return line_count;
}
//----------------------------------------------------------------------------------------------------------------------------
double load_cship(int i)
{
std::ifstream infile(file_cship);
double a;
char b[12];
double numbers[30];      // adapt number if you add variables ! 
int indx=1;
while(infile >> a >> b ) 
{
numbers[indx-1]= a;
indx = indx + 1;
}
double cdata = numbers[i];
return cdata;
}
//----------------------------------------------------------------------------------------------------------------------------
//
//                                                     - Load INPUT -  
//
//----------------------------------------------------------------------------------------------------------------------------
double load_from_file(int i , int lines, char file[50])
{
std::ifstream infile(file);
double x;
double numbers[lines];
int indx=1;
while(infile >> x ) 
{
numbers[indx-1]= x;
indx = indx + 1;
}
double cdata = numbers[i];
return cdata;
}
//----------------------------------------------------------------------------------------------------------------------------
// 
//                                                main
//
//----------------------------------------------------------------------------------------------------------------------------
int main(){
//                Line counter
int    lines_interim = line_counter(file_interim);
int    lines_c_mass  = line_counter(file_c_mass);
//--------------------------------------------------------------------------------------------------------------------------- 
// INPUT from interim2.dat
//---------------------------------------------------------------------------------------------------------------------------
double Mfuel = load_from_file(3, lines_interim, file_interim);
double T_TDS = load_from_file(2, lines_interim, file_interim);
//---------------------------------------------------------------------------------------------------------------------------
// INPUT from c_mass_loop.dat
//---------------------------------------------------------------------------------------------------------------------------
double M_payload = load_from_file(0, lines_c_mass, file_c_mass);
double Phi_shield = load_from_file(1, lines_c_mass, file_c_mass); 
double SafMar = load_from_file(2, lines_c_mass, file_c_mass);
double qmax = 2000;
int FuComb = load_from_file(3, lines_c_mass, file_c_mass);
FuComb = FuComb - 1; 
double P_tank_fuel = load_from_file(4, lines_c_mass, file_c_mass);
double Cd = load_from_file(5, lines_c_mass, file_c_mass);
double BC = load_cship(1);
// ---------------------------------------------------------
double g_0 = 9.81;
double MfacT = 5000;
double Dens_OX[2] = { 1141 ,  1442 } ;
double Dens_FU[2] = { 424 ,  880 } ;
double OF[2] = { 3.5 ,  2.0 } ; 
//----------------------------------------------------------
//                    Structure
//----------------------------------------------------------
double f_structure = (0.0232 * pow( qmax, 0.1708 ) ) ; //Structure mass fraction 
double f_backshell = 0.14 ; // Backshell mass fraction
//----------------------------------------------------------
//                   Reaction Control System
//----------------------------------------------------------
double f_eng_rcs = 0.005 ; // RCS engines
double f_prop_rcs = 0.0101 ; // RCS propellants
//----------------------------------------------------------
// Dependend mass: Sum factors:
double SumFdep = f_structure + f_backshell + f_eng_rcs + f_prop_rcs;
//----------------------------------------------------------------------------------------
// Independend mass: Thrusted Descend System
//----------------------------------------------------------------------------------------
double dens_ox = Dens_OX[FuComb];   				 // Density oxidizer 
double dens_fu = Dens_FU[FuComb];   				 // Density fuel
double of = OF[FuComb];              				 // Oxidizer to fuel ratio
double m_fu = Mfuel / (of + 1 ) ;   				 // fuel mass
double m_ox = of * m_fu;            				 // oxidizer mass
double V_fu = m_fu / dens_fu;       				 // fuel volume
double V_ox = m_ox / dens_ox;      				 // oxidizer volume
double m_eng_tds = ( 0.00144 * T_TDS ) + 49.6 ;  		 // TDS engine mass
double m_tank_fuel_tds = P_tank_fuel * V_fu / ( g_0 * MfacT ) ;  // fuel tank mass
double m_tank_ox_tds = P_tank_fuel * V_ox / ( g_0 * MfacT ) ;    // oxidizer tank mass
//------------------------------------------------------------------------------------------
double M_Indep = m_eng_tds + m_tank_fuel_tds + m_tank_ox_tds + m_fu + m_ox;
double M0 = (M_payload + M_Indep * SafMar ) / ( 1 - SumFdep * SafMar - Phi_shield / ( BC * Cd ) ) ; 
//----------------------------------------------------------------------------------------------------------------------------
double structure = ( f_structure + f_backshell ) * M0;
double rcs = ( f_eng_rcs + f_prop_rcs ) * M0;
double mfs = Phi_shield /(BC * Cd ) ;
cout << " M0          = " << M0 << endl;
cout << "----------------------------------" << endl;
cout << " M_payload   = " << M_payload << endl;
cout << " M_structure = " << structure*SafMar << endl;
cout << " M_RCS       = " << rcs*SafMar << endl ; 
cout << " M_TDS       = " << M_Indep*SafMar << endl;
cout << " M_Heatshield= " << mfs*M0  << endl;
cout << "----------------------------------" << endl;
cout << " p_payload   = " << M_payload/M0*100*SafMar << endl;
cout << " p_structure = " << structure/M0*100*SafMar << endl;
cout << " p_RCS       = " << rcs/M0*100*SafMar << endl ; 
cout << " p_TDS       = " << M_Indep/M0*100*SafMar << endl;
cout << " p_Heatshield= " << mfs*100  << endl;
double total = (M_payload/M0*100 + structure/M0*100 + rcs/M0*100 + M_Indep/M0*100)*SafMar + mfs*100 ;
cout << ".................................." << endl;
cout << " p_total     = " << total << endl;
}
