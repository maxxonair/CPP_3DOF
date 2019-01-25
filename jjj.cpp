//----------------------------------------------------------------------------------------------------------------------------
//
//                                                      TrajectSim V1.03
//
//              3 degree of freedom flight dynamic simulation for atmospheric entry descent and landing trajectories
// 
//  Author: Max Braun
//  Date:  27/10/2016
//
//                             - Requires the free C++ boost library (https://www.boost.org)! 
//                             - Input files: 
//                                              * atm_*** -> atmosphere data (temperature, pressure/density, gas constant,
//                                                           heat capacity ratio) - 1D / altitude depending
//                                              
//                                              * c_***   -> control files ***_initial -> initial values for postion and 
//                                                                                        velocity vector
//                                                                         ***_integ   -> for numberical integration(t0,t1,
//                                                                                        step)
//                                                                         ***_planet  -> target planet data
//                                                                         ***_ship    -> entry vehicle data
//                                                                         ***_event   -> EDL event control file
//
//                                              * i_***   -> input files   ***_Cd    -> Drag coefficient data (continuum flow)
//                                                                                      Mach/AoA depending
//                                                                         ***_bank  -> Bank angle -> time depending 
//                                                                         ***_thrust-> Thrust profile, time depending/ t0 at
//                                                                                      ignition time!
//
double DIV = 1.0;
double tconst = 1000;
//----------------------------------------------------------------------------------------------------------------------------
#include <iostream>
#include <boost/array.hpp>
#include <fstream>
#include <cmath>
#include <boost/numeric/odeint.hpp>
# define lines(array)  ((sizeof(array))/(sizeof(*array)))
using namespace std;
using namespace boost::numeric::odeint;
//----------------------------------------------------------------------------------------------------------------------------
    long D ;                               // Drag Froce                               [N]
    long L ;                               // Lift Force                               [N]
    long Ty;                               // Side Force                               [N]
    double bank = 0;                       // Bank angle                               [rad]
    const double PI = 3.14159265359;       // PI                                       [-]
    const double kB = 1.380650424e-23;     // Boltzmann constant                       [SI]
    const double sigma = 1.6311e-9;        // Average collision diameter (CO2)         [m]
//---------------------------------------------------------------
//
//                    INPUT - Files
//
//--------------------------------------------------------------
    char file_rho[50]     = "INPUT/atm_rho";
    char file_T[50]       = "INPUT/atm_T";
    char file_R[50]       = "INPUT/atm_R";
    char file_gamma[50]   = "INPUT/atm_gamma";
    char file_cinteg[50]  = "INPUT/c_integ.dat";
    char file_cinitial[50]= "INPUT/c_initial.dat";
    char file_cplanet[50] = "INPUT/c_planet.dat";
    char file_cship[50]   = "INPUT/c_ship.dat";
    char file_cevent[50]  = "INPUT/c_event.dat";
    char file_thrust[50]  = "INPUT/i_thrust.dat";
    char file_cd[50]      = "INPUT/i_Cd.dat";
    char file_bank[50]    = "INPUT/i_bank.dat";
    char file_mheatshield[50] = "INPUT/interim_mheatshield.dat"; 
//----------------------------------------------------------------------------------------------------------------------------
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
// cout << line_count << endl;
return line_count;
}
//----------------------------------------------------------------------------------------------------------------------------
double LinearInterpolate( double atm_x[] , double atm_y[] , double xx, int lines)
{
int ii;
double yvalue;
double y1,y2,x1,x2;
//............................................
for(ii=lines;ii>0;ii--)
{
if(atm_x[ii]>xx){
y1 = atm_y[ii];
x1 = atm_x[ii];
}
}
//............................................
for(ii=0;ii<lines;ii++)
{
if(atm_x[ii]<xx){
y2 = atm_y[ii];
x2 = atm_x[ii];
}
}
//...........................................
if(xx<atm_x[0]){
yvalue = atm_y[0];
}
else{
yvalue = y1 + ( y2 - y1 ) * ( xx - x1 ) / ( x2 - x1 ) ; 
}
return yvalue;
}
//----------------------------------------------------------------------------------------------------------------------------
//
//                                    Load INPUT - data from control files ( c_*** ) 
//
//----------------------------------------------------------------------------------------------------------------------------
double mheatshield = 0;
int mstop =1; 
double load_mheatshield()
{
std::ifstream infile(file_mheatshield);
double a;
double numbers[1];
char b[4];
int indx=1;
while(infile >> a ) 
{
numbers[indx-1]= a;
indx = indx + 1;
}
double cdata = numbers[0];
return cdata;
}
//-----------------------------------------------------------------------------------------------------------------------------
double load_cinteg(int i)
{
std::ifstream infile(file_cinteg);
double a;
char b[4];
double numbers[10];
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
double load_cinitial( int i )
{
std::ifstream infile(file_cinitial);
double a;
char b[12];
double numbers[6];
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
double load_cplanet( int i )
{
std::ifstream infile(file_cplanet);
double a;
char b[12];
double numbers[3];
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
double m = load_cship(0);
double Lt = load_cship(2);
double Spar = (load_cship(3)/2) *(load_cship(3)/2) * PI; 
double LD = load_cship(4);
int    AoA = load_cship(5);
int    Bank_on_off = load_cship(6);
double ISP = load_cship(7);
double TTWin = load_cship(8);
const double BC_set = load_cship(1); 
const double S = m /( 1.5 *  BC_set ) ;             // Entry vehicles surface area [m²]
//---------------------------------------------------------------------------------------------------------------------------
double load_cevent(int i)
{
std::ifstream infile(file_cevent);
double a;
char b[14];
double numbers[10];
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
//............................................                                       .........................................
//
//                                     Load local atmosphere properties from atm_*** files 
//
//----------------------------------------------------------------------------------------------------------------------------
const int lines_rho    = line_counter(file_rho);
const int lines_T      = line_counter(file_T);
const int lines_R      = line_counter(file_R);
const int lines_gamma  = line_counter(file_gamma);
//------------------------ ------------ ------------  -----------   -----------------    ------------     ------------
double get_density(double h)                            // Density [kg/m³]
{
std::ifstream infile(file_rho);
double a,b;
double atm_rho[lines_rho][2];
int indx=1;
while (infile >> a >> b)
{
atm_rho[indx-1][0]=a;
atm_rho[indx-1][1]=b;
indx=indx+1;
}
indx=indx-1;
int i;
double density;
double atm_x[lines_rho] ;
double atm_y[lines_rho] ;
for(i=0;i<lines_rho;i++){
atm_x[i] = atm_rho[i][0];
atm_y[i] = atm_rho[i][1];
}
density = LinearInterpolate(atm_x,atm_y,h,lines_rho);
return density ; 
}
//---------------------------------------------------------------------------------------------------------------------------
double get_temperature(double h)                     // Temperature [K]
{
std::ifstream infile(file_T);
double a,b;
double atm_T[lines_T][2];
int indx=1;
while (infile >> a >> b)
{
atm_T[indx-1][0]=a;
atm_T[indx-1][1]=b;
indx=indx+1;
}
indx=indx-1;
int i;
double temperature;
double atm_x[lines_T] ;
double atm_y[lines_T] ;
for(i=0;i<lines_T;i++){
atm_x[i] = atm_T[i][0];
atm_y[i] = atm_T[i][1];
}
temperature = LinearInterpolate(atm_x,atm_y,h, lines_T);
return temperature ; 
}
//---------------------------------------------------------------------------------------------------------------------------
double get_R(double h)                     // Gas constant
{
std::ifstream infile(file_R);
double a,b;
double atm_R[lines_R][2];
int indx=1;
while (infile >> a >> b)
{
atm_R[indx-1][0]=a;
atm_R[indx-1][1]=b;
indx=indx+1;
}
indx=indx-1;
int i;
double R;
double atm_x[lines_R] ;
double atm_y[lines_R] ;
for(i=0;i<lines_R;i++){
atm_x[i] = atm_R[i][0];
atm_y[i] = atm_R[i][1];
}
R = LinearInterpolate(atm_x,atm_y,h, lines_R);
return R ; 
}
//---------------------------------------------------------------------------------------------------------------------------
double get_gamma(double h)                      // Heat capacity ratio
{
std::ifstream infile(file_gamma);
double a,b;
double atm_gamma[lines_gamma][2];
int indx=1;
while (infile >> a >> b)
{
atm_gamma[indx-1][0]=a;
atm_gamma[indx-1][1]=b;
indx=indx+1;
}
indx=indx-1;
int i;
double atm_x[lines_gamma] ;
double atm_y[lines_gamma] ;
for(i=0;i<lines_gamma;i++){
atm_x[i] = atm_gamma[i][0];
atm_y[i] = atm_gamma[i][1];
}
double gamma = LinearInterpolate(atm_x,atm_y,h, lines_gamma);
return gamma ; 
}
//----------------------------------------------------------------------------------------------------------------------------
//
//                                         Calculate local gravity vector ( 2D ) 
//
//----------------------------------------------------------------------------------------------------------------------------
const long rm = load_cplanet(0);                             // Planets average radius [m]
double get_gr(double r, double lat)                          // Gravity acceleration in radial direction [m/s²]
{
long long  mu = 4.28384E13;   // Standard gravitational constant (Mars) in SI
double phi = 3.1415/2-lat;    // Phi [rad]
double J2 = 1.9554537E-3;     // First Jeofreys? Constant
double J3 = 0;                // Second ...
double J4 = 0;                // Third ... 
double gr;                    // Gravitational acceleration in radial -> _r and north-south -> _n direction
gr = mu * ( 1 - 1.5 * J2 * ( 3 * cos(phi) * cos(phi) - 1) * (rm/r) * (rm/r) - 2 * J3 * cos(phi) * (5 * cos(phi) * cos(phi) - 3) * (rm/r) * (rm/r) * (rm/r) - (5/8) * J4 * (35 * cos(phi) * cos(phi) * cos(phi) * cos(phi) - 30 * cos(phi) * cos(phi) + 3) * (rm/r) * (rm/r) * (rm/r) * (rm/r) )/(r * r); 
return gr;
}
// ......................................................
double get_gn(double r, double lat)        // Gravity in north-south direction [m/s²]
{
long long mu = 4.28384E13;    // Standard gravitational constant (Mars) in SI
double phi = 3.1415/2-lat;    // Phi [rad]
double J2 = 1.9554537E-3;     // First Jeofreys? Constant
double J3 = 0;                // Second ...
double J4 = 0;                // Third ... 
double gn;                    // Gravitational acceleration in radial -> _r and north-south -> _n direction
gn = -3 * mu * sin(phi) * cos(phi) * (rm/r) * (rm/r) * (J2 + 0.5 * J3 * ( 5*cos(phi) *cos(phi) -1) * (rm/r)/cos(phi) + (5/6) * J4 * ( 7 * cos(phi) * cos(phi) - 1) * (rm/r) * (rm/r) ) /(r * r);
//......................................................
return gn;
}
//----------------------------------------------------------------------------------------------------------------------------
//
//                                          Calculate Drag with three flowzone approach 
//                                        Free molecular -> transitional -> Contiuum flow
// 
//----------------------------------------------------------------------------------------------------------------------------
double load_Drag(double vel, double h, double P, double T, double CdC)
{
double CD;
double Kn = kB * T / ( sqrt(2) * PI * sigma * sigma * P * Lt );
if(Kn<0.1){
//                 Continuum flow        <---------------
CD=CdC;
}
if(Kn>0.1 && Kn<10){
//                 Transtional zone      <---------------
double S = vel / sqrt(2 * get_R(h) * T);
double Cdfm= 1.75 + sqrt(PI)/(2 * S);
CD= CdC + ( Cdfm - CdC ) * ( 1/3 * log10( Kn / sin( PI / 6 ) ) * 0.5113 ) ;
}
if(Kn>10){
//                 Free molecular zone   <---------------
double S = vel / sqrt(2 * get_R(h) * T);
CD= 1.75 + sqrt(PI)/(2 * S);
}
return CD;
}
//----------------------------------------------------------------------------------------------------------------------------
int calc_flowzone( double vel, double h , double P , double T)
{
double Lt = load_cship(2);
double Kn = kB * T / ( sqrt(2) * PI * sigma * sigma * P * Lt );
int flowzone;
if(Kn<0.1){
flowzone = 1;
}
if(Kn>0.1 && Kn<10){
flowzone = 2;
}
if(Kn>10){
flowzone = 3;
} 
return flowzone; 
}
//----------------------------------------------------------------------------------------------------------------------------
//
//                                                EDL - event management 
//
//----------------------------------------------------------------------------------------------------------------------------
double load_Cdpar(double vel, double qinf, double Ma, double h )
{
double Cdpar = 0 ;
double var_event   = load_cevent(0);
double var_trigger = load_cevent(1);
double val_trig1   = load_cevent(3);
double x_trig1, x_trig2;
double CDPAR_value = 0.58;                   // Average value
//...........................................
if(var_trigger == 1){
x_trig1 = vel;
}
if(var_trigger == 2){
x_trig1 = qinf;
}
if(var_trigger == 3){
x_trig1 = Ma;
}
if(var_trigger == 4){
x_trig1 = h; 
}
//...........................................
if(var_event == 4){
double var_trigger2 = load_cevent(2);
double val_trig2 = load_cevent(4);
if(var_trigger2 == 1 ){
x_trig2 = vel;
}
if(var_trigger2 == 2){
x_trig2 = qinf ;
}
if(var_trigger2 == 3){
x_trig2 = Ma;
}
if(var_trigger2 == 4){
x_trig2 = h;
}
if(x_trig1 < val_trig1 && x_trig2 > val_trig2){
Cdpar = CDPAR_value;
}
}
//...........................................
if( var_event == 1 ){
Cdpar = 0;
}
if( var_event == 2 ){
if(x_trig1 < val_trig1 ){
Cdpar = CDPAR_value;
}
}
if(var_event == 3 ){
Cdpar = 0;
}
return Cdpar;
}
double mm = load_cship(0);
//----------------------------------------------------------------------------------------------------------------------------
double get_Thrust(double m, double phi, double v)
{
double thrust;
double v_vert = v * ( - sin(phi));
double TTWout = ( 1/5 ) * TTWin;
thrust = TTWin * mm - TTWout * m  ;
return thrust;
}
//.........................................................
double load_Thrust(double vel, double qinf, double Ma, double h , double t, double m, double phi, double v)
{
double Thrust=0;
double x_trig1,x_trig2;
double var_event = load_cevent(0);
double var_trigger = load_cevent(1);
double val_trig1   = load_cevent(3);
if(var_trigger == 1){
x_trig1 = vel;
}
if(var_trigger == 2){
x_trig1 = qinf;
}
if(var_trigger == 3){
x_trig1 = Ma;
}
if(var_trigger == 4){
x_trig1 = h; 
}
//.................................
if(var_event == 1 ){
Thrust = 0;
}
if(var_event == 2 ){
Thrust = 0;
}
if(var_event == 3 ){
if(x_trig1 < val_trig1 ){
Thrust = get_Thrust(m, phi, v); 
}
}
if(var_event == 4){
double var_trigger2 = load_cevent(2);
double val_trig2 = load_cevent(4);
if(var_trigger2 == 1 ){
x_trig2 = vel;
}
if(var_trigger2 == 2){
x_trig2 = qinf ;
}
if(var_trigger2 == 3){
x_trig2 = Ma;
}
if(var_trigger2 == 4){
x_trig2 = h;
}
if( x_trig2 < val_trig2){
Thrust = get_Thrust(m, phi, v);
}
}
return Thrust;
}
//----------------------------------------------------------------------------------------------------------------------------
double Th_t0 = 0;
double Th_t1 = 0;
double find_Th_t0t(double vel, double qinf, double Ma, double h , double time)
{
double x_trig1,x_trig2;
int Th_t0 = 0;
double var_event = load_cevent(0);
double var_trigger = load_cevent(1);
double val_trig1   = load_cevent(3);
if(var_trigger == 1){
x_trig1 = vel;
}
if(var_trigger == 2){
x_trig1 = qinf;
}
if(var_trigger == 3){
x_trig1 = Ma;
}
if(var_trigger == 4){
x_trig1 = h; 
}
//.................................
if(var_event == 1 ){
}
if(var_event == 2 ){
}
if(var_event == 3 ){
if(x_trig1 < val_trig1 ){
Th_t0 = time;
}
}
if(var_event == 4){
double var_trigger2 = load_cevent(2);
double val_trig2 = load_cevent(4);
if(var_trigger2 == 1 ){
x_trig2 = vel;
}
if(var_trigger2 == 2){
x_trig2 = qinf ;
}
if(var_trigger2 == 3){
x_trig2 = Ma;
}
if(var_trigger2 == 4){
x_trig2 = h;
}
if( x_trig2 < val_trig2){
Th_t0 = time; 
}
}
return Th_t0;
}
//----------------------------------------------------------------------------------------------------------------------------//
//                        Contituum flow Drag coefficient - depending on Angle of Attack and Mach number
//                                              Bank angle - time depending 
//
//----------------------------------------------------------------------------------------------------------------------------
double get_CdC(double Ma, int AoA)
{
std::ifstream infile(file_cd);
double a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r;
int lin = line_counter(file_cd);
double i_Cd[lin][18];
if(AoA>18){
AoA = 16;
}
if(AoA < 0){
AoA = 0 ;
}
int indx=1;
while (infile >> a >> b >> c >> d >> e >> f >> g >> h >> i >> j >> k >> l >> m >> n >> o >> p >> q >> r )
{
i_Cd[indx-1][0]=a;
i_Cd[indx-1][1]=b;
i_Cd[indx-1][2]=c;
i_Cd[indx-1][3]=d;
i_Cd[indx-1][4]=e;
i_Cd[indx-1][5]=f;
i_Cd[indx-1][6]=g;
i_Cd[indx-1][7]=h;
i_Cd[indx-1][8]=i;
i_Cd[indx-1][9]=j;
i_Cd[indx-1][10]=k;
i_Cd[indx-1][11]=l;
i_Cd[indx-1][12]=m;
i_Cd[indx-1][13]=n;
i_Cd[indx-1][14]=o;
i_Cd[indx-1][15]=p;
i_Cd[indx-1][16]=q;
i_Cd[indx-1][17]=r;
indx=indx+1;
}
indx=indx-1;
int ii;
double dat_x[lin] ;
double dat_y[lin] ;
for(ii=0;ii<lin;ii++){
dat_x[ii] = i_Cd[ii][0];
dat_y[ii] = i_Cd[ii][AoA + 1];
}
double CdC = LinearInterpolate(dat_x,dat_y,Ma, lin);
return CdC ; 
}
//----------------------------------------------------------------------------------------------------------------------------
double get_bank(double time, int Bank_on_off)
{
std::ifstream infile(file_bank);
double a,b;
int lin = line_counter(file_bank);
if(Bank_on_off == 1){
double i_bank[lin][2];
int indx=1;
while (infile >> a >> b )
{
i_bank[indx-1][0]=a;
i_bank[indx-1][1]=b;
indx=indx+1;
}
indx=indx-1;
int ii;
double dat_x[lin] ;
double dat_y[lin] ;
for(ii=0;ii<lin;ii++){
dat_x[ii] = i_bank[ii][0];
dat_y[ii] = i_bank[ii][1];
}
bank = LinearInterpolate(dat_x,dat_y,time, lin) * PI/180;
}
else{
bank = 0;
}
return bank ; 
}
//----------------------------------------------------------------------------------------------------------------------------
typedef boost::array< double, 6 > state_type;
//----------------------------------------------------------------------------------------------------------------------------
//
//                                           -->    Equation of motion     <--
//
//----------------------------------------------------------------------------------------------------------------------------
const double omega = load_cplanet(1);                      // Planets rotational speed [rad/s]
double gr, gn, rho, qinf, T, R, Ma, CdPar, Thrust, CdC, Cd, Cl;
double vminone = load_cinitial(3); 
double XminOne = 0;
double X_horizontal = 0;
//--------------------------------------
void eqom( const state_type &x, state_type &dxdt , double t )
{
//............................................................................................................................
 gr = get_gr(x[2], x[1]);                          		  // gravity acceleration in radial direction [m/s²]
 gn = get_gn(x[2], x[1]);                            		  // gravity acceleration in north-south direction [m/s²]
//............................................................................................................................
 rho = get_density( x[2] - rm ) ;                    		  // density [kg/m³]
//............................................................................................................................
  qinf = 0.5 * rho * ( x[3] * x[3] ) ;               		  // Dynamic pressure [Pa]
  T = get_temperature(x[2] - rm) ;                   		  // Temperature [K]
double   gamma = get_gamma( x[2] - rm) ;              		  // 
  R = get_R( x[2] - rm ) ;                           		  // Gas Constant [Si]
  Ma = x[3] / sqrt( T * gamma * R);                  		  // Mach number [-]
  CdPar  = load_Cdpar( x[3], qinf, Ma, x[2] - rm);   		  // Parachute Drag coefficient [-]
  Thrust = load_Thrust(x[3], qinf, Ma, x[2] - rm, Th_t1, m , x[4], x[3]);        // Thrust [N]
  CdC = DIV  * get_CdC(Ma, AoA);                                         // Continuum flow drag coefficient [-]
  Cd = load_Drag( x[3], x[2] - rm , get_density(x[2] -rm) * get_R(x[2] -rm ) * get_temperature(x[2] -rm) , get_temperature(x[2] - rm), CdC ); 
if(Thrust == 0){
  Cl = LD * Cd; 
}
else
{ 
Cl = 0;
mheatshield = load_mheatshield();
}                                                  // Lift coefficient [-]
  bank = get_bank(t, Bank_on_off);                                // Bank angle [rad]
//-----------------------------------------------------------------------------------------------
  D  = - qinf * S * Cd - qinf * Spar * CdPar - Thrust;            // Drag Force [N]
  L  =   qinf * S * Cl * cos( bank ) ;                            // Lift Force [N]
  Ty =   qinf * S * Cl * sin( bank ) ;                            // Side Force [N]
//----------------------------------------------------------------------------------------------
if(Thrust == 0  ){

}
else
{
if(mstop ==1){
mheatshield = load_mheatshield();
mstop=2;
}
else
{
mheatshield=0;
}
} 
  m  =   m - Thrust/(ISP*9.81)*load_cinteg(2) - mheatshield;                    // EDL vehicle mass [kg]
//--------------------------------------------------------
if(Th_t0 == 0){
 Th_t0 =  find_Th_t0t(x[3], qinf,  Ma, x[2]-rm , t);
}
if(Th_t0 != 0 ){
 Th_t1 = t - Th_t0;
}
// .........................................................................................................................  
    dxdt[0] = x[3] * cos( x[4] ) * sin( x[5] ) / ( x[2] * cos( x[1] ) );
    dxdt[1] = x[3] * cos( x[4] ) * cos( x[5] ) / x[2] ;
    dxdt[2] = x[3] * sin( x[4] );
    dxdt[3] = -gr * sin( x[4] ) + gn * cos( x[5] ) * cos( x[4] ) + D / m + omega * omega * x[2] * cos( x[1] ) * ( sin( x[4] )               * cos( x[1] ) - cos( x[1] ) * cos( x[5] ) * sin( x[1] ) ) ;
    dxdt[4] = ( x[3] / x[2] - gr/ x[3] ) * cos( x[4] ) - gn * cos( x[5] ) * sin( x[4] ) / x[3] + L / ( x[3] * m ) + 2 * omega               * sin( x[5] ) * cos( x[1] ) + omega * omega * x[2] * cos( x[1] ) * ( cos( x[4] ) * cos( x[1] ) + sin( x[4] ) *                cos( x[5] ) * sin( x[1] ) ) / x[3] ;
    dxdt[5] = x[3] * sin( x[5] ) * tan( x[1] ) * cos( x[4] ) / x[2] - gn * sin( x[5] ) / x[3] - Ty / ( x[3] - cos( x[4] ) * m               ) - 2 * omega * ( tan( x[4] ) * cos( x[5] ) * cos( x[1] ) - sin( x[1] ) ) + omega * omega * x[2] * sin( x[5] ) *              sin( x[1] ) * cos( x[1] ) / ( x[3] * cos( x[4] ) ) ;
}
//----------------------------------------------------------------------------------------------------------------------------
//
//                                                   Write results
//
//----------------------------------------------------------------------------------------------------------------------------
double step = load_cinteg(2);
void write_res( const state_type &x , const double t )
{
double   g = sqrt( gr * gr + gn * gn );                       // Avarage gravity aceleration [m/s²]
double   P = rho * R * T ;                                    // Static ambient pressure [Pa]
int      flowzone = calc_flowzone( x[3], x[2] -rm , P, T);    // Flowzone variable 
double   gamma = get_gamma( x[2] - rm) ;                      // Heat capacity ratio [Si]
double   acc = ((x[3] - vminone) / step ) / 9.81;             // Normalized deceleration [-]
vminone = x[3];
double   BC = m/(Cd * S);                         // Ballistic Coefficient [kg/m²]
//................................
int event=0;
if(Thrust!=0 && CdPar ==0){
event=1;
}
if(Thrust == 0 && CdPar != 0){
event=2;
}
//...............................
// Horizontal track computation: 
double V_horizontal = x[3] * cos(x[4]);
X_horizontal = XminOne + V_horizontal * step;
XminOne = X_horizontal;
//...............................
double heating = 1.87e-4 * sqrt(rho/1.5) * x[3] * x[3] * x[3] ;
//.................................................................
    cout << t << '\t' << x[0] << '\t' << x[1] << '\t' << x[2]-rm << '\t' << x[3] << '\t' << x[4] << '\t' << x[5] <<  '\t' << rho << '\t' << D << '\t' << L << '\t' << Ty << '\t' << gr << '\t' << gn << '\t' << g << '\t' << T << '\t' <<Ma << '\t' << gamma << '\t' << R << '\t' << P << '\t' << Cd << '\t' << Cl << '\t' << heating << '\t' << flowzone << '\t' << qinf << '\t' << CdPar <<  '\t' << Thrust << '\t' << CdC << '\t' << bank << '\t' << acc << '\t' << m << '\t' << event << '\t' << BC << '\t' << X_horizontal << '\t' << S << endl;
}
//----------------------------------------------------------------------------------------------------------------------------
//
//                                                     Integrator
//
//----------------------------------------------------------------------------------------------------------------------------
int main(int argc, char **argv)
{
double t0   = load_cinteg(0);                 // Start time [s]
double t1   = load_cinteg(1);                 // End time   [s]
double step = load_cinteg(2);                 // Step time  [s]
int i;
double init[6];
for(i=0;i<6;i++)
{
init[i] = load_cinitial(i);
}
//----------------------------------------------------------------------------------------------------------------------------
    state_type x = {    init[0]   ,   init[1]      , (init[2]+rm) ,   init[3]      ,  init[4]    ,       init[5]     }; 
//                    Long [rad]  ;   Lat [rad]    ; Radius [m]   ; Velocity [m/s] ;   FPA [rad] ; Local azimuth [rad] 
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//double stepper;
runge_kutta4< state_type > stepper;
    integrate_const(stepper,  eqom , x , t0 , t1 , step , write_res );  // Integrator <------------------ !!!
}
//----------------------------------------------------------------------------------------------------------------------------
