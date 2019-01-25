//                                                  
//                                                        Write results 
//  Author: Max Braun
//  Date: 29/10/2016 
//
//
//----------------------------------------------------------------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;
//----------------------------------------------------------------------------------------------------------------------------
    char file_cship[50]      = "INPUT/c_ship.dat";
//----------------------------------------------------------------------------------------------------------------------------
    char file_interim[50]    = "interim.dat";
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
//
//                                                     - Load INPUT -  
//
//----------------------------------------------------------------------------------------------------------------------------
double load_interim(int i, int  variable , int lines )
{
std::ifstream infile(file_interim);
double x[column];
double numbers[lines];
int indx=1;
while(infile >> x[0] >> x[1] >> x[2] >> x[3]) 
{
numbers[indx-1]= x[variable];
indx = indx + 1;
}
double cdata = numbers[i];
return cdata;
}
//-------------------------------------------------------------------------------------------------------------------------------
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
//                                                main
//
//----------------------------------------------------------------------------------------------------------------------------
int main(){
//----------------------------------------------------------------------------------------------------------------------------
double M0 = load_cship(0);
double BC = load_cship(1);
int lines = line_counter(file_interim);
double Vign = load_interim(0,1, lines);
cout << BC << '\t' << M0 << '\t' << Vign << endl;
}
