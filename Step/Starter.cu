#include <iostream>
#include <fstream>
#include <cmath>
#include "ausmPlus.h"

using namespace std;

double gammma;
double mu;
double k;
double R;	

int main()
{
	gammma=1.4;
	mu=1.789*pow(10,-5);
	k=0.0257;
	R=286.9;
	
	double initial[4];
	double temp=300;
	double speed=1000;
	double speed_sound=sqrt(gammma*R*temp);
	double mach=speed/speed_sound;
	double pressure=101325/pow((1+(gammma-1)/2*mach*mach),gammma/(gammma-1));
	double density = pressure/(R*temp);
	double Re=(pressure/R/temp)*speed*13.85/mu;
	double Prandle=mu*1004/k;
	//Rho
	initial[0]=1;
	//Rho*U
	initial[1]=-1;
	//Rho *V
	initial[2]=0;
	//Rho*E, E is the internal energy including the kinetic energy(i.e. total intenal energy)
	initial[3]=1;
	//Time steps and delta_t
	double timesteps=750000;
	double deltat=0.01	;

	cout<<"pressure="<<pressure<<endl;
	cout<<" density="<<density<<endl;
	cout<<"static energy="<<R/(gammma-1)*temp*pow(101325/pressure,1/gammma)*initial[0]<<endl;
	cout<<"mach = "<<mach<<" speed = "<<speed<<endl;
	cout<<"Reynolds = "<<Re<<" Prandl = "<<Prandle<<endl;

	ausmplus(initial,timesteps,deltat,1,1,Re,Prandle,mach);

	return 0;
}