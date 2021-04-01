#include <iostream>
#include <cstdlib>

//Include files from the library:

//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"

//For certain lipid models (liposome, flat membrane, etc...)
#include "../models/lipidModels.h"
using namespace std;

#define CUTOFF 2.0
#define RMIN 1.0

int main(int argc, char **argv)
{
	if(argc<10)
	{
		cout << "Usage: " << argv[0] << " name seed Umin_Anchor_Head fraction_solvent_removed nLipids arialDensity density nMonomers nTess" << endl;
		return 0;
	}

	char *name=argv[1];
	double removeSolventFraction=atof(argv[4]);
	int nLipids=atoi(argv[5]);
	double arialDensity=atof(argv[6]);
	double density=atof(argv[7]);
	int nMonomers=atoi(argv[8]);
	int nTess=atoi(argv[9]);
	
	
	///the variables for the simulation
	
	Blob<double> System;
	
	System.setGamma(1.0);
	System.setNTypes(6);
	System.setSeed(atoi(argv[2]));
	
	threeVector<bool> wrap;//periodic boundaries flag, be forwarned, this doesn't work yet, it is always on
	wrap.x=true;
	wrap.y=true;
	wrap.z=true;
	System.setPeriodic(wrap);
	System.setCutoff(2.0);
	
	System.setInitialTime(0);
	System.setFinalTime(50000);
	System.setDeltaT(0.02);
	System.setStoreInterval(100);
	System.setMeasureInterval(10);
	//System.setDeltaLXY(0.01);//if you don't set it, it doesn't activate
	System.setInitialTemp(3.0);
	System.setFinalTemp(3.0);
	
	///initialize constants (force and potential constants)
	//generate force constants
	double *Umin=new double[System.readNTypes()*System.readNTypes()];
	double *Umax=new double[System.readNTypes()*System.readNTypes()];
	double *lbond=new double[System.readNTypes()*System.readNTypes()];
	double *kbond=new double[System.readNTypes()*System.readNTypes()*System.readNTypes()];
	double *abend=new double[System.readNTypes()*System.readNTypes()*System.readNTypes()];
	double *kbend=new double[System.readNTypes()*System.readNTypes()*System.readNTypes()];
	
	//could be one loop, but it's good for an example
	//of how to set constants
	for(int i=0;i<System.readNTypes();i++)
	{
		for(int j=0;j<System.readNTypes();j++)
		{
			//i=first type, j=second type, indexed grid
			int k=i+j*System.readNTypes();
			Umin[k]=0;
			Umax[k]=100;
			lbond[k]=0.7;
			kbond[k]=100;
		}
	}
	
	//umin and umax exceptions, tail types
	Umin[TAIL+TAIL*System.readNTypes()]=-6;
	Umax[TAIL+TAIL*System.readNTypes()]=200;
	
	//for polymer strand interactions
	Umin[MONOMER+MONOMER*System.readNTypes()]=0;
	Umax[MONOMER+MONOMER*System.readNTypes()]=100;
	
	Umin[HEAD+MONOMER*System.readNTypes()]=0;
	Umax[HEAD+MONOMER*System.readNTypes()]=100;
	
	Umin[MONOMER+HEAD*System.readNTypes()]=Umin[HEAD+MONOMER*System.readNTypes()];
	Umax[MONOMER+HEAD*System.readNTypes()]=Umax[HEAD+MONOMER*System.readNTypes()];

	Umin[TAIL+MONOMER*System.readNTypes()]=0;
	Umax[TAIL+MONOMER*System.readNTypes()]=100;

	Umin[MONOMER+TAIL*System.readNTypes()]=0;
	Umax[MONOMER+TAIL*System.readNTypes()]=100;
	
	Umin[HEAD+ANCHOR*System.readNTypes()]=atof(argv[3]);
	Umax[HEAD+ANCHOR*System.readNTypes()]=100;
	
	Umin[ANCHOR+HEAD*System.readNTypes()]=Umin[HEAD+ANCHOR*System.readNTypes()];
	Umax[ANCHOR+HEAD*System.readNTypes()]=Umax[HEAD+ANCHOR*System.readNTypes()];
	
//	Umin[HEAD+CYTO*System.readNTypes()]=-6;
//	Umax[HEAD+CYTO*System.readNTypes()]=200;
	
//	Umin[CYTO+HEAD*System.readNTypes()]=Umin[HEAD+CYTO*System.readNTypes()];
//	Umax[CYTO+HEAD*System.readNTypes()]=Umax[HEAD+CYTO*System.readNTypes()];
	
	//how to do it with one loop
	for(int i=0;i<System.readNTypes()*System.readNTypes()*System.readNTypes();i++)
	{
		kbend[i]=100;
		abend[i]=1;
	}
	
	//two body constants
	for(int i=0;i<System.readNTypes();i++)
	{
		for(int j=0;j<System.readNTypes();j++)
		{
			int g=i+j*System.readNTypes();
			//U[r<=rmin]=((Umax-Umin)*(rmin-r)^2/rmin^2)+Umin
			//U[rmin<r<=rc]=(-2*Umin*(rc-r)^3/(rc-rmin)^3)+(3*Umin*(rc-r)^2/(rc-rm)^2)
			//F[r<=rmin]=(-2*(Umax-Umin)/rmin^2)*(rmin-r)/r
			//F[rmin<r<=rc]=(6*Umin/(rc-rmin)^3)*(rc-r)^2/r-(6*Umin/(rc-rmin)^2)*(rc-r)/r
			
			//constants[U[rmin<r<=rc]]: 0:-2*Umin/(rc-rmin)^3 1:3*Umin/(rc-rmin)^2
			//constants[F[rmin<r<=rc]]: 2:-6*Umin/(rc-rmin)^3  3:6*Umin/(rc-rmin)^2
			
			//constants[U[r<=rmin]]: 4:(Umax-Umin)/rmin^2  5:Umin
			//constants[F[r<=rmin]]: 6:2*(Umax-Umin)/rmin^2
			//constants[general]: 7:rc 8:rmin 9:rc^2 10:rmin^2
			
			//F constants, force constants
			System.addTwoBodyFconst(RMIN);//C8
			System.addTwoBodyFconst((2.0*(Umax[g]-Umin[g]))/(RMIN*RMIN));//C6
			System.addTwoBodyFconst(0);//part of index trick
			System.addTwoBodyFconst(CUTOFF);//C7
			System.addTwoBodyFconst((6.0*Umin[g])/((CUTOFF-RMIN)*(CUTOFF-RMIN)));//C3
			System.addTwoBodyFconst((6.0*Umin[g])/((CUTOFF-RMIN)*(CUTOFF-RMIN)*(CUTOFF-RMIN)));//C2
			
			//U constants, potential constants
			System.addTwoBodyUconst(RMIN);//C8
			System.addTwoBodyUconst((Umax[g]-Umin[g])/(RMIN*RMIN));//C4
			System.addTwoBodyUconst(Umin[g]);//C5,no index trick
			System.addTwoBodyUconst(CUTOFF);//C7
			System.addTwoBodyUconst((3.0*Umin[g])/((CUTOFF-RMIN)*(CUTOFF-RMIN)));//C1
			System.addTwoBodyUconst((2.0*Umin[g])/((CUTOFF-RMIN)*(CUTOFF-RMIN)*(CUTOFF-RMIN)));//C0
		}
	}
	
	//just place it far away
	threeVector<double> pos;
	pos.x=100;
	pos.y=100;
	pos.z=100;
	
	double lipidConstants[4];
	lipidConstants[0]=0.7;
	lipidConstants[1]=100;
	lipidConstants[2]=1.0;
	lipidConstants[3]=100;
	
	double cytoBondLength=1.0;
	double cytoConstants[4];
	cytoConstants[0]=cytoBondLength;
	cytoConstants[1]=100;
	cytoConstants[2]=1.0;
	cytoConstants[3]=100;
	
	//Add liposome
	double radius=liposome(System, nLipids, 3, pos, 0.7, arialDensity, lipidConstants, 4);
	
	//Add cytoskeleton
	sphericalCyto(System, nMonomers, pos, radius, nTess, cytoConstants, 4);
	/*
	//Adjust size to accomidate solvent
	position<double> *p=System.getPositions();
	threeVector<double> lower,upper;
	lower=System.readSize();
	upper=0;
	for(int i=0;i<System.readNParticles();i++)
	{
		lower.x=(lower.x>p[i].x)?p[i].x:lower.x;
		lower.y=(lower.y>p[i].y)?p[i].y:lower.y;
		lower.z=(lower.z>p[i].z)?p[i].z:lower.z;
		
		upper.x=(upper.x<p[i].x)?p[i].x:upper.x;
		upper.y=(upper.y<p[i].y)?p[i].y:upper.y;
		upper.z=(upper.z<p[i].z)?p[i].z:upper.z;
	}
	
	for(int i=0;i<System.readNParticles();i++)
	{
		p[i].x+=((upper.x-lower.x)*1.0/8.0-lower.x);
		p[i].y+=((upper.y-lower.y)*1.0/8.0-lower.y);
		p[i].z+=((upper.z-lower.z)*1.0/8.0-lower.z);
		if(p[i].x<0 || p[i].x>(upper.x-lower.x)*(10.0/8.0) ||
			p[i].y<0 || p[i].y>(upper.y-lower.y)*(10.0/8.0) ||
			p[i].z<0 || p[i].z>(upper.z-lower.z)*(10.0/8.0))
		{
			std::cout << "Particle out of box!\n";
			std::cout << i << '\n';
			std::cout << lower.x << ' ' << lower.y << ' ' << lower.z << '\n';
			std::cout << upper.x << ' ' << upper.y << ' ' << upper.z << '\n';
			std::cout << p[i].x << ' ' << p[i].y << ' ' << p[i].z << '\n';
			return 0;
		}
	}
	
	upper.x=(upper.x-lower.x)*(10.0/8.0);
	upper.y=(upper.y-lower.y)*(10.0/8.0);
	upper.z=(upper.z-lower.z)*(10.0/8.0);
	
	System.setSize(upper);
	*/
	//Add solvent
	//solventFill<double> (System,density);
	
	//remove solvent fraction
	//System.setRemoveSolvent(removeSolventFraction);
	
	///Write configuration
	Script<double,Blob <double> > output(name,ios::out,&System);
	output.write();
	
	delete Umin,Umax,lbond,kbond,kbend,abend;
	return 0;
}


