/**
 * \brief Simple implicit molecular dynamic system. Just a bilayer in various geometries.
 * The steps to make a system are pretty simple, especially given the lipidModels.h header.
 * Just put together the parameters for initial conditions, generate the geometry, and
 * then run the system. That's it! If you want more, less, or something different, just
 * modify this program. There is a section that takes the command line arguments (just
 * after the main function definition) and checks that there are enough. After that,
 * atoi() and atof() (see some cstdlib documentation) convert the text to values. Those
 * are basically your initial conditions. The Blob class template is also useful for this.
 * Just take a parameter you want to modify and either use setParameter, addParameter, or
 * delParameter members in Blob to modify it. I've used the variable definition
 * 'Blob\<double\> System' for convienience. Functions in lipidModels.h accept System as
 * a reference, and can modify parameters as one would in main. For convienience, this
 * program outputs a name.xyz file to quickly view the initial geometry; this is useful
 * for debugging.
 */

#include <iostream>
#include <cstdlib>

//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"
#include "../models/lipidModels.h"
using namespace std;

#define CUTOFF 2.0
#define RMIN 1.0

#define NANOPARTICLE 4

int main(int argc, char **argv)
{
	if(argc<14)
	{
		cout << "Usage: " << argv[0] << " name seed nLipids arealDensity akbend bkbend iRatio oRatio lipidLength initialTemp finalTemp tempStepInterval finalTime" << endl;
		return 0;
	}

	char *name=argv[1];
	int nLipids=atoi(argv[3]);
	double arealDensity=atof(argv[4]);
	double akbend=atof(argv[5]);
	double bkbend=atof(argv[6]);
	double iRatio=atof(argv[7]);
	double oRatio=atof(argv[8]);
	int lipidLength=atoi(argv[9]);

	//the variables for the simulation
	Blob<double> System;
	
	System.setGamma(1.0);
	System.setNTypes(5);
	System.setSeed(atoi(argv[2]));
	
	threeVector<bool> wrap;//periodic boundaries flag, this doesn't work yet, it is always on
	wrap.x=true;
	wrap.y=true;
	wrap.z=true;
	System.setPeriodic(wrap);
	System.setCutoff(2.0);
	
	System.setInitialTime(0);
	System.setFinalTime(atof(argv[13]));
	System.setDeltaT(0.02);
	System.setStoreInterval(10);
	System.setMeasureInterval(1);
	//System.setDeltaLXY(0.01);//if you don't set it, it doesn't activate
	System.setInitialTemp(atof(argv[10]));
	System.setFinalTemp(atof(argv[11]));
	System.setTempStepInterval(atof(argv[12]));
	
	//initialize constants (force and potential constants)
	//generate force constants
	double *Umin=new double[System.readNTypes()*System.readNTypes()];
	double *Umax=new double[System.readNTypes()*System.readNTypes()];
	
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
		}
	}
	
	//umin and umax exceptions, tail types
	Umin[TAILA+TAILA*System.readNTypes()]=-6;
	Umax[TAILA+TAILA*System.readNTypes()]=200;
	
	Umin[TAILB+TAILB*System.readNTypes()]=-6;
	Umax[TAILB+TAILB*System.readNTypes()]=200;
	
	Umin[TAILA+TAILB*System.readNTypes()]=-6;
	Umax[TAILA+TAILB*System.readNTypes()]=200;
	
	Umin[TAILB+TAILA*System.readNTypes()]=Umin[TAILA+TAILB*System.readNTypes()];
	Umax[TAILB+TAILA*System.readNTypes()]=Umax[TAILA+TAILB*System.readNTypes()];
	
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
	
	/*
	MTRand *randNum=new MTRand(System.readSeed());
	//Velocities
	T Vrms=sqrt(3.0*System.readInitialTemp());
	
	for(int i=0;i<20;i++)
	{
		position<double> p;
		threeVector<double> a,v;
		//inner monomers to outer monomers
		p.x=0.01;
		p.y=(double)i;
		p.z=7.5;
		p.type=MONOMER;
		//velocity
		double theta=M_PI*randNum->rand53();
		double phi=M_PI*2.0*randNum->rand53();
		v.x=Vrms*cos(phi)*sin(theta);
		v.y=Vrms*sin(phi)*sin(theta);
		v.z=Vrms*cos(theta);
		System.addParticle(p,v,a);
	}
	delete randNum;
	*/
	
	double radius=sqrt(((double)nLipids/arealDensity)/(4.0*M_PI));//density adjusted
	threeVector<double>size;
	size=radius*4.0;
	System.setSize(size);
	
	threeVector<double> pos;
	pos.x=System.readSize().x/2.0;
	pos.y=System.readSize().y/2.0;
	pos.z=System.readSize().z/2.0;
	double bondLength=0.7;
	
	double *aConstants=new double[4];
	double *bConstants=new double[4];
	aConstants[0]=bondLength;
	aConstants[1]=100;
	aConstants[2]=1;
	aConstants[3]=akbend;
	
	bConstants[0]=bondLength;
	bConstants[1]=100;
	bConstants[2]=1;
	bConstants[3]=bkbend;
	
	int *types=new int[4];
	types[0]=TAILA;
	types[1]=HEADA;
	types[2]=TAILB;
	types[3]=HEADB;
	
	//liposome with cytoskeleton
	radius=liposome<double>(System, nLipids, lipidLength, pos, 0.7, arealDensity, iRatio, oRatio, types, aConstants, bConstants, 4);
	//sphericalCyto<double>(System, nMonomers, pos, radius, nTess, constants, 4, System.readNMolecules()-1);
	
	//add solvent
	//solventFill<double>(System, 2.0, SOLVENT_FLAG);
	
	//Write configuration
	Script<double,Blob <double> > output(name,ios::out,&System);
	output.write();
	
	//This stores the xyz file to check new configuration
	position<double> *p=System.getPositions();
	int nParticles=System.readNParticles();
	
	std::string newName("");
	newName+=name;
	newName+=".xyz";
	
	xyzFormat<double> xyzFile(p, nParticles);
	xyzFile.open(newName.c_str(), std::ios::out);
	xyzFile.store();
	xyzFile.close();
	
	delete aConstants,bConstants;
	delete Umin,Umax;
	delete types;
	return 0;
}


