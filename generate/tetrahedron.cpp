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

int main(int argc, char **argv)
{
	if(argc<3)
	{
		cout << "Usage: " << argv[0] << " name seed\n";
		return 0;
	}

	char *name=argv[1];
	
	//the variables for the simulation
	Blob<double> System;
	
	System.setGamma(1.0);
	System.setNTypes(2);
	System.setSeed(atoi(argv[2]));
	
	threeVector<bool> wrap;//periodic boundaries flag, be forwarned, this doesn't work yet, it is always on
	wrap.x=true;
	wrap.y=true;
	wrap.z=true;
	System.setPeriodic(wrap);
	System.setCutoff(2.0);
	
	System.setInitialTime(0);
	System.setFinalTime(1000);
	System.setDeltaT(0.02);
	System.setStoreInterval(0.02);
	System.setMeasureInterval(0.02);
	//System.setDeltaLXY(0.01);//if you don't set it, it doesn't activate
	System.setInitialTemp(3.0);
	System.setFinalTemp(3.0);
	
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
	
	threeVector<double>size;
	size.x=50;
	size.y=50;
	size.z=50;
	System.setSize(size);
	
	//Adding a tetrahedron
	MTRand randNum(System.readSeed());
	double Vrms=sqrt(3.0*System.readInitialTemp());
	fourVector<int> bond;
	molecule<double,fourVector<int> > m;
	position<double> vertex;
	vertex.type=1;
	threeVector<double> v,a;
	a=0;
	//in plane
	vertex.x=25;
	vertex.y=25;
	vertex.z=25;
	double theta=M_PI*randNum.rand53();
	double phi=M_PI*2.0*randNum.rand53();
	v.x=Vrms*cos(phi)*sin(theta);
	v.y=Vrms*sin(phi)*sin(theta);
	v.z=Vrms*cos(theta);
	System.addParticle(vertex,v,a);
	vertex.x=25;
	vertex.y=26;
	vertex.z=25;
	theta=M_PI*randNum.rand53();
	phi=M_PI*2.0*randNum.rand53();
	v.x=Vrms*cos(phi)*sin(theta);
	v.y=Vrms*sin(phi)*sin(theta);
	v.z=Vrms*cos(theta);
	System.addParticle(vertex,v,a);
	vertex.x=26;
	vertex.y=25;
	vertex.z=25;
	theta=M_PI*randNum.rand53();
	phi=M_PI*2.0*randNum.rand53();
	v.x=Vrms*cos(phi)*sin(theta);
	v.y=Vrms*sin(phi)*sin(theta);
	v.z=Vrms*cos(theta);
	System.addParticle(vertex,v,a);
	//out of plane
	vertex.x=25.5;
	vertex.y=25.5;
	vertex.z=25.5;
	theta=M_PI*randNum.rand53();
	phi=M_PI*2.0*randNum.rand53();
	v.x=Vrms*cos(phi)*sin(theta);
	v.y=Vrms*sin(phi)*sin(theta);
	v.z=Vrms*cos(theta);
	System.addParticle(vertex,v,a);
	
	m.addConstant(1);
	m.addConstant(100);
	m.setType(BOND);
	bond.s[0]=0;
	bond.s[1]=1;
	m.addBond(bond);
	bond.s[0]=0;
	bond.s[1]=2;
	m.addBond(bond);
	bond.s[0]=0;
	bond.s[1]=3;
	m.addBond(bond);
	bond.s[0]=1;
	bond.s[1]=2;
	m.addBond(bond);
	bond.s[0]=1;
	bond.s[1]=3;
	m.addBond(bond);
	bond.s[0]=2;
	bond.s[1]=3;
	m.addBond(bond);
	System.addMolecule(m);
	
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
	
	delete Umin,Umax;
	return 0;
}
