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

//Include files from the library:

//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"

//For certain lipid models
#include "../models/lipidModels.h"

//For certain nanoparticle models
#include "../models/nanoModels.h"

//Because I'm lazy
using namespace std;

#define CUTOFF 2.0
#define RMIN 1.0

#define HEAD 2
#define TAIL 3
#define HEAD2 4
#define TAIL2 5

int main(int argc, char **argv)
{
	if(argc<7)
	{
		//9, bilayer
		//cout << "Usage: " << argv[0] << " name seed nLipids arealDensity nMonomers nAnchors.x nAnchors.y anchorHeadUmin" << endl;
		
		//7, vesicle
		//cout << "Usage: " << argv[0] << " name seed anchorHeadUmin nLipids arealDensity nMonomers nTess\n";
		
		//7, bilayer, nanoparticle
		cout << "Usage: " << argv[0] << " name seed anchorHeadUmin nLipids lipidArealDensity nanoRadius nanoHeadUmin\n";
		return 0;
	}

	char *name=argv[1];
	double anchorHeadUmin=atof(argv[3]);
	int nLipids=atoi(argv[4]);
	double arealDensity=atof(argv[5]);
	//int nMonomers=atof(argv[6]);
	//int nTess=atof(argv[7]);
	double radius=atof(argv[6]);
	double UminNanoHead=atof(argv[7]);
	//twoVector<int> nAnchors;
	//nAnchors.x=atof(argv[6]);
	//nAnchors.y=atof(argv[7]);
	
	//the variables for the simulation
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
	System.setFinalTime(100);
	System.setDeltaT(0.01);
	System.setStoreInterval(10);
	System.setMeasureInterval(1);
	//System.setDeltaLXY(0.01);//if you don't set it, it doesn't activate
	System.setInitialTemp(3.0);
	System.setFinalTemp(3.0);
	//System.setSolventGamma(5.0);
	
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
	
	//tail types hydrophobicity
	Umin[TAIL+TAIL*System.readNTypes()]=-6;
	Umax[TAIL+TAIL*System.readNTypes()]=200;
	
	Umin[4+HEAD*System.readNTypes()]=UminNanoHead;
	Umax[4+HEAD*System.readNTypes()]=200;
	
	Umin[HEAD+4*System.readNTypes()]=Umin[4+HEAD*System.readNTypes()];
	Umax[HEAD+4*System.readNTypes()]=Umax[4+HEAD*System.readNTypes()];
	
	
	//Umin[TAIL2+TAIL2*System.readNTypes()]=-6;
	//Umax[TAIL2+TAIL2*System.readNTypes()]=200;
	
	//Umin[TAIL2+TAIL*System.readNTypes()]=-5.5;
	//Umax[TAIL2+TAIL*System.readNTypes()]=200;
	
	//Umin[TAIL+TAIL2*System.readNTypes()]=-5.5;
	//Umax[TAIL+TAIL2*System.readNTypes()]=200;
	
	/*
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
	
	//to prevent tangling of anchors
	Umin[ANCHOR+ANCHOR*System.readNTypes()]=0;
	Umax[ANCHOR+ANCHOR*System.readNTypes()]=300;
	
	//anchor attachement without transmembrane proteins
	Umin[HEAD+ANCHOR*System.readNTypes()]=anchorHeadUmin;
	Umax[HEAD+ANCHOR*System.readNTypes()]=200;
	
	Umin[ANCHOR+HEAD*System.readNTypes()]=Umin[HEAD+ANCHOR*System.readNTypes()];
	Umax[ANCHOR+HEAD*System.readNTypes()]=Umax[HEAD+ANCHOR*System.readNTypes()];
	*/
	
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
	
	//double radius=sqrt(((double)nLipids/arealDensity)/(4.0*M_PI));//density adjusted
	threeVector<double>size;
	size=0;
	size.z=40;
	System.setSize(size);
	
	threeVector<double> pos;
	pos.x=0;//System.readSize().x/2.0;
	pos.y=0;//System.readSize().y/2.0;
	pos.z=10;//System.readSize().z/2.0;
	double bondLength=0.7;
	
	double *constants=new double[4];
	constants[0]=bondLength;
	constants[1]=100;
	constants[2]=1;
	constants[3]=100;
	
	//bilayer with cytoskeleton, regular hexagonal cytoskeleton has an aspect ratio sqrt(3/2)=y/x
	bilayer<double>(System, nLipids, 3, pos, bondLength, arealDensity, constants, 4, 1.0);
	//bilayer<double>(System, nLipids, 3, pos, bondLength, arealDensity, constants, 4, (double(nAnchors.x)/double(nAnchors.y))/HEXAGONAL_ASPECT_RATIO);
	//pos.z-=bondLength;
	//flatHexagonalCyto<double>(System, nAnchors, nMonomers, pos, constants, 4);
	
	//liposome with cytoskeleton
	//radius=liposome<double>(System, nLipids, 3, pos, bondLength, arealDensity, constants,4);
	//constants[0]=1.4;//we are adjusting this to prevent blebbing for a while
	//sphericalCyto<double>(System, nMonomers, pos, radius, nTess, constants, 4, System.readNMolecules()-1);
	
	//add solvent
	//solventFill<double>(System, 2.0, SOLVENT_FLAG);
	
	
	//Add a nanoparticle
	int nanoType=4;
	pos.x=System.readSize().x/2.0;
	pos.y=System.readSize().y/2.0;
	pos.z+=radius+1.0+3.0*2.0*0.7;
	
	position<double> *unitCell=new position<double>[4];
	unitCell[0].x=0;
	unitCell[0].y=0;
	unitCell[0].z=0;
	unitCell[0].type=nanoType;
	unitCell[1].x=0.5;
	unitCell[1].y=0.5;
	unitCell[1].z=0;
	unitCell[1].type=nanoType;
	unitCell[2].x=0.5;
	unitCell[2].y=0;
	unitCell[2].z=0.5;
	unitCell[2].type=nanoType;
	unitCell[3].x=0;
	unitCell[3].y=0.5;
	unitCell[3].z=0.5;
	unitCell[3].type=nanoType;
	
	position<double> *bonded=new position<double>[12];
	bonded[0].x=0.5;
	bonded[0].y=0.5;
	bonded[0].z=0;
	bonded[0].type=nanoType;
	bonded[1].x=0.5;
	bonded[1].y=0;
	bonded[1].z=0.5;
	bonded[1].type=nanoType;
	bonded[2].x=0;
	bonded[2].y=0.5;
	bonded[2].z=0.5;
	bonded[2].type=nanoType;
	bonded[3].x=-0.5;
	bonded[3].y=-0.5;
	bonded[3].z=0;
	bonded[3].type=nanoType;
	bonded[4].x=-0.5;
	bonded[4].y=0;
	bonded[4].z=-0.5;
	bonded[4].type=nanoType;
	bonded[5].x=0;
	bonded[5].y=-0.5;
	bonded[5].z=-0.5;
	bonded[5].type=nanoType;
	bonded[6].x=0.5;
	bonded[6].y=-0.5;
	bonded[6].z=0;
	bonded[6].type=nanoType;
	bonded[7].x=0.5;
	bonded[7].y=0;
	bonded[7].z=-0.5;
	bonded[7].type=nanoType;
	bonded[8].x=0;
	bonded[8].y=0.5;
	bonded[8].z=-0.5;
	bonded[8].type=nanoType;
	bonded[9].x=-0.5;
	bonded[9].y=0.5;
	bonded[9].z=0;
	bonded[9].type=nanoType;
	bonded[10].x=-0.5;
	bonded[10].y=0;
	bonded[10].z=0.5;
	bonded[10].type=nanoType;
	bonded[11].x=0;
	bonded[11].y=-0.5;
	bonded[11].z=0.5;
	bonded[11].type=nanoType;
	
	constants[0]=sqrt(2.0)/4.0;//bondlength=latticeLength*sqrt(2.0)/2.0 where latticeLength=0.5
	constants[1]=2800;//kbond
	
	nanoSphere<double>(System, pos, radius, 0.5*1.31345933134, unitCell, 4, bonded, 12, constants, 2);
	
	std::cout << "Storing configuration...";
	
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
	
	std::cout << "Done.\nExiting...\n";
	
	delete constants;
	delete Umin,Umax;
	delete bonded,unitCell;
	return 0;
}


