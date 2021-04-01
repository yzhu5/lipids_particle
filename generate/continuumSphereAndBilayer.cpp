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
#include "../models/nanoModels.h"
using namespace std;

#define CUTOFF 2.0
#define RMIN 1.0

#define HEAD 2
#define TAIL 3
#define HEAD2 6
#define TAIL2 7

#define NANOTYPE 4
//#define SINGLEBEAD 5

int main(int argc, char **argv)
{
	if(argc!=14 && argc!=10 && argc!=6 && argc!=11 && argc!=13)
	{
		std::cout << "nArgs: " << argc << std::endl;
		//14, bilayer with cytoskeleton and nanoparticle
		cout << "Usage: " << argv[0] << " name seed nLipids lipidArealDensity ";
		cout << "nanoRadius nanoHeadUmin nanoHeadUmax sigma nNanoparticles";
		cout << "nMonomers nAnchors.x nAnchors.y anchorHeadUmin " << endl;
		cout << "\tFor bilayer with cytoskeleton and nanoparticles near surface." << endl << endl;
		
		//6, nanoparticles only
		cout << "Usage: " << argv[0] << " name seed nNanoparticles nanoRadius nanoDensity\n";
		cout << "\tFor nanoparticles only." << endl << endl;
		
		//10, bilayer, nanoparticle
		cout << "Usage: " << argv[0] << " name seed nLipids lipidArealDensity ";
		cout << "nanoRadius nanoHeadUmin nanoHeadUmax sigma nNanoparticles" << endl;
		cout << "\tFor bilayer and nanoparticles." << endl << endl;
		
		//11, bilayer, nanoparticle, specific spacing
		cout << "Usage: " << argv[0] << " name seed nLipids lipidArealDensity ";
		cout << "nanoRadius nanoHeadUmin nanoHeadUmax sigma nNanoparticles dNanoparticles" << endl;
		cout << "\tFor bilayer and nanoparticles with a specific spacing (chain configuration)." << endl << endl;
		
		//13, bilayer, nanoparticle, specific spacing
		cout << "Usage: " << argv[0] << " name seed nLipids lipidArealDensity ";
		cout << "nanoRadius nanoHeadUmin nanoHeadUmax sigma nNanoparticles dNanoparticles";
		cout << " boundaryWidth boundaryK" << endl;
		cout << "\tFor bilayer and nanoparticles with a specific spacing (linear configuration)." << endl << endl;
		return 0;
	}
	
	char *name=argv[1];
	std::stringstream cmdArg;
	for(int i=2;i<argc;i++)
		cmdArg << argv[i] << '\t';
	
	int seed, nLipids=0, nMonomers=0, nNanoparticles;
	double lipidArealDensity=0, anchorHeadUmin=0, nanoRadius=0, nanoDensity=0, dNanoparticles=0, nanoHeadUmin=0;
	double nanoNanoUmin=0, nanoHeadUmax=0, sigma=0;
	double boundaryWidth=0, boundaryK=0;
	twoVector<int> nAnchors;
	
	if(argc>=3)
		cmdArg >> seed;
	if(argc==6)
		cmdArg >> nNanoparticles >> nanoRadius >> nanoDensity;
	else if(argc>=8)
	{
		cmdArg >> nLipids >> lipidArealDensity >> nanoRadius >> nanoHeadUmin >> nanoHeadUmax >> sigma >> nNanoparticles;
		if(argc==11 || argc==13)
			cmdArg >> dNanoparticles;
		if(argc==13)
			cmdArg >> boundaryWidth >> boundaryK;
		if(argc==14)
			cmdArg >> nMonomers >> nAnchors.x >> nAnchors.y >> anchorHeadUmin;
	}
	
	//the variables for the simulation
	Blob<double> System;
	
	System.setGamma(1.0);
	System.setNTypes(8);
	System.setSeed(seed);
	
	threeVector<bool> wrap;//periodic boundaries flag, be forwarned, this doesn't work yet, it is always on
	wrap.x=true;
	wrap.y=true;
	wrap.z=true;
	System.setPeriodic(wrap);
	System.setCutoff(2.0);
	
	System.setInitialTime(0);
	System.setFinalTime(200000);
	System.setDeltaT(0.02);
	System.setStoreInterval(1000);
	System.setMeasureInterval(100);
	System.setDeltaLXY(0.02);//if you don't set it, it doesn't activate
	//System.setTension(1.0);
	System.setInitialTemp(3.0);
	System.setFinalTemp(3.0);
	//System.setSolventGamma(5.0);
	
	//initialize constants (force and potential constants)
	//generate force constants
        std::vector<double> Umin(System.readNTypes()*System.readNTypes(),0.0);
        std::vector<double> Umax(System.readNTypes()*System.readNTypes(),0.0);
	//double *Umin=new double[System.readNTypes()*System.readNTypes()];
	//double *Umax=new double[System.readNTypes()*System.readNTypes()];
	
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
        
     /*   
	//Nanoparticle-nanoparticle attraction
	for(int i=4;i<4+nNanoparticles;i++)
	{
		for(int j=4;j<4+nNanoparticles;j++)
		{
			if(i!=j)
			{
				Umin[i+j*System.readNTypes()]=nanoHeadUmin;
				Umax[i+j*System.readNTypes()]=50;
			}
		}
	
		Umin[i+HEAD*System.readNTypes()]=nanoHeadUmin;
		Umax[i+HEAD*System.readNTypes()]=50;
		
		Umin[HEAD+i*System.readNTypes()]=Umin[i+HEAD*System.readNTypes()];
		Umax[HEAD+i*System.readNTypes()]=Umax[i+HEAD*System.readNTypes()];
	}	
	*/
	Umin[TAIL2+TAIL2*System.readNTypes()]=-6;
	Umax[TAIL2+TAIL2*System.readNTypes()]=200;
	
	Umin[TAIL2+TAIL*System.readNTypes()]=-6;
	Umax[TAIL2+TAIL*System.readNTypes()]=200;
	
	Umin[TAIL+TAIL2*System.readNTypes()]=Umin[TAIL2+TAIL*System.readNTypes()];
	Umax[TAIL+TAIL2*System.readNTypes()]=Umax[TAIL2+TAIL*System.readNTypes()];
	
	Umin[HEAD2+HEAD*System.readNTypes()]=-6;
	Umax[HEAD2+HEAD*System.readNTypes()]=200;
	
	Umin[HEAD+HEAD2*System.readNTypes()]=Umin[HEAD2+HEAD*System.readNTypes()];
	Umax[HEAD+HEAD2*System.readNTypes()]=Umax[HEAD2+HEAD*System.readNTypes()];
	
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
			
			//this might be changed to modify the cutoff and such later
			//#ifndef SINGLEBEAD
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
			//#else
			//	//F constants, force constants
			//	System.addTwoBodyFconst(RMIN);//C8
			//	System.addTwoBodyFconst((2.0*(Umax[g]-Umin[g]))/(RMIN*RMIN));//C6
			//	System.addTwoBodyFconst(0);//part of index trick
			//	System.addTwoBodyFconst(CUTOFF);//C7
			//	System.addTwoBodyFconst((6.0*Umin[g])/((CUTOFF-RMIN)*(CUTOFF-RMIN)));//C3
			//	System.addTwoBodyFconst((6.0*Umin[g])/((CUTOFF-RMIN)*(CUTOFF-RMIN)*(CUTOFF-RMIN)));//C2
			//	
			//	//U constants, potential constants
			//	System.addTwoBodyUconst(RMIN);//C8
			//	System.addTwoBodyUconst((Umax[g]-Umin[g])/(RMIN*RMIN));//C4
			//	System.addTwoBodyUconst(Umin[g]);//C5,no index trick
			//	System.addTwoBodyUconst(CUTOFF);//C7
			//	System.addTwoBodyUconst((3.0*Umin[g])/((CUTOFF-RMIN)*(CUTOFF-RMIN)));//C1
			//	System.addTwoBodyUconst((2.0*Umin[g])/((CUTOFF-RMIN)*(CUTOFF-RMIN)*(CUTOFF-RMIN)));//C0
			//#endif
			
		}
	}
	
	//double radius=sqrt(((double)nLipids/lipidArealDensity)/(4.0*M_PI));//density adjusted
	threeVector<double>size;
	size=0;//80;
	size.z=160;
	
	if(nanoDensity>0)
	{
		double volume=nNanoparticles/nanoDensity;
		size.x=pow(volume,1.0/3.0);
		size.y=pow(volume,1.0/3.0);
		size.z=pow(volume,1.0/3.0);
		std::cerr << "Setting system size: " << size.x << '\t' << size.y << '\t' << size.z << std::endl;
	}
	
	System.setSize(size);
	
	threeVector<double> pos;
	pos.x=0;//System.readSize().x/2.0;
	pos.y=0;//System.readSize().y/2.0;
	pos.z=System.readSize().z/2.0;
	double bondLength=0.7;
	
	double *constants=new double[4];
	constants[0]=bondLength;
	constants[1]=100;
	constants[2]=1;
	constants[3]=100;
	
	//bilayer with cytoskeleton, regular hexagonal cytoskeleton has an aspect ratio sqrt(3/2)=y/x
	if(nLipids>0)
	{
		if(nAnchors.x<=0 || nAnchors.y<=0 || nMonomers<=0)
			bilayer<double>(System, nLipids, 3, pos, bondLength, lipidArealDensity, constants, 4, 1.0);
		else
			bilayer<double>(System, nLipids, 3, pos, bondLength, lipidArealDensity, constants, 4, (double(nAnchors.x)/double(nAnchors.y))/HEXAGONAL_ASPECT_RATIO);
	}
	//pos.z-=bondLength;
	if(nAnchors.x>0 && nAnchors.y>0 && nMonomers>0)
		flatHexagonalCyto<double>(System, nAnchors, nMonomers, pos, constants, 4);
	
	//liposome with cytoskeleton
	//radius=liposome<double>(System, nLipids, 3, pos, bondLength, lipidArealDensity, constants,4);
	//constants[0]=1.4;//we are adjusting this to prevent blebbing for a while
	//sphericalCyto<double>(System, nMonomers, pos, radius, nTess, constants, 4, System.readNMolecules()-1);
	
	//add solvent
	//solventFill<double>(System, 2.0, SOLVENT_FLAG);
	
	
	//Add a nanoparticle
	//int NANOTYPE=4;//could use different nanotypes at different lattice positions
	int alterType=5;
	
	//double sigma=73.0;
	
	//constants for bead type
	std::vector<double> C;
	//defaults, 200=Umax, 0=Umin, rc=2, rmin=1, sigma=46
	double rminNP=0.1;
	double cutoffNP=0.2;
	for(int typeA=0;typeA<System.readNTypes();typeA++)
	{
		for(int typeB=0;typeB<System.readNTypes();typeB++)
		{
			//double rmin=1.0;
			//double rc=2.0*rmin;
			double UminNano=0.0;
			double UmaxNano=10.0;
			
			//bead-particle constants
			C.push_back(cutoffNP+nanoRadius);//0
			C.push_back(-7.0/4.0*rminNP);//1
			C.push_back(2.0*rminNP*rminNP);//2
			C.push_back(UminNano*M_PI*nanoRadius*sigma/(rminNP*rminNP*rminNP));//3
			C.push_back(nanoRadius);//4
			C.push_back(rminNP);//5
			double D=sigma*M_PI*nanoRadius;
			double A=UmaxNano-UminNano;
			
			C.push_back(-D*A/(2.0*rminNP*rminNP));//6,0@B^4
			C.push_back(2.0*D*A/(3.0*rminNP));//7,1,@B^3
			C.push_back(-D*UminNano);//8,2,@B^2
			C.push_back(2.0*D*UminNano*rminNP);//9,3,@B^1
			C.push_back(D*1.3*UminNano*rminNP*rminNP);//10,4
			
			//bead-bead constants
			D=pow(2.0*M_PI*sigma*nanoRadius,2.0)/(rminNP*rminNP*rminNP);
			
			C.push_back(2.0*nanoRadius+rminNP);//0
			C.push_back(2.0*nanoRadius+cutoffNP);//1
			C.push_back(UminNano*D*2.0/30.0);//2, x^6
			C.push_back(-UminNano*D*7.0*rminNP/20.0);//3, x^5
			C.push_back(UminNano*D*rminNP*rminNP/2.0);//4, x^4
			
			C.push_back(-D*A*rminNP/20.0);//5,0@B^5
			C.push_back(D*A*rminNP*rminNP/12.0);//6,1,@B^4
			C.push_back(-D*UminNano*pow(rminNP,3.0)/6.0);//7,2,@B^3
			C.push_back(D*UminNano*pow(rminNP,4.0)/2.0);//8,3,@B^2
			C.push_back(D*1.3*UminNano*pow(rminNP,5.0)/2.0);//9,@B
			C.push_back(D*13.0*UminNano*pow(rminNP,6.0)/60.0);//10
		}
	}
	
	//exceptions
	//nanoparticle to head
	if(nNanoparticles>0)
	{
		int cindex=nBEADCONST*(HEAD+NANOTYPE*System.readNTypes());
		int cindexMirror=nBEADCONST*(NANOTYPE+HEAD*System.readNTypes());
		
		double UminNano=nanoHeadUmin;
		double UmaxNano=nanoHeadUmax;
		
		//bead-particle constants
		C[cindex+0]=(cutoffNP+nanoRadius);//c0
		C[cindex+1]=(-7.0/4.0*rminNP);//c1
		C[cindex+2]=(2.0*rminNP*rminNP);//c2
		C[cindex+3]=(UminNano*M_PI*nanoRadius*sigma/(rminNP*rminNP*rminNP));//c3
		C[cindex+4]=(nanoRadius);//c4
		C[cindex+5]=(rminNP);//c5
		
		double D=sigma*M_PI*nanoRadius;
		double A=UmaxNano-UminNano;
			
		C[cindex+6]=-D*A/(2.0*rminNP*rminNP);//c6
		C[cindex+7]=2.0*D*A/(3.0*rminNP);//c7
		C[cindex+8]=-D*UminNano;
		C[cindex+9]=2.0*D*UminNano*rminNP;
		C[cindex+10]=D*1.3*UminNano*rminNP*rminNP;
		
		//bead-bead constants
		D=pow(2.0*M_PI*sigma*nanoRadius,2.0)/(rminNP*rminNP*rminNP);
		
		C[cindex+11]=(2.0*nanoRadius+rminNP);//0
		C[cindex+12]=(2.0*nanoRadius+cutoffNP);//1
		C[cindex+13]=(UminNano*D*2.0/30.0);//2, x^6
		C[cindex+14]=(-UminNano*D*7.0*rminNP/20.0);//3, x^5
		C[cindex+15]=(UminNano*D*rminNP*rminNP/2.0);//4, x^4
		
		C[cindex+16]=(-D*A*rminNP/20.0);//5,0@B^5
		C[cindex+17]=(D*A*rminNP*rminNP/12.0);//6,1,@B^4
		C[cindex+18]=(-D*UminNano*pow(rminNP,3.0)/6.0);//7,2,@B^3
		C[cindex+19]=(D*UminNano*pow(rminNP,4.0)/2.0);//8,3,@B^2
		C[cindex+20]=(D*1.3*UminNano*pow(rminNP,5.0)/2.0);//9,@B
		C[cindex+21]=(D*13.0*UminNano*pow(rminNP,6.0)/60.0);//10
		
		for(int i=0;i<nBEADCONST;i++)
			C[cindexMirror+i]=C[cindex+i];
	}
	
	std::vector< threeVector<double> > nanoPos;
	MTRand randNum(System.readSeed());
	
	int nanoOffset=System.readNParticles();
	
	for(int i=0;i<nNanoparticles;i++)
	{
		std::cout << "Placing nanoparticle " << i << "!" << std::endl;
		threeVector<double> toPlace;
		bool overlap;
		do
		{
			overlap=false;
			if(dNanoparticles<=0)//no chain
			{
				toPlace.x=(System.readSize().x-2.0*nanoRadius)*randNum.rand53()+nanoRadius;
				toPlace.y=(System.readSize().y-2.0*nanoRadius)*randNum.rand53()+nanoRadius;
				if(nanoDensity<=0)//bilayer
					toPlace.z=pos.z+nanoRadius+rminNP+3.0*2.0*0.7;
				else//no bilayer
					toPlace.z=(System.readSize().y-2.0*nanoRadius)*randNum.rand53()+nanoRadius;
			}
			else if(boundaryWidth<=0)//chain placement
			{
				double xPos, yPos;
				bool outOfBounds;
				if(i!=0)
				{
					do
					{
						outOfBounds=false;
						double theta=2.0*M_PI*randNum.rand53();
						xPos=sin(theta)*(dNanoparticles+2.0*nanoRadius);
						yPos=cos(theta)*(dNanoparticles+2.0*nanoRadius);
						if(nanoPos[i-1].x+xPos>System.readSize().x-nanoRadius)
							outOfBounds=true;
						if(nanoPos[i-1].y+yPos>System.readSize().y-nanoRadius)
							outOfBounds=true;
						if(nanoPos[i-1].x+xPos<nanoRadius)
							outOfBounds=true;
						if(nanoPos[i-1].y+yPos<nanoRadius)
							outOfBounds=true;
					} while(outOfBounds);
					
					toPlace.x=xPos+nanoPos[i-1].x;
					toPlace.y=yPos+nanoPos[i-1].y;
					toPlace.z=pos.z+nanoRadius+1.0+3.0*2.0*0.7;
				}
				else
				{
					toPlace.x=(System.readSize().x-2.0*nanoRadius)*randNum.rand53()+nanoRadius;
					toPlace.y=(System.readSize().y-2.0*nanoRadius)*randNum.rand53()+nanoRadius;
					toPlace.z=pos.z+nanoRadius+1.0+3.0*2.0*0.7;
				}
			}
			else//linear placement, center of system
			{
				if(nanoPos.size()>0)
					toPlace.x=2.0*nanoRadius+dNanoparticles+nanoPos[i-1].x;
				else
					toPlace.x=System.readSize().x/2.0-(dNanoparticles/2.0+nanoRadius)*static_cast<double>(nNanoparticles);
				toPlace.y=System.readSize().y/2.0;
				toPlace.z=pos.z+nanoRadius+1.0+3.0*2.0*0.7;
			}
			for(int j=0;j<nanoPos.size();j++)
			{
				threeVector<double> d;
				d.x=nanoPos[j].x-toPlace.x;
				d.y=nanoPos[j].y-toPlace.y;
				d.z=nanoPos[j].z-toPlace.z;
				d.x-=(d.x>System.readSize().x/2.0)?System.readSize().x:0;
				d.x+=(d.x<-System.readSize().x/2.0)?System.readSize().x:0;
				
				d.y-=(d.y>System.readSize().y/2.0)?System.readSize().y:0;
				d.y+=(d.y<-System.readSize().y/2.0)?System.readSize().y:0;
				
				d.z-=(d.z>System.readSize().z/2.0)?System.readSize().z:0;
				d.z+=(d.z<-System.readSize().z/2.0)?System.readSize().z:0;
				
				if(sqrt(d.x*d.x+d.y*d.y+d.z*d.z)<2.0*nanoRadius+0.5)
					overlap=true;
			}
		} while(overlap);
		std::cerr << i << '\t' << toPlace.x << '\t' << toPlace.y << '\t' << toPlace.z << endl;
		nanoPos.push_back(toPlace);
		
		//continuumSphere<double>(System, nanoPos[i], nanoRadius, &C[0], System.readNTypes()*System.readNTypes()*nBEADCONST, NANOTYPE);
	}
	continuumSphere<double>(System, &(nanoPos[0]), nanoPos.size(), nanoRadius, &C[0], 
				System.readNTypes()*System.readNTypes()*nBEADCONST, NANOTYPE);
	if(boundaryWidth>0)
	{
		molecule< double, fourVector<int> > boundaryY;
		boundaryY.setType(BOUNDARY);
		for(int i=nanoOffset;i<nanoOffset+nNanoparticles;i++)
		{
			fourVector<int> buf;
			buf.s[0]=i;
			boundaryY.addBond(buf);
		}
		boundaryY.addConstant(1);//along y
		boundaryY.addConstant(0.5);//center of y
		boundaryY.addConstant(boundaryWidth/2.0);//
		boundaryY.addConstant(boundaryK/2.0);
		
		System.addMolecule(boundaryY);
	}
	
	std::cerr << "Storing configuration...";
	
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
	
	std::cerr << "Done.\nExiting...\n";
	return 0;
}


