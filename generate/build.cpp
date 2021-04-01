#include <iostream>
#include <cstdlib>
//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"
using namespace std;

#define CUTOFF 2.0
#define RMIN 1.0

#define NANOPARTICLE 4

int main(int argc, char **argv)
{
	if(argc<4)
	{
		cout << "Usage: " << argv[0] << " nLipids arial_density name" << endl;
		return 0;
	}

	char *name=argv[3];
	
	//load parameters
	int nLipids=atoi(argv[1]);
	double arialDensity=atof(argv[2]);
	
	///the variables for the simulation
	
	Blob<double> System;
	
	System.setGamma(1.0);
	System.setNTypes(5);
	System.setSeed(1234);
	
	threeVector<bool> wrap;//periodic boundaries flag, be forwarned, this doesn't work yet, it is always on
	wrap.x=true;
	wrap.y=true;
	wrap.z=true;
	System.setPeriodic(wrap);
	System.setCutoff(2.0);
	
	threeVector<double> size;   //component sizes of system
	size.x=sqrt(nLipids/arialDensity);//pow((double)particles/density,1.0/3.0);
	size.y=size.x;
	size.z=100.0;//size.x;
	System.setSize(size);
	
	System.setInitialTime(0);
	System.setFinalTime(10);
	System.setDeltaT(0.02);
	System.setStoreInterval(10);
	System.setMeasureInterval(1);
	//System.setDeltaLXY(0.01);//if you don't set it, it doesn't activate
	System.setInitialTemp(3.0);
	System.setFinalTemp(3.0);
	
	threeVector<int> latticeSize;
	threeVector<double> latticeLength;
	
	latticeSize.x=sqrt((double)nLipids/2.0)+1;//int(pow((double)particles,1.0/3.0))+1;
	latticeSize.y=latticeSize.x;
	latticeSize.z=latticeSize.y;
	
	latticeLength.x=size.x/latticeSize.x;
	latticeLength.y=size.y/latticeSize.y;
	latticeLength.z=size.z/latticeSize.z;
	
	//creating a bunch of molecules with specialized properties
	
	fourVector<int> bond;
	molecule<double,fourVector<int> > m;
	
	m.setType(CHAIN);
	//constants for bond
	m.addConstant(0.7);
	m.addConstant(100);
	m.addConstant(1.0);
	m.addConstant(100);
	//For a chain type bond
	bond.s[START]=0;
	bond.s[NCHAINS]=nLipids;
	bond.s[CHAINLENGTH]=3;
	m.addBond(bond);
	
	System.addMolecule(m);
	
	m.~molecule();
	
	MTRand *randNum=new MTRand(System.readSeed());
	
	//Velocities
	double Vrms=sqrt(3*System.readInitialTemp());
	
	///This is important to reduce excessive allocations
	System.allocParticle(System.readNParticles()+nLipids*3);
	
	//initialize lipid positions
	for(int i=0;i<latticeSize.x;i++)
	{
		for(int j=0;j<latticeSize.y;j++)
		{
			int lipidIndex=2*(j*latticeSize.x+i);
			double theta,phi;
			
			position<double> p;
			threeVector<double> v;
			threeVector<double> a;
			//set all acceleration to 0 initially
			a.x=0;
			a.y=0;
			a.z=0;
			
			if(lipidIndex<nLipids)
			{
				//position
				p.x=i*latticeLength.x+0.001;
				p.y=j*latticeLength.y+0.001;
				p.z=0.001+size.z/2.0-2.501;
				p.type=HEAD;
				//velocity
				theta=M_PI*randNum->rand53();
				phi=M_PI*2*randNum->rand53();
				v.x=Vrms*cos(phi)*sin(theta);
				v.y=Vrms*sin(phi)*sin(theta);
				v.z=Vrms*cos(theta);
				//save it
				System.addParticle(p,v,a);
				
				//position
				p.x=i*latticeLength.x+0.001;
				p.y=j*latticeLength.y+0.001;
				p.z=1.001+size.z/2.0-2.501;//+latticeLength.z;
				p.type=TAIL;
				//velocity
				theta=M_PI*randNum->rand53();
				phi=M_PI*2*randNum->rand53();
				v.x=Vrms*cos(phi)*sin(theta);
				v.y=Vrms*sin(phi)*sin(theta);
				v.z=Vrms*cos(theta);
				//save it
				System.addParticle(p,v,a);
				
				//position
				p.x=i*latticeLength.x+0.001;
				p.y=j*latticeLength.y+0.001;
				p.z=2.001+size.z/2.0-2.501;//+latticeLength.z*2.0;
				p.type=TAIL;
				//velocity
				theta=M_PI*randNum->rand53();
				phi=M_PI*2*randNum->rand53();
				v.x=Vrms*cos(phi)*sin(theta);
				v.y=Vrms*sin(phi)*sin(theta);
				v.z=Vrms*cos(theta);
				//save it
				System.addParticle(p,v,a);
			}
			if(lipidIndex+1<nLipids)
			{
				//position
				p.x=i*latticeLength.x+0.001;
				p.y=j*latticeLength.y+0.001;
				p.z=3.001+size.z/2.0-2.501;
				p.type=TAIL;
				//velocity
				theta=M_PI*randNum->rand53();
				phi=M_PI*2*randNum->rand53();
				v.x=Vrms*cos(phi)*sin(theta);
				v.y=Vrms*sin(phi)*sin(theta);
				v.z=Vrms*cos(theta);
				//save it
				System.addParticle(p,v,a);
				
				//position
				p.x=i*latticeLength.x+0.001;
				p.y=j*latticeLength.y+0.001;
				p.z=4.001+size.z/2.0-2.501;//+latticeLength.z;
				p.type=TAIL;
				//velocity
				theta=M_PI*randNum->rand53();
				phi=M_PI*2*randNum->rand53();
				v.x=Vrms*cos(phi)*sin(theta);
				v.y=Vrms*sin(phi)*sin(theta);
				v.z=Vrms*cos(theta);
				//save it
				System.addParticle(p,v,a);
				
				//position
				p.x=i*latticeLength.x+0.001;
				p.y=j*latticeLength.y+0.001;
				p.z=5.001+size.z/2.0-2.501;//+latticeLength.z*2.0;
				p.type=HEAD;
				//velocity
				theta=M_PI*randNum->rand53();
				phi=M_PI*2*randNum->rand53();
				v.x=Vrms*cos(phi)*sin(theta);
				v.y=Vrms*sin(phi)*sin(theta);
				v.z=Vrms*cos(theta);
				//save it
				System.addParticle(p,v,a);
			}
		}
	}
	
	///Making a cube
	/*{
		int nPCube=0;
		double latticeUnit=0.5;
		//Position particles
		for(int i=0;i<10;i++)
		{
			for(int j=0;j<10;j++)
			{
				for(int k=0;k<10;k++)
				{
					p[nLipids*3+nPCube].type=NANOPARTICLE;
					p[nLipids*3+nPCube].x=latticeUnit*i+2.0;
					p[nLipids*3+nPCube].y=latticeUnit*j+2.0;
					p[nLipids*3+nPCube].z=latticeUnit*k+2.0;
					nPCube++;
				}
			}
		}
		
		m[1].type=BONDED;
		m[1].nBonded=0;
		m[1].mC=new double[nBONDEDCONST];
		m[1].mC[ABOND]=latticeUnit;
		m[1].mC[KBOND]=200;
		
		m[2].type=BENDED;
		m[2].nBonded=0;
		m[2].mC=new double[nBENDEDCONST];
		m[2].mC[ABOND]=0.0;
		m[2].mC[KBOND]=200;
		
		m[3].type=BENDED;
		m[3].nBonded=0;
		m[3].mC=new double[nBENDEDCONST];
		m[3].mC[ABOND]=1.0;
		m[3].mC[KBOND]=200;
		
		//Count bond neighbors
		for(int i=nLipids*3;i<nPCube+nLipids*3;i++)
		{
			for(int j=i+1;j<nPCube+nLipids*3;j++)
			{
				double dxa=p[i].x-p[j].x;
				double dya=p[i].y-p[j].y;
				double dza=p[i].z-p[j].z;
				double dra=sqrt(dxa*dxa+dya*dya+dza*dza);
				if(dra<=latticeUnit*1.1)
				{
					m[1].nBonded++;
					for(int k=j+1;k<nPCube+nLipids*3;k++)
					{
						double dxb=p[j].x-p[k].x;
						double dyb=p[j].y-p[k].y;
						double dzb=p[j].z-p[k].z;
						double drb=sqrt(dxb*dxb+dyb*dyb+dzb*dzb);
						if(drb<=latticeUnit*1.1)
						{
							dxa/=dra;
							dya/=dra;
							dza/=dra;
							
							dxb/=drb;
							dyb/=drb;
							dzb/=drb;
							
							double costheta=(dxa*dxb)+(dya*dyb)+(dza*dzb);
							if(costheta<m[2].mC[ABOND]+0.1 && costheta>m[2].mC[ABOND]-0.1)
							{
								m[2].nBonded++;
								
							}
							std::cout << costheta << '\n';
							if(costheta<m[3].mC[ABOND]+0.1 && costheta>m[3].mC[ABOND]-0.1)
							{
								m[3].nBonded++;
								if(i==nLipids*3+545 || i==nLipids*3+672)
								{
								//	std::cout << 4 << '\t' << p[i].x << '\t' << p[i].y << '\t' << p[i].z << '\n';
								//	std::cout << 5 << '\t' << p[j].x << '\t' << p[j].y << '\t' << p[j].z << '\n';
								//	std::cout << 5 << '\t' << p[k].x << '\t' << p[k].y << '\t' << p[k].z << '\n';
								}
							}
						}
					}
				}
			}
		}
		
		m[3].bond=new genericVector<int,4>[m[3].nBonded];
		m[2].bond=new genericVector<int,4>[m[2].nBonded];
		m[1].bond=new genericVector<int,4>[m[1].nBonded];
		//std::cout << m[1].nBonded << "\nbonds\n";
		
		int pair=0;
		int ninetyTriplet=0;
		int oneEightyTriplet=0;
		//Locate bond neighbors
		for(int i=nLipids*3;i<nPCube+nLipids*3;i++)
		{
			for(int j=i+1;j<nPCube+nLipids*3;j++)
			{
				double dxa=p[i].x-p[j].x;
				double dya=p[i].y-p[j].y;
				double dza=p[i].z-p[j].z;
				double dra=sqrt(dxa*dxa+dya*dya+dza*dza);
				if(dra<=latticeUnit*1.1)
				{
					//std::cout << "1\t";
					//std::cout << (p[i].x+p[j].x)/2.0 << '\t';
					//std::cout << (p[i].y+p[j].y)/2.0 << '\t';
					//std::cout << (p[i].z+p[j].z)/2.0 << '\n';
					m[1].bond[pair].s[0]=i;
					m[1].bond[pair++].s[1]=j;
					for(int k=j+1;k<nPCube+nLipids*3;k++)
					{
						double dxb=p[j].x-p[k].x;
						double dyb=p[j].y-p[k].y;
						double dzb=p[j].z-p[k].z;
						double drb=sqrt(dxb*dxb+dyb*dyb+dzb*dzb);
						if(drb<=latticeUnit*1.1)
						{
							dxa/=dra;
							dya/=dra;
							dza/=dra;
							
							dxb/=drb;
							dyb/=drb;
							dzb/=drb;
							
							double costheta=(dxa*dxb)+(dya*dyb)+(dza*dzb);
							if(costheta<m[2].mC[ABOND]+0.1 && costheta>m[2].mC[ABOND]-0.1)
							{
								m[2].bond[ninetyTriplet].s[0]=i;
								m[2].bond[ninetyTriplet].s[1]=j;
								m[2].bond[ninetyTriplet++].s[2]=k;
							}
							if (costheta<m[3].mC[ABOND]+0.1 && costheta>m[3].mC[ABOND]-0.1)
							{
								m[3].bond[oneEightyTriplet].s[0]=i;
								m[3].bond[oneEightyTriplet].s[1]=j;
								m[3].bond[oneEightyTriplet++].s[2]=k;
							}
						}
					}
				}
			}
		}
	}
	*/
	/*
	///this is really all that's needed to generate a simple molecule
	{
		//creating a molecule
		m[1].bond=new genericVector<int,4>[1];
		m[1].type=CHAIN;
		m[1].bond[0].s[START]=nLipids*3;
		m[1].bond[0].s[NCHAINS]=1;
		m[1].bond[0].s[LENGTH]=20;
	
		for(int i=0;i<20;i++)
		{
			p[i+m[1].bond[0].s[START]].x=(double)i+0.001;
			p[i+m[1].bond[0].s[START]].y=size.y/2.0;
			p[i+m[1].bond[0].s[START]].z=size.z/2.0;
			p[i+m[1].bond[0].s[START]].type=HEAD;
		}
	}
	*/
	
	
	
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
	
	Umin[NANOPARTICLE+NANOPARTICLE*System.readNTypes()]=0;
	Umax[NANOPARTICLE+NANOPARTICLE*System.readNTypes()]=0;
	
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
	
	Script<double,Blob <double> > output(name,ios::out,&System);
	
	output.write();
	
	delete randNum;
	delete Umin,Umax,lbond,kbond,kbend,abend;
	return 0;
}


