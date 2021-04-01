/**
 * \brief Various model geometries associated with lipid systems.
 * We have flat and spherical cytoskeletons, solvent fills, and flat and spherical lipid bilayers. 
 * Many options for connectivity as well. Assumes your system blob has some functions to utilize 
 * with molecules, particles, and temperature.
 */

#include "../include/algorithms/functions.h"

//aspectRatio is x/y
template <typename T>
bool brush(Blob<T> &System, int nBrushes, int brushLength, threeVector<T> pos, T bondLength, T arealDensity, T *constants, int nConstants, double aspectRatio, int bottomType, int chainType)
{
	//Error checking
	if(nConstants==4)
		std::cerr << "nConstants is 4, assuming CHAIN type (brush)...\n";
	if(nConstants==2)
		std::cerr << "nConstants is 2, assuming BOND type (brush)...\n";
	
	if(nConstants!=2 && nConstants!=4)
	{
		std::cerr << "Error: nConstants is not 2 or 4![bool brush()]\n";
		return true;
	}
	
	//Initialize molecule
	fourVector<int> bond;
	molecule<T,fourVector<int> > m;
	
	if(nConstants==4)
	{
		m.setType(CHAIN);
		
		//bond for a chain type
		bond.s[START]=System.readNParticles();
		bond.s[NCHAINS]=0;//nBrushes;
		bond.s[CHAINLENGTH]=brushLength;
		//m.addBond(bond);
	}
	
	if(nConstants==2)
	{
		m.setType(BOND);
		m.allocBonds(nBrushes*brushLength);
	}
	
	for(int i=0;i<nConstants;i++)
		m.addConstant(constants[i]);
	
		
	
	//adjust size of system if it is out of bounds
	threeVector<T> s;
	
	s.x=sqrt((T)nBrushes/arealDensity)*sqrt(aspectRatio);
	s.y=sqrt((T)nBrushes/arealDensity)/sqrt(aspectRatio);
	s.z=pos.z;
	std::cerr << "Prefered system size is (brush): " << s.x << '\t' << s.y << '\t' << s.z << '\n';
	
	threeVector<T> size=System.readSize();
	
	std::cerr << "Current system size is (brush): " << size.x << '\t' << size.y << '\t' << size.z << '\n';
	
	//Adjust number of brushes to match grafting density
	
	size.x=(size.x<s.x)?s.x:size.x;
	size.y=(size.y<s.y)?s.y:size.y;
	size.z=(size.z<s.z)?s.z:size.z;
	
	
	MTRand *randNum=new MTRand(System.readSeed());
	//Velocities
	T Vrms=sqrt(3.0*System.readInitialTemp());
	
	//This is important to reduce excessive reallocations
	System.allocParticle(System.readNParticles()+2*nBrushes*brushLength);
	
	threeVector<int> latticeSize;
	threeVector<T> latticeLength;
	
	latticeSize.x=static_cast<int>(sqrt((arealDensity*size.x*size.x)/aspectRatio));
	latticeSize.y=static_cast<int>(sqrt((arealDensity*size.y*size.y)*aspectRatio));
	
	int mBrushes=latticeSize.x*latticeSize.y;
	
	size.x=sqrt((T)mBrushes/arealDensity)*sqrt(aspectRatio);
	size.y=sqrt((T)mBrushes/arealDensity)/sqrt(aspectRatio);
	
	latticeLength.x=size.x/static_cast<T>(latticeSize.x);
	latticeLength.y=size.y/static_cast<T>(latticeSize.y);
		
	System.setSize(size);
	
	std::cerr << "Adjusted system size is (brush): " << size.x << '\t' << size.y << '\t' << size.z << '\n';	
	std::cerr << "Adjusted density is (brush): " << static_cast<T>(latticeSize.x*latticeSize.y)/(size.x*size.y) << '=';
	std::cerr << latticeSize.x << '*' << latticeSize.y << "/(" << size.x << '*' << size.y << ')' << std::endl;
	//latticeLength.x=size.x/static_cast<double>(latticeSize.x);
	//latticeLength.y=size.y/static_cast<double>(latticeSize.y);
	
	//initialize brush positions
	for(int i=0;i<latticeSize.x;i++)
	{
		for(int j=0;j<latticeSize.y;j++)
		{
			int brushIndex=(j*latticeSize.x+i);
			//if(brushIndex%1000==0)
			//	std::cerr << brushIndex << std::endl;
			T theta,phi;
			
			position<T> p;
			threeVector<T> v;
			threeVector<T> a;
			//set all acceleration to 0 initially
			a.x=0;
			a.y=0;
			a.z=0;
			
			for(int k=0;k<brushLength;k++)
			{
				//position
				p.x=i*latticeLength.x+0.001+size.x/2.0;
				p.y=j*latticeLength.y+0.001+size.y/2.0;
				p.z=pos.z+bondLength*(k);
				p.type=(k==0)?bottomType:chainType;
				
				//velocity
				theta=M_PI*randNum->rand53();
				phi=M_PI*2*randNum->rand53();
				v.x=Vrms*cos(phi)*sin(theta);
				v.y=Vrms*sin(phi)*sin(theta);
				v.z=Vrms*cos(theta);
				
				//save it
				System.addParticle(p,v,a);
				
				//if it is of BOND type and not the last monomer
				if(nConstants==2 && k!=brushLength-1)
				{
					//bond
					bond.s[0]=System.readNParticles()-1;
					bond.s[1]=System.readNParticles();
					m.addBond(bond);
				}
			}
			if(nConstants==4)
			{
				bond.s[NCHAINS]++;
			}
		}
	}	
	if(nConstants==4)
	{
		m.addBond(bond);
	}
	System.addMolecule(m);
	
	return true;
}
