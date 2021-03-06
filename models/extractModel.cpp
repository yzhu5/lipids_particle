//This is a very short and simple program that shows you how to modify a variable already present in an mpd script.
//The reason this is used is because some variables (forces, potentials, geometry, etc...) are not easily modified
// in the mpd script. Using a file that opens the mpd script, modifies, and then saves it (as below) is far easier
// than modifying the actual variable. Most of it looks like the run.cpp file, it just replaces the "run loop"
// with the modification.

#define NO_REQUIRED_COMMANDS

#include "../include/systemMD.h"

int main(int argc, char* argv[])
{
	if(argc<=1)
	{
		//so simple, this can't possibly mess it up
		std::cout << "Usage (To List Molecules): " << argv[0] << " name\n";
		std::cout << "Usage (To Extract a Model): " << argv[0] << " name modelName molecule1 molecule2 molecule3 ...\n";
		return 0;
	}
	
	
	//read in options
	char *name=argv[1];
	
	
	//Configuration variables
	Blob<double> System;
	
	
	//load old variables, then initialize them
	Script<double, Blob <double> > fileIO(name,std::ios::in,&System);
	fileIO.read();
	fileIO.close();
	
	
	//Check which molecules are available
	if(argc==2)
	{
		for(int i=0;i<System.readNMolecules();i++)
		{
			molecule< double,fourVector<int> > *m=System.getMolecule();
			int nBonded=m[i].readNBond();
			fourVector<int> *bond=m[i].getBonds();
			int type=m[i].readType();
			switch(type)
			{
				case CHAIN:
					std::cout << "Molecule " << i << " is a chain type.\n";
					std::cout << "\tIt has " << nBonded << " chain lists.\n";
					for(int j=0;j<nBonded;j++)
					{
						std::cout << "\t\tList " << j << " has " << bond[j].s[NCHAINS] << " chains of length " << bond[j].s[LENGTH] 
						<< " starting at " << bond[j].s[START] << '\n';
						std::cout << "\t\t\tStarting index, " << bond[j].s[START] 
						<< ", is type " << System.getPositions()[bond[j].s[START]].type << '\n';
					}
					break;
				case BOND:
					std::cout << "Molecule " << i << " is a bond type.\n";
					std::cout << "\tIt has " << nBonded << " bonds.\n";
					std::cout << "\t\tStarting pairs, " << bond[0].s[0] << " and " << bond[0].s[1] 
					<< ", are types " << System.getPositions()[bond[0].s[0]].type << " and "  
					<< System.getPositions()[bond[0].s[1]].type << '\n';
					break;
				case BEND:
					std::cout << "Molecule " << i << " is a bend type.\n";
					std::cout << "\tIt has " << nBonded << " bends.\n";
					std::cout << "\t\tStarting triplets , " << bond[0].s[0] << " and " << bond[0].s[1] << " and " << bond[0].s[2] 
					<< ", are types " << System.getPositions()[bond[0].s[0]].type << " and "  
					<< System.getPositions()[bond[0].s[1]].type << " and " << System.getPositions()[bond[0].s[2]].type << '\n';
					break;
				default:
					std::cout << "Molecule " << i << " isn't recognized.\n";
					break;
			}
		}
	}
	//Extract molecules and save them new mpd files
	else 
	{
		//read in the rest of the options
		char *modelName=argv[2];
		
		std::stringstream cmdIn;
		for(int i=3;i<argc;i++)
			cmdIn << argv[i] << ' ';
		cmdIn << '\n';
		std::vector<int> molIndices;
		//this is zero if argc==3, so nothing will happen
		for(int i=0;i<argc-3;i++)
		{
			int mol;
			cmdIn >> mol;
			molIndices.push_back(mol);
		}
		
		//set up model structure
		Blob<double> model;
		model.allocParticle(System.readNParticles());
		model.setSize(System.readSize());
		
		//make a bijective map of particles found in molecule
		int *particleMap=new int[System.readNParticles()];
		for(int i=0;i<System.readNParticles();i++)
			particleMap[i]=-1;
		
		//this map is a new set starting with zero
		int newParticleIndex=0;
		
		//efficiency: O(N*log(N))????????
		for(int j=0;j<molIndices.size();j++)
		{
			//shorter names
			int i=molIndices[j];
			molecule< double,fourVector<int> > *m=System.getMolecule();
			int nBonded=m[i].readNBond();
			fourVector<int> *bond=m[i].getBonds();
			
			//our new molecule
			molecule<double,fourVector<int> > newMolecule;
			newMolecule=m[i];
			
			switch(m[i].readType())
			{
				case CHAIN:
					std::cout << "Copying molecule " << i << ", a chain type.\n";
					//This sets up constants, types, and bonds all at once
					
					
					//remap all particles within molecule
					for(int k=0;k<nBonded;k++)
					{
						//range over list
						for(int index=bond[k].s[START];index<bond[k].s[START]+bond[k].s[NCHAINS]*bond[k].s[LENGTH];index++)
						{
							//if it isn't already a value
							if(particleMap[index]==-1)
							{
								particleMap[index]=newParticleIndex++;
								model.addParticle(System.getPositions()[index],System.getVelocities()[index],System.getAccelerations()[index]);
							}
						}
						
						fourVector<int> newBond;
						newBond.s[START]=particleMap[bond[k].s[START]];
						newBond.s[NCHAINS]=bond[k].s[NCHAINS];
						newBond.s[LENGTH]=bond[k].s[LENGTH];
						
						//Actual remapped bond
						//This doesn't work
						newMolecule.setBond(k,newBond);
					}
					
					
					break;
				case BOND:
					std::cout << "Copying molecule " << i << ", a bond type.\n";
					
					//map all particles within molecule
					for(int k=0;k<nBonded;k++)
					{
						//if it isn't already a value
						if(particleMap[bond[k].s[0]]==-1)
						{
							particleMap[bond[k].s[0]]=newParticleIndex++;
							int index=bond[k].s[0];
							model.addParticle(System.getPositions()[index],System.getVelocities()[index],System.getAccelerations()[index]);
						}
						//second particle in bond, check if it isn't already a value
						if(particleMap[bond[k].s[1]]==-1)
						{
							particleMap[bond[k].s[1]]=newParticleIndex++;
							int index=bond[k].s[1];
							model.addParticle(System.getPositions()[index],System.getVelocities()[index],System.getAccelerations()[index]);
						}
						
						fourVector<int> newBond;
						newBond.s[0]=particleMap[bond[k].s[0]];
						newBond.s[1]=particleMap[bond[k].s[1]];
						
						//Actual remapped bond
						newMolecule.setBond(k,newBond);
					}
					break;
				case BEND:
					std::cout << "Copying molecule " << i << ", a bend type.\n";
					//map all particles within molecule
					for(int k=0;k<nBonded;k++)
					{
						//if it isn't already a value
						if(particleMap[bond[k].s[0]]==-1)
						{
							particleMap[bond[k].s[0]]=newParticleIndex++;
							int index=bond[k].s[0];
							model.addParticle(System.getPositions()[index],System.getVelocities()[index],System.getAccelerations()[index]);
						}
						//second particle in bend, check if it isn't already a value
						if(particleMap[bond[k].s[1]]==-1)
						{
							particleMap[bond[k].s[1]]=newParticleIndex++;
							int index=bond[k].s[1];
							model.addParticle(System.getPositions()[index],System.getVelocities()[index],System.getAccelerations()[index]);
						}
						//third particle in bend, check if it isn't already a value
						if(particleMap[bond[k].s[2]]==-1)
						{
							particleMap[bond[k].s[2]]=newParticleIndex++;
							int index=bond[k].s[2];
							model.addParticle(System.getPositions()[index],System.getVelocities()[index],System.getAccelerations()[index]);
						}
						
						fourVector<int> newBond;
						newBond.s[0]=particleMap[bond[k].s[0]];
						newBond.s[1]=particleMap[bond[k].s[1]];
						newBond.s[2]=particleMap[bond[k].s[2]];
						
						//Actual remapped bond
						newMolecule.setBond(k,newBond);
					}
					break;
				default:
					std::cout << "Molecule " << i << " isn't recognized.\n";
					break;
			}
			model.addMolecule(newMolecule);
		}
		
		delete particleMap;
		//write out model's mpd file
		Script<double, Blob <double> > modelIO(modelName,std::ios::out,&model);
		modelIO.write();
		modelIO.close();
	}
	return 0;
}

