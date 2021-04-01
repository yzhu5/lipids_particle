//This is a very short and simple program that shows you how to modify a variable already present in an mpd script.
//The reason this is used is because some variables (forces, potentials, geometry, etc...) are not easily modified
// in the mpd script. Using a file that opens the mpd script, modifies, and then saves it (as below) is far easier
// than modifying the actual variable. Most of it looks like the run.cpp file, it just replaces the "run loop"
// with the modification.

#include "../include/MD.h"
#include "../include/system.h"
#include "../include/algorithms/tesselaSphere.h"


#define CUTOFF 2.0
#define RMIN 1.0

int main(int argc, char* argv[])
{
	if(argc!=14)
	{
		//so simple, this can't possibly mess it up
		std::cout << "Usage: " << argv[0] << " oldName newName vBond cBond kbend particlebond ratio1 ratio2 particledistance vesicleradius particleradius tesselation Umin\n";
		return 0;
	}
	
	
	//read in options to parameters
	char *oldName=argv[1];
	char *newName=argv[2];

	
	threeVector<double> pos1, pos2;
	double vBond, cBond, kbend, ratio1, ratio2, particledistance,vesicleradius, pbond, particleradius,Umin;
	int tesselation;
	std::stringstream cmdIn;
	for(int i=3;i<argc;i++)
		cmdIn << argv[i] << ' ';
	cmdIn >> vBond >> cBond >> kbend >> pbond >> ratio1 >> ratio2 >> particledistance >> vesicleradius >> particleradius >> tesselation >> Umin;
	
	std::cout<< "test1" << std::endl;
	//some lists
	std::vector<position<double>> vertices;
	std::vector<std::vector<fourVector<int>>> edges;
	std::vector<double> lengths; //since the edge have different lengths,creat a list to save them
	std::vector<std::vector<double>> dMap; //adjanct map that save the distance between any two edges
	std::vector<double> cos_O;//different cos_O
	std::vector<std::vector<fourVector<int>>> classfiedcos_O;

	//std::vector<fourVector<int>> faces;
	double Radius;

	std::cout<< "test2" << std::endl;
	//Creat Sphere
	position<double> sphereCenter; //should always set the center at 0, because the relocation does not work right now
    sphereCenter.x = 0;
    sphereCenter.y = 0;
    sphereCenter.z = 0;

    //particle radius
    //particleradius = 19.9766; 

	//creat the first Sphere
	tesselaSphere<double> FirstSphere(sphereCenter, particleradius, tesselation);
	FirstSphere.setcenterpointtype(5);
	FirstSphere.settype1(4);
	FirstSphere.settype2(5);
	FirstSphere.JanusRatioSet(ratio1);

    //creat the second Sphere
    tesselaSphere<double> SecondSphere(sphereCenter, particleradius, tesselation);
    SecondSphere.setcenterpointtype(5);
    SecondSphere.settype1(4);
    SecondSphere.settype2(5);
    SecondSphere.JanusRatioSet(ratio2);


	vertices = FirstSphere.ReturnVertex();
	//Do not forget the center point of the sphere
	vertices.push_back(FirstSphere.ReturnCenterPoint());	
	std::cout<< "test21" << std::endl;
	lengths = FirstSphere.ReturnEdgeLengths();
	edges.resize(lengths.size());
	dMap = FirstSphere.ReturndistanceMap();
		std::cout<< "test22" << std::endl;
	edges = FirstSphere.ReturnClassfiedEdges();
	Radius = FirstSphere.ReturnRadius();
	cos_O = FirstSphere.Returncos_Os();
	classfiedcos_O = FirstSphere.ReturnClassfiedangleIndexList();

	std::cout<< "vertice number: "<<vertices.size()-1 << std::endl;
	
	//Configuration variables
	Blob<double> System;
	
	//load old variables, then initialize them
	Script<double, Blob <double> > fileIO(oldName,std::ios::in,&System);
	fileIO.read();
	fileIO.close();
	
	//the offset for inserting particles, needed to remap molecules from model
	int newParticleOffset=System.readNParticles();
	
	//get the new number of particles
	int newNParticles=System.readNParticles()+vertices.size();
	System.allocParticle(newNParticles);
	
    //get the center of the system
    threeVector<double> center = System.readSize()/2;
    
    std::cout << "system size " << center.x << " " << center.y <<  " " << center.z << std::endl;

    //place the center of two particle
    //the distance between center of vesicle and center of particle
    //tempradius = the radius of vesicle(from center to inner layer of vesicle)
    //             + 5 * bond length(6 beads) + distance between outer layper to the particle(smaller than 1)
    //             + the radius of particle
    double tempradius = vesicleradius + 5*0.7 + 0.8 + particleradius;
    double tempx = sqrt(tempradius*tempradius - (0.5*particledistance)*(0.5*particledistance));
    
    //first particle center 
    pos1.x = center.x + tempx;
    pos1.y = center.y + 0.5*particledistance;
    pos1.z = center.z;

    //second particle center
    pos2.x = center.x + tempx;
    pos2.y = center.y - 0.5*particledistance;
    pos2.z = center.z;


	//Locate center of mass, we are placing at pos
	position<double> comValue=com< position<double> >(&vertices[0], vertices.size());
	
	//Velocities
	MTRand *randNum=new MTRand(System.readSeed());
	double Vrms=sqrt(3.0*System.readInitialTemp());
	
		std::cout<< "test4" << std::endl;
	threeVector<double> zero=0;
	//put vertices in system, repositioning center of mass
	for(int i=0;i<vertices.size();i++)
	{
		vertices[i].x+=pos1.x-comValue.x;
		vertices[i].y+=pos1.y-comValue.y;
		vertices[i].z+=pos1.z-comValue.z;
		threeVector<double> vel;
		double theta=M_PI*randNum->rand53();
		double phi=M_PI*2*randNum->rand53();
		vel.x=Vrms*cos(phi)*sin(theta);
		vel.y=Vrms*sin(phi)*sin(theta);
		vel.z=Vrms*cos(theta);
		
		System.addParticle(vertices[i],vel,zero);
	}


	std::cout<< "test5" << std::endl;



//set constant for interaction between JanusParticle and head
//head = 2; attract part on Janus particle is 4
///initialize constants (force and potential constants)
	//generate force constants
	double *UminS=new double[System.readNTypes()*System.readNTypes()];
	double *UmaxS=new double[System.readNTypes()*System.readNTypes()];
	
	//could be one loop, but it's good for an example
	//of how to set constants
	for(int i=0;i<System.readNTypes();i++)
	{
		for(int j=0;j<System.readNTypes();j++)
		{
			//i=first type, j=second type, indexed grid
			int k=i+j*System.readNTypes();
			UminS[k]=0;
			UmaxS[k]=100;
		}
	}
	
	//umin and umax exceptions, tail types
	UminS[TAIL+TAIL*System.readNTypes()]=-6;
	UmaxS[TAIL+TAIL*System.readNTypes()]=200;

	UminS[HEAD+4*System.readNTypes()] = Umin;
	UmaxS[HEAD+4*System.readNTypes()] = 200;
	UminS[4+HEAD*System.readNTypes()] = UminS[HEAD+4*System.readNTypes()];
	UmaxS[4+HEAD*System.readNTypes()] = UmaxS[HEAD+4*System.readNTypes()];

	UminS[HEAD+5*System.readNTypes()] = 0;
	UmaxS[HEAD+5*System.readNTypes()] = 100;
	UminS[5+HEAD*System.readNTypes()] = UminS[HEAD+5*System.readNTypes()];
	UmaxS[5+HEAD*System.readNTypes()] = UmaxS[HEAD+5*System.readNTypes()];


	
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
			//System.addTwoBodyFconst(RMIN);//C8
			//System.addTwoBodyFconst((2.0*(Umax[g]-Umin[g]))/(RMIN*RMIN));//C6
			//System.addTwoBodyFconst(0);//part of index trick
			//System.addTwoBodyFconst(CUTOFF);//C7
			//System.addTwoBodyFconst((6.0*Umin[g])/((CUTOFF-RMIN)*(CUTOFF-RMIN)));//C3
			//System.addTwoBodyFconst((6.0*Umin[g])/((CUTOFF-RMIN)*(CUTOFF-RMIN)*(CUTOFF-RMIN)));//C2
			System.setTwoBodyFconst(g*6, RMIN);
			System.setTwoBodyFconst(g*6+1, (2.0*(UmaxS[g]-UminS[g]))/(RMIN*RMIN));
			System.setTwoBodyFconst(g*6+2,0);
			System.setTwoBodyFconst(g*6+3, CUTOFF);
			System.setTwoBodyFconst(g*6+4,(6.0*UminS[g])/((CUTOFF-RMIN)*(CUTOFF-RMIN)));
			System.setTwoBodyFconst(g*6+5, (6.0*UminS[g])/((CUTOFF-RMIN)*(CUTOFF-RMIN)*(CUTOFF-RMIN)));



			//U constants, potential constants
			//System.addTwoBodyUconst(RMIN);//C8
			//System.addTwoBodyUconst((Umax[g]-Umin[g])/(RMIN*RMIN));//C4
			//System.addTwoBodyUconst(Umin[g]);//C5,no index trick
			//System.addTwoBodyUconst(CUTOFF);//C7
			//System.addTwoBodyUconst((3.0*Umin[g])/((CUTOFF-RMIN)*(CUTOFF-RMIN)));//C1
			//System.addTwoBodyUconst((2.0*Umin[g])/((CUTOFF-RMIN)*(CUTOFF-RMIN)*(CUTOFF-RMIN)));//C0
			
			System.setTwoBodyUconst(g*6, RMIN);
			System.setTwoBodyUconst(g*6+1, (UmaxS[g]-UminS[g])/(RMIN*RMIN));
			System.setTwoBodyUconst(g*6+2, UminS[g]);
			System.setTwoBodyUconst(g*6+3, CUTOFF);
			System.setTwoBodyUconst(g*6+4, (3.0*UminS[g])/((CUTOFF-RMIN)*(CUTOFF-RMIN)));
			System.setTwoBodyUconst(g*6+5, (2.0*UminS[g])/((CUTOFF-RMIN)*(CUTOFF-RMIN)*(CUTOFF-RMIN)));
		}
	}








	//add the edge bonds into the system for the first system
	for(int i = 0; i < lengths.size();i++)
	{
		molecule<double,fourVector<int>> mBonds;
		mBonds.setType(BOND);
		mBonds.addConstant(lengths[i]);
		mBonds.addConstant(vBond);
		for(int j=0;j<edges[i].size();j++)
	{
		edges[i][j].s[0]+=newParticleOffset;
		edges[i][j].s[1]+=newParticleOffset;
		mBonds.addBond(edges[i][j]);
	}
	
	System.addMolecule(mBonds);
	}


	//add the edge bonds between vertices and center point into the system for the first particle 
	molecule<double,fourVector<int>> cBonds;
	cBonds.setType(BOND);
	cBonds.addConstant(Radius);
	cBonds.addConstant(cBond);
	for(int i = 0; i < vertices.size()-1; i++)
	{
		fourVector<int> tempedge;
		tempedge.s[0] = i + newParticleOffset;
		tempedge.s[1] = vertices.size()-1 + newParticleOffset;
		cBonds.addBond(tempedge);
	}
	System.addMolecule(cBonds);
	
    int firstParticleIndex = vertices.size() - 1 + newParticleOffset;


//////////////////////////////////////////////////////////////////////////////////////

/*
	//add the rigid force between center and vertices
	double Umin=0;
	double Umax=70;
	Umax = cBond;
	double rm=1;
	double rc=2;
	double R=Radius-rm;
	double D=5.88*M_PI*R;
	double D2=pow(2*M_PI*5.88*R, 2);
	double A=Umax-Umin;
*/

//this is not necessary, just for refference
/*
	double C0=rc+R;
	double C1=-7/4*rm;
	double C2=2*rm*rm;
	double C3=Umin*M_PI*R*5.88/(rm*rm*rm);
	double C4=R;
	double C5=rm;
	double C6=-D*A/(2*rm*rm);
	double C7=2*D*A/(3*rm);
	double C8=-D*Umin;
	double C9=2*D*Umin*rm;
	double C10=D*1.3*Umin*rm*rm;
	double C11=2*R+rm;
	double C12=2*R+rc;
	double C13=Umin*D2*2/30;
	double C14=-Umin*D2*7*rm/20;
	double C15=Umin*D2*rm/2;
	double C16=-D2*A*rm/20;
	double C17=D2*A*rm*rm/12;
	double C18=-D2*Umin*pow(rm,3)/6;
	double C19=D2*Umin*pow(rm,4)/2;
	double C20=D2*1.3*Umin*pow(rm,5)/2;
	double C21=D2*13*Umin*pow(rm,6)/60;
*/

/*
	//Make a list of constants
std::vector<double> C;
for(int i=0;i<System.readNTypes();i++)
{
    for(int j=0;j<System.readNTypes();j++)
    {
        C.push_back(rc+R); 
        C.push_back(-7/4*rm);
        C.push_back(2*rm*rm);
        C.push_back(Umin*M_PI*R*5.88/(rm*rm*rm));
        C.push_back(R);
        C.push_back(rm);
        C.push_back(-D*A/(2*rm*rm));
        C.push_back(2*D*A/(3*rm));
        C.push_back(-D*Umin);
        C.push_back(2*D*Umin*rm);
        C.push_back(D*1.3*Umin*rm*rm);
        C.push_back(2*R+rm);
        C.push_back(2*R+rc);
        C.push_back(Umin*D2*2/30);
        C.push_back(-Umin*D2*7*rm/20);
        C.push_back(Umin*D2*rm/2);
        C.push_back(-D2*A*rm/20);
        C.push_back(D2*A*rm*rm/12);
        C.push_back(-D2*Umin*pow(rm,3)/6);
        C.push_back(D2*Umin*pow(rm,4)/2);
        C.push_back(D2*1.3*Umin*pow(rm,5)/2);
        C.push_back(D2*13*Umin*pow(rm,6)/60);
    }
}

*/

/*

		//create and add the rigid interaction
		molecule<double,fourVector<int>> m;
		m.setType(9);//or BEAD

		//add constants from our constants list
		for(int i=0;i<C.size();i++)
			m.addConstant(C[i]);
		fourVector<int> bond;
		bond.s[0]=System.readNParticles()-1;
		m.addBond(bond);
		System.addMolecule(m);

*/
//////////////////////////////////////////////////////////////////////////



	//add the bend force for the first particle
	for(int i = 0; i < cos_O.size(); i++)
	{
		molecule<double,fourVector<int>> Bends;
		Bends.setType(BEND);
		Bends.addConstant(-cos_O[i]);
		Bends.addConstant(kbend);
		for(int j = 0; j< classfiedcos_O[i].size(); j++)
		{
			classfiedcos_O[i][j].s[0] += newParticleOffset;
			classfiedcos_O[i][j].s[1] += newParticleOffset;
			classfiedcos_O[i][j].s[2] += newParticleOffset;
			Bends.addBond(classfiedcos_O[i][j]);
		}
		System.addMolecule(Bends);
	}
	


	
    ////////////////////////////////////////////////////////////////////////////////
    //locate the second particle

    //locate the important parameter 
    vertices.clear();
    vertices = SecondSphere.ReturnVertex();
	//Do not forget the center point of the sphere
	vertices.push_back(SecondSphere.ReturnCenterPoint());	
	std::cout<< "test21" << std::endl;
    lengths.clear();
	lengths = SecondSphere.ReturnEdgeLengths();
    edges.clear();
	edges.resize(lengths.size());
	dMap = SecondSphere.ReturndistanceMap();
	std::cout<< "test22" << std::endl;
	edges = SecondSphere.ReturnClassfiedEdges();
	Radius = SecondSphere.ReturnRadius();
    cos_O.clear();
	cos_O = SecondSphere.Returncos_Os();
    classfiedcos_O.clear();
	classfiedcos_O = SecondSphere.ReturnClassfiedangleIndexList();
	std::cout<< "test3" << std::endl;



    //redo this two step for second particle
    //the offset for inserting second particle, needed to remap molecules from model
	newParticleOffset=System.readNParticles();
	
	//get the new number of particles
	newNParticles=System.readNParticles()+vertices.size();
	System.allocParticle(newNParticles);



    //Locate center of mass, we are placing at pos
	comValue=com< position<double> >(&vertices[0], vertices.size());
	
	
	std::cout<< "test4" << std::endl;
	//put vertices of second particle in system, repositioning center of mass
	for(int i=0;i<vertices.size();i++)
	{
		vertices[i].x+=pos2.x-comValue.x;
		vertices[i].y+=pos2.y-comValue.y;
		vertices[i].z+=pos2.z-comValue.z;
		threeVector<double> vel;
		double theta=M_PI*randNum->rand53();
		double phi=M_PI*2*randNum->rand53();
		vel.x=Vrms*cos(phi)*sin(theta);
		vel.y=Vrms*sin(phi)*sin(theta);
		vel.z=Vrms*cos(theta);
		
		System.addParticle(vertices[i],vel,zero);
	}


	std::cout<< "test5" << std::endl;


//add the edge bonds(betweem verticle) of second particle into the system
	for(int i = 0; i < lengths.size();i++)
	{
		molecule<double,fourVector<int>> mBonds;
		mBonds.setType(BOND);
		mBonds.addConstant(lengths[i]);
		mBonds.addConstant(vBond);
		for(int j=0;j<edges[i].size();j++)
	{
		edges[i][j].s[0]+=newParticleOffset;
		edges[i][j].s[1]+=newParticleOffset;
		mBonds.addBond(edges[i][j]);
	}
	
	System.addMolecule(mBonds);
	}


	//add the edge bonds between vertices and center point into the system for the second particle 
	molecule<double,fourVector<int>> cBonds2;
	cBonds2.setType(BOND);
	cBonds2.addConstant(Radius);
	cBonds2.addConstant(cBond);
	for(int i = 0; i < vertices.size()-1; i++)
	{
		fourVector<int> tempedge;
		tempedge.s[0] = i + newParticleOffset;
		tempedge.s[1] = vertices.size()-1 + newParticleOffset;
		cBonds2.addBond(tempedge);
	}
	System.addMolecule(cBonds2);
	
    int secondParticleIndex = vertices.size()-1 + newParticleOffset;



	//add the bend force for the second particle
	for(int i = 0; i < cos_O.size(); i++)
	{
		molecule<double,fourVector<int>> Bends;
		Bends.setType(BEND);
		Bends.addConstant(-cos_O[i]);
		Bends.addConstant(kbend);
		for(int j = 0; j< classfiedcos_O[i].size(); j++)
		{
			classfiedcos_O[i][j].s[0] += newParticleOffset;
			classfiedcos_O[i][j].s[1] += newParticleOffset;
			classfiedcos_O[i][j].s[2] += newParticleOffset;
			Bends.addBond(classfiedcos_O[i][j]);
		}
		System.addMolecule(Bends);
	}
	



    //add the Elastic force between two particle center
	molecule<double,fourVector<int>> particlebond;
    particlebond.setType(BOND);
    particlebond.addConstant(particledistance);
    particlebond.addConstant(pbond);
    fourVector<int> twoparticleedge;
    twoparticleedge.s[0] = firstParticleIndex;
    twoparticleedge.s[1] = secondParticleIndex;
    particlebond.addBond(twoparticleedge);
    System.addMolecule(particlebond);






    //this comment just for format refferece, do not un-comment 
	/*
	
	//add the edge bonds into the system
	molecule<double,fourVector<int>> mBonds;
	mBonds.setType(BOND);
	mBonds.addConstant(lBond);
	mBonds.addConstant(kBond);
	//for(int i=0;i<edges.size();i++)
	for(int i=0;i<1920;i++)
	{
		edges[i].s[0]+=newParticleOffset;
		edges[i].s[1]+=newParticleOffset;
		mBonds.addBond(edges[i]);
	}
	
	System.addMolecule(mBonds);
	
	
	//add the center edge bonds into the system
	molecule<double,fourVector<int>> cBonds;
	cBonds.setType(BOND);
	cBonds.addConstant(10);
	cBonds.addConstant(10);
	for(int i=1920;i<edges.size();i++)
	{
		edges[i].s[0]+=newParticleOffset;
		edges[i].s[1]+=newParticleOffset;
		cBonds.addBond(edges[i]);
	}
	
	System.addMolecule(cBonds);
	
	*/
	//molecule<double,fourVector<int>> mBends;
	//mBends.setType(BEND);
	//mBends.addConstant(cosThetaBend);
	//mBends.addConstant(kBend);
	//for(int i=0;i<faces.size();i++)
	//{
	//	faces[i].s[0]+=newParticleOffset;
	//	faces[i].s[1]+=newParticleOffset;
	//	faces[i].s[2]+=newParticleOffset;
	//	mBends.addBond(faces[i]);
	//}
	
	//Save system
	fileIO.open(newName,std::ios::out);
	fileIO.write();
	fileIO.close();
	
	return 0;
}

