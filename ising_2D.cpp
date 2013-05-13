#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <time.h>

using namespace std;

//File to throw data into
//ofstream DATA("output.dat",ios::out);

//struct for 2D lattice
struct lat_coord{
	int x;
	int y;
};

const int size	 = 1000;		//lattice size
const int lsize	 = size-1;		//array size for lattice
const int n	 = size*size;		//number of spin points on sq latice
float T		 = 0.001;		//starting point for temperature
const float minT = 0.001;		//minimum temperature
float change	 = 0.01;		//size of steps of temp loop
int lat[size+1][size+1];		//2D lattice for holding spins
long unsigned int mcs = 100000;		//number of Monte Carlo steps
int transient	 = 10000;		//number of transient steps
double norm	 = (1.0/float(mcs*n));  //normalization const for avgs

//initializer
void initialize(){
	srand(time(NULL));		//random seed
	for(int i=size;i>=1;i--){
		for(int j=1;j<=size;j++){
			if((rand()%100)/100.>=0.5)
				lat[i][j]=1;
			else
				lat[i][j]=-1;
		}
	}
}

//print out lattice on screen for funsies
void output(){
	for(int i=size; i>=1;i--){
		for(int j=1;j<=size;j++){
			if(lat[i][j]<0)
				cout << " - ";
			else
				cout << " + ";
		}
		cout << endl;
	}
}

//function for choosing random position on lattice
void choose_rand_pos_lat(lat_coord &pos){
	pos.x = (int)ceil(rand()%size);
	pos.y = (int)ceil(rand()%size);
	if(pos.x>size || pos.y>size){
		cout << "Something went wrong";
		exit;
	}
}

//Calculate energy at particular position on the lattice
int energy_coord(lat_coord &pos){
	//P.B.C. 
	int up,down,left,right,E;
	if(pos.y==size)
		up=1;
	else
		up=pos.y+1;
	if(pos.y==1)
		down=size;
	else
		down=pos.y-1;
	if(pos.x==1)
		left=size;
	else
		left=pos.x-1;
	if(pos.x==size)
		right=1;
	else
		right=pos.x+1;

	//energy
	E=-1*lat[pos.x][pos.y]*(lat[left][pos.y]+lat[right][pos.y]+lat[pos.x][up]+lat[pos.x][down]);
	return E;
}

//function for testing the validity of flipping a spin at a position
bool test_flip(lat_coord &pos, int &dE){
	dE=-2*energy_coord(pos);	//change in energy for specific spin
	if(dE<0)			//flip due to lower energy
		return true;
	else if(rand()%100/100.<exp(-dE/T))	//flip due to heat bath
		return true;
	else
		return false;		//no flip :/
}

//flip spin at given position
void flipIt(lat_coord pos){
	lat[pos.x][pos.y] = -lat[pos.x][pos.y];

}

//function for disregarding the transient results
void transient_result(){
	lat_coord pos;
	int dE = 0;
	for(int i=1; i<=transient; i++){
		for(int j=1; j<=n; j++){
			choose_rand_pos_lat(pos);
			if(test_flip(pos,dE))
				flipIt(pos);
		}
	}
}

//figure out the total magnetization of lattice
int total_magnetization(){
	int m = 0;
	for(int i=size; i>=1; i--){
		for(int j=1; j<=size; j++){
			m+=lat[i][j];
		}
	}
	return m;
}

//figure out the total energy lattice
int total_energy(){
	lat_coord pos;
	int E = 0;
	for(int y=size; y>=1; y--){
		pos.y=y;
		for(int x=1; x<=size; x++){
			pos.x=x;
			E+=energy_coord(pos);
		}
	}
	return E;
}

int main(){
	cout << "HI" << endl;
	//set up file for data
	ofstream myfile;
	myfile.open("Ising2D_Data.dat");
	//declare variables
	double E=0, Esq=0, Esq_avg=0, E_avg=0, etot=0, etotsq=0;
	double M=0, Msq=0, Msq_avg=0, M_avg=0, mtot=0, mtotsq=0;
	double Mabs=0, Mabs_avg=0, Mq_avg=0, mabstot=0, mqtot=0;
	int de = 0;
	lat_coord pos;

	//initialize lattice with random config
	initialize();
	//start a temperature loop
	for(;T>=minT;T=T-change){
		//transient phase
		transient_result();
		//observables adopt equilibrated lattice config values
		M=total_magnetization();
		Mabs-abs(total_magnetization());
		E=total_energy();
		//initialize summation variables at each temp step
		etot=0;
		etotsq=0;
		mtot=0;
		mtotsq=0;
		mabstot=0;
		mqtot=0;

		//Monte Carlo loop
		for(int i=1; i<mcs; i++){
			//Metropolis loop
			for(int j=1; j<=n; j++){
				choose_rand_pos_lat(pos);
				if(test_flip(pos, de)){
					flipIt(pos);
					//adjust observables
					E+=2*de;
					M+=2*lat[pos.x][pos.y];
					Mabs+=abs(lat[pos.x][pos.y]);
				}
			}

		//keep summation over observables
		etot+=E/2.0;	//don't double count spins
		etotsq+=E/2.0*E/2.0;
		mtot+=M;
		mtotsq+=M*M;
		mqtot+=M*M*M*M;
		mabstot+=sqrt(M*M);
		}
	//average out and normalize the observables
	E_avg=etot*norm;
	Esq_avg=etotsq*norm;
	M_avg=mtot*norm;
	Msq_avg=mtotsq*norm;
	Mabs_avg=mabstot*norm;
	Mq_avg=mqtot*norm;

	//Output the data to a file
	myfile << T <<		 				//temperature
	"\t" << M_avg << "\t" << Mabs_avg << "\t" << Msq_avg << //<M>; <|M|>;<M^2> per
	"\t" << (Msq_avg-(M_avg*M_avg*n))/T <<		//susceptibility per spin
	"\t" << (Msq_avg-(Mabs_avg*Mabs_avg*n))/T <<	//susceptibility avg
	"\t" << E_avg << "\t" << Esq_avg <<		//<E>;<E^2> per spin
	"\t" << (Esq_avg-(E_avg*E_avg*n))/(T*T) <<	//heat capacity per spin
	"\t" << 1-((Mq_avg)/(3*Msq_avg))<<endl;		//cumulant (U_L)

	}
	//Close file b/c we are nice ppl
	myfile.close();
	return 0;
}
