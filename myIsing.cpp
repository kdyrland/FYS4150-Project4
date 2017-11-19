#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <armadillo>
#include <string>

using namespace std;
using namespace arma;

ofstream ofile;

inline int PeriodicBoundary(int i, int limit, int add){
	return (i+limit+add) % (limit);
}

//Initialize the lattize matrix (energy and magnetization)
void InitLattice(int, mat &, double&, double&);
//Metropolis Algorithm
void Metropolis(int, int, double, vec&);
//Write the results
void WritetoFile(int, int, double,vec);

//-------------------------------------------------------------------------------------


int main(int argc, char* argv[]){

	// Initialize the arguments defined by the program
	string filename;
	int NSpin, MCcycles;
	double Tempi, Tempf, Tempstep;
	if (argc <= 5) {
		cout << "Bad Usage:" << argv[0] << "read output file, number of spins, MC cycles, intial temp, final temp, temp step." << endl;
		exit(1);
	}

	// Saving arguments
	if (argc > 1) {
		filename=argv[1];
		NSpin=atoi(argv[2]);	// number of spins - 2 for 4b
		MCcycles=atoi(argv[3]);	// monte carlo cycles - 100...
		Tempi=atof(argv[4]);	// initial temperature - 1...
		Tempf=atof(argv[5]);	// final temperature - 1...
		Tempstep=atof(argv[6]);	// temp steps - 1...
	}


	// Create a new file name with the lattice size added, then number of MCcycles
	string fileout = filename;
	string argument = to_string(NSpin);
	string argument2 = to_string(MCcycles);
	fileout.append(argument);
	fileout += "_";
	fileout.append(argument2);
	fileout += ".txt";
	ofile.open(fileout);


	// Print arguments to the terminal
	cout << "File name:" << fileout << endl;
	cout << "Number of spins:" << NSpin << endl;
	cout << "Number of MC cycles:" << MCcycles << endl;
	cout << "Initial Temperature:" << Tempi << endl;
	cout << "Final Temperature:" << Tempf << endl;
	cout << "Temperature Step Size:" << Tempstep << endl;

	// Temperature Loop
	for (double Temperature = Tempi; Temperature <= Tempf; Temperature+=Tempstep){
		vec ExpectationValues = zeros<mat>(5);
		// Time to Monte the Carlo
		Metropolis(NSpin,MCcycles,Temperature,ExpectationValues);
		// Then write the expectation values found to the file
		WritetoFile(NSpin,MCcycles,Temperature,ExpectationValues);
	}
	ofile.close();
 	return 0;
 }

// The Metropolis algorithm to filp a spin and find the new energy and magnetization of the system
void Metropolis(int NSpin, int MCcycles, double Temperature, vec &ExpectationValues){
	
	// Preparing to use the Mersenne randomizer algorithm
	std::random_device rd;
	std::mt19937_64 gen(rd());

	// Using a uniform distribution for x in [0,1]
	std::uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);

	// Initialize lattice spin values
	mat SpinMat = zeros<mat>(NSpin,NSpin);

	// Initialize energy and magnetization
	double energy = 0; double Magnet = 0;

	// Initialize expectation values
	InitLattice(NSpin,SpinMat,energy,Magnet);

	// Possible energy differences are precalculated to save computation time
	vec EnergyDiff = zeros<mat>(17);
	for (int de=-8; de<=8; de+=4) EnergyDiff(de+8) = exp(-de/Temperature);
	// EnergyDiff = [e^-8/Ti 0 0 0 e^-4/Ti 0 0 0 0 0 0 0 e^4/Ti 0 0 0 e^8/Ti]
	
	//Bring in Monte!
	int count = 0;
	for (int cycles=1; cycles <= MCcycles; cycles++){
		for (int x=0;x<NSpin; x++){
			for (int y=0;y<NSpin;y++){
				// pick a # bw 0 and one and turn into an integer 0 or 1:
				int ix = (int) (RandomNumberGenerator(gen) * (double)NSpin);
				int iy = (int) (RandomNumberGenerator(gen) * (double)NSpin);
				int deltaE = 2 * SpinMat(ix,iy) * (			// dE = 2s_k*sum(s_l)
					SpinMat(ix,PeriodicBoundary(iy,NSpin,-1)) +
					SpinMat(PeriodicBoundary(ix,NSpin,-1),iy) +
					SpinMat(ix,PeriodicBoundary(iy,NSpin,1)) +
					SpinMat(PeriodicBoundary(ix,NSpin,1),iy));
				if (RandomNumberGenerator(gen) <= EnergyDiff(deltaE+8)){
					SpinMat(ix,iy) *= -1.0; 				//flip one spin and accept
					Magnet += (double) 2 * SpinMat(ix,iy); 	//2 bc one for initial spin and one for final
					energy += (double) deltaE;
					count += 1;
				}
			}
		}

		// Keeping all of the expectation values
		ExpectationValues(0) += energy;
		ExpectationValues(1) += energy*energy;
		ExpectationValues(2) += Magnet;
		ExpectationValues(3) += Magnet*Magnet;
		ExpectationValues(4) += fabs(Magnet);
	}
cout << "Accepted configs: " << count << endl;
}


// Initiallize the Ising spin lattice (Matrix)
void InitLattice(int NSpin, mat &SpinMat, double &energy, double &Magnet){
	for (int x=0;x<NSpin; x++){
		for (int y=0;y<NSpin;y++){
			std::random_device rd;
			std::mt19937_64 gen(rd());
			std::uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);
			if (RandomNumberGenerator(gen) <= 0.5) {
				SpinMat(x,y) = 1.0;
			}
			else {
				SpinMat(x,y) = -1.0;
			}
			//SpinMat(x,y)=1.0; //Everything in the ground state = ordered up
			Magnet += (double)SpinMat(x,y);
		}
	}
	for (int x=0;x<NSpin; x++){
		for (int y=0;y<NSpin;y++){
			energy -= (double)SpinMat(x,y)*
			(SpinMat(PeriodicBoundary(x,NSpin,-1),y)+
			 SpinMat(x,PeriodicBoundary(y,NSpin,-1)));
		}
	}
}

// Program to write the results of the Ising Model to a file
void WritetoFile(int NSpin, int MCcycles, double Temperature, vec ExpectationValues){
	double norm = 1.0/((double)(MCcycles)); 				//Normalizing factor (1/N)
	double expectationE = ExpectationValues(0)*norm;
	double expectationE2 = ExpectationValues(1)*norm;
	double expectationM = ExpectationValues(2)*norm;
	double expectationM2 = ExpectationValues(3)*norm;
	double expectationMabs = ExpectationValues(4)*norm;
	// but the true expectation values is per spins...
	double Evariance = (expectationE2 - expectationE*expectationE);
	double Mvariance = (expectationM2 - expectationMabs*expectationMabs);
	ofile << setiosflags(ios::showpoint | ios::uppercase);
	ofile << setw(15) << setprecision(8) << Temperature;	// temp
	ofile << setw(15) << setprecision(8) << expectationE;	// energy
	ofile << setw(15) << setprecision(8) << Evariance/Temperature/Temperature; // specific heat
	ofile << setw(15) << setprecision(8) << expectationM;	// magnetization
	ofile << setw(15) << setprecision(8) << Mvariance/Temperature; 		// susceptibility
	ofile << setw(15) << setprecision(8) << expectationMabs << endl; 	// abs magnetization
}//end of the output function