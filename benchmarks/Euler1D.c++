/*
 * Euler1D.c++
 *
 *  Created on: Aug 28, 2009
 *      Author: meithan
 */

#import <math.h>
#import <string>

// GLOBAL VARIABLES (EVIL)

// Simulation parameters
const int NX = 5000;		// Number of grid points
const double L = 1.0;		// Grid physical length (arbitrary units)
const double NT = 0.5;	// Final integration time (arbitrary units)
const double CP = 0.5;	// Courant Parameter
const double GAM = 1.4;	// Heat capacities ratio
const int NOUT = 20;		// Number of outputs

// Benchmarking flags
// quiet: non-verbose progress report
// bench: no data outputs + quiet flag
const bool quiet = true;
const bool bench = true;

// Declare main physical variables
double U[NX][3];
double UP[NX][3];
double P[NX][3];
double F[NX][3];

// Other simulation variables
int IT = 0;
double T = 0;
double DT = 0;
double DX = L / NX;
int NDUMP = 0;

//Calculate integration variables from primitives
void getvars(){
	for (int X=0; X<NX; X++){
		U[X][0] = P[X][1];
		U[X][1] = P[X][1]*P[X][0];
		U[X][2] = 0.5*P[X][1]*P[X][0]*P[X][0] + P[X][2]/(GAM-1);
	}
}

// Load initial conditions
void initcond(){
	// Sod Shock Tube initial conditions
	// This is a standard test for hydro codes
	// Velocity = 0 everywhere
	// Density = 1.0 on left half, 0.125 on right half
	// Pressure = 1.0 on left, 0.1 on right
	for (int X=0; X<NX; X++){
		if (X+1<=NX/2){
			P[X][0] = 0.0;
			P[X][1] = 1.0;
			P[X][2] = 1.0;
		} else {
			P[X][0] = 0.0;
			P[X][1] = 0.125;
			P[X][2] = 0.1;
		}
	}
	getvars();
}

//  Calculate primitives from integration variables
void getprims(){
	for (int X=0; X<NX; X++){
		P[X][0] = U[X][1]/U[X][0];
		P[X][1] = U[X][0];
		P[X][2] = (GAM-1)*(U[X][2]-0.5*U[X][1]*U[X][1]/U[X][0]);
	}
}

// Boundary Conditions
void boundary(){
	// Reflection boundary condition
	UP[0][0] = UP[1][0];
	UP[NX-1][0] = UP[NX-2][0];
	UP[0][1] = (-1.0)*UP[1][1];
	UP[NX-1][1] = (-1.0)*UP[NX-2][1];
	UP[0][2] = UP[1][2];
	UP[NX-1][2] = UP[NX-2][2];
}


// Calculate Timestep
double getstep(){
	// The timestep is calculated following the
	// Courant-Friedrichs-Lewis criterion:
	// DT = DX / UMAX * CP
	// where DX is the grid spacing, CP is the
	// Courant parameter and UMAX is the maximum
	// value of local velocity + local sound speed
	double UMAX = 0;
	double ULOC = 0;
	for (int X=0; X<NX; X++){
		ULOC = fabs(P[X][0]) + sqrt(GAM*P[X][2]/P[X][1]);
		if (ULOC > UMAX){
			UMAX = ULOC;
		}
	}
	return DX / UMAX * CP;
}

// Lax Method
void lax(){
	// Calculate fluxes
	for (int X=0; X<NX; X++){
		F[X][0] = P[X][1]*P[X][0];
		F[X][1] = P[X][1]*P[X][0]*P[X][0] + P[X][2];
		F[X][2] = P[X][0]*(0.5*P[X][1]*P[X][0]*P[X][0] + GAM/(GAM-1)*P[X][2]);
	}
	// Calculate timestepped variables for non-boundary points
	for (int X=1; X<NX-1; X++){
		for (int V=0; V<3; V++){
			UP[X][V] =  0.5*(U[X+1][V] + U[X-1][V]) - 0.5 * DT/DX * (F[X+1][V] - F[X-1][V]);
		}
	}
	// Treats boundary grid points
	boundary();
}


// Advance Simulation
void advance(){
	for (int X=0; X<NX; X++){
		for (int V=0; V<3; V++){
			U[X][V] = UP[X][V];
		}
	}
	T = T + DT;
	IT = IT + 1;
}

// Dump primitives to output data file
static void output(){
//		// Build output filename
//		String number = String.format("%04d", NDUMP);
//		String fname = "./data/Sod." + number + ".dat";
//		// Open output file
//		FileWriter fstream = new FileWriter(fname);
//		BufferedWriter datafile = new BufferedWriter(fstream);
//		// Report dump, if not quiet
//		if (!quiet) System.out.println("Dumping Output " +  NDUMP + " to " + fname);
//		// Write Header
//		datafile.write("# Dumped on XXX\n");
//		datafile.write("# Using method XXX\n");
//		datafile.write("# IT=" + IT + ", T=" + T + " (" + T/NT + "%)\n");
//		// Write Data
//		for (int X=0; X<NX; X++){
//			String buffer = String.format("%5f\t%5f\t%5f\t%5f\n",X*DX,P[X][0],P[X][1],P[X][2]);
//			datafile.write(buffer);
//		}
//		// Close file
//		datafile.close();
//
}

// Report progress and dump output if needed
void progress(){

	if (!quiet) {
		printf("IT=%5d, T=%7.4f (%7.4f%%)\n",IT,T,T/NT*100);
	}

	if (T/NT >= (float)NDUMP/NOUT) {
		printf("Progress: %7.4f%% (IT=%5d)\n",T/NT*100,IT);
		if (!bench) output();
		NDUMP += 1;
	}
}


int main(){

	// Run the program
	initcond();
	output();

	// Main Loop
	while(T<=NT){
		DT = getstep();  // Get timestep
		lax();         	 // Lax Solver
		advance();     	 // Advance simulation
		getprims();    	 // Update primitives
		progress();   	 // Show progress; dump output if needed
	}

	return 0;

}
