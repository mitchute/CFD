// C++ headers
#include <iostream>
#include <fstream>
#include <algorithm>

// My headers
#include "header.hh"

// Number of points in N x N domain
int const Nx = 80;
int const Ny = 80;

double const L = 1;
double const H = 1;

double const dx = L / Nx;
double const dy = H / Ny;

double const alpha = 1;

double const errorTolerance = 0.0001;

int counter = 0;
double maxError = 0.0;

// Create cell array
cell cells[Nx+1][Ny+1];

int main() { // Start of program
	 explicitMethod();
	 implicitMethod();
}

void initDomain() { // Initialize domain correctly for each run.

	for ( int j = 0; j <= Ny; ++j ) {
		double y = double(j) * dy;
		for ( int i = 0; i <= Nx; ++i ) {
			double x = double(i) * dx;
			cells[i][j].x_index = i;
			cells[i][j].y_index = j;
			cells[i][j].x_loc = x;
			cells[i][j].y_loc = y;

			if ( i == 0 || i == Nx || j == 0 || j == Ny ) { // All edges set to x + y = u(x,y)
				cells[i][j].temp = x + y;
				cells[i][j].temp_iter = x + y;
				cells[i][j].temp_prev_ts = x + y;
			} else { // Everywhere else initialized to 0.0
				cells[i][j].temp = 0.0;
				cells[i][j].temp_iter = 100; // initialized big
				cells[i][j].temp_prev_ts = 0.0;
				//cells[i][j].temp = x + y; // for testing the contour plot only
			}
		}
	}
}

void reportResults( std::string const fileName ) { // Report results to file

	std::ofstream file( fileName + ".csv", std::ofstream::out );

	for ( int j = 0; j <= Ny; ++j ) {
		if ( j == 0 ) { // print x locations
			file << ","; // pad space to get correct alignment
			for ( int i = 0; i <= Nx; ++i ) {
				if ( i == Nx ) {
					file << cells[i][j].x_loc << std::endl;
				} else {
					file << cells[i][j].x_loc << ",";
				}
			}
		}

		for ( int i = 0; i <= Nx; ++i ) {
			if ( i == 0 ) { // print y locations
				file << cells[i][j].y_loc << ",";
			}

			if ( i == Nx ) {
				file << cells[i][j].temp;
			} else {
				file << cells[i][j].temp << ",";
			}
		}
		file << std::endl;
	}

	file.close();
}

void reportStatistics( std::string const fileName ) {	// Report statistics from run

	calcMaxError();

	std::ofstream file( fileName + ".txt", std::ofstream::out );

	file << "Num. Iterations: " << counter << std::endl;
	file << "Abs. Error: " << maxError << std::endl;

	file.close();
}

void explicitMethod() { // Simulate explicit method
	initDomain();
	simulateExplicit();
	reportResults( "explicit" );
	reportStatistics( "explicit" );
}

void implicitMethod() { // Simulate implicit method
	initDomain();
	simulateImplicit();
	reportResults( "implicit" );
	reportStatistics( "implicit" );
}

void simulateExplicit() {

	bool converged = false;
	int counter = 0;

	while ( !converged )
	{
		fieldUpdateExplicit();
		converged = isConverged_ts();
		shiftTempsForNewTimeStep();
	}
}

void fieldUpdateExplicit() {

	double dt = std::pow( L / Nx, 2.0 ) * 0.25 / alpha;

	for ( int j = 1; j < Ny; ++j ) { // looping from i/j = 1 to N-1
		for ( int i = 1; i < Nx; ++i ) {
			cells[i][j].temp = alpha * dt * ( ( cells[i+1][j].temp_prev_ts - 2 * cells[i][j].temp_prev_ts + cells[i-1][j].temp_prev_ts ) / std::pow( dx, 2.0 ) \
											+ ( cells[i][j+1].temp_prev_ts - 2 * cells[i][j].temp_prev_ts + cells[i][j-1].temp_prev_ts ) / std::pow( dy, 2.0 ) ) + cells[i][j].temp_prev_ts;
		}
	}
}

bool isConverged_ts() {

	bool converged = true;

	for ( int j = 1; j < Ny; ++j ) { // looping from i/j = 1 to N-1
		for ( int i = 1; i < Nx; ++i ) {
			auto & thisCell( cells[i][j]);
			if ( abs( thisCell.temp - thisCell.temp_prev_ts ) > errorTolerance ) {
				converged = false;
				return converged;
			}
		}
	}

	return converged;
}

bool isConverged_iter() {

	bool converged = true;

	for ( int j = 1; j < Ny; ++j ) { // looping from i/j = 1 to N-1
		for ( int i = 1; i < Nx; ++i ) {
			auto & thisCell( cells[i][j]);
			if ( abs( thisCell.temp - thisCell.temp_iter ) > errorTolerance ) {
				converged = false;
				return converged;
			}
		}
	}

	return converged;
}

void shiftTempsForNewTimeStep() {
	for ( int j = 1; j < Ny; ++j ) { // looping from i/j = 1 to N-1
		for ( int i = 1; i < Nx; ++i ) {
			cells[i][j].temp_prev_ts = cells[i][j].temp;
		}
	}
}

void shiftTempsForNewIteration() {
	for ( int j = 1; j < Ny; ++j ) { // looping from i/j = 1 to N-1
		for ( int i = 1; i < Nx; ++i ) {
			cells[i][j].temp = cells[i][j].temp_iter;
		}
	}
}

void simulateImplicit() {

	bool isConverged = false;

	while( !isConverged ) { // Iteration loop
		++counter;
		fieldUpdateImplicit();
		//shiftTempsForNewIteration();
		isConverged = isConverged_ts();
		shiftTempsForNewTimeStep();
	}
}

void fieldUpdateImplicit() {

	double dt = std::pow( L / Nx, 2.0 ) * 0.25 / alpha;

	for ( int j = 1; j < Ny; ++j ) { // looping from i/j = 1 to N-1
		for ( int i = 1; i < Nx; ++i ) {
			cells[i][j].temp = alpha * dt * ( ( cells[i+1][j].temp - 2 * cells[i][j].temp + cells[i-1][j].temp ) / std::pow( dx, 2.0 ) + \
												( cells[i][j+1].temp - 2 * cells[i][j].temp + cells[i][j-1].temp ) / std::pow( dy, 2.0 ) ) + cells[i][j].temp_prev_ts;
		}
	}
}

void calcMaxError() {

	for ( int j = 1; j < Ny; ++j ) { // looping from i/j = 1 to N-1
		for ( int i = 1; i < Nx; ++i ) {
			auto & thisCell( cells[i][j] );
			maxError = std::max( maxError, std::abs( thisCell.temp - ( thisCell.x_loc + thisCell.y_loc ) ) );
		}
	}
}
