// C++ headers
#include <algorithm>
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdio.h>


// My headers
#include "header.hh"

// Number of points in N x N domain
int const Nx = 100;
int const Ny = 100;

double const L = 1;
double const H = 1;

double const dx = L / Nx;
double const dy = H / Ny;

double const alpha = 1;

double const errorTolerance = 0.0001;

double LNorm = 0.0;

// Create cell array
cell cells[Nx+1][Ny+1];

int main() { // Start of program
    /**
    Test a comment here.
    This is where the program starts
    */
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

	double l1 = calcLNorm( 1.0 );
	double l2 = calcLNorm( 2.0 );

	std::ofstream file( fileName + ".txt", std::ofstream::out );

	file << "L1 Norm: " << l1 << std::endl;
	file << "L2 Norm: " << l2 << std::endl;

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

	int counter = 0;

	double dt = 0.00001;

	for ( double t = 0; t < 1; t = t + dt )
	{
		fieldUpdateExplicit( dt );
		shiftTempsForNewTimeStep();
	}
}

void fieldUpdateExplicit(
	double dt
)
{
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
			cells[i][j].temp_iter = cells[i][j].temp;
		}
	}
}

void simulateImplicit() {

	double dt = 0.00001;

	for ( double t = 0; t < 1; t = t + dt ) {

		bool isConverged = false;

		while( !isConverged ) { // Iteration loop
			shiftTempsForNewIteration();
			fieldUpdateImplicit( dt );
			isConverged = isConverged_iter();
		}

		shiftTempsForNewTimeStep();

	}
}

void fieldUpdateImplicit(
	double dt
)
{
	for ( int j = 1; j < Ny; ++j ) { // looping from i/j = 1 to N-1
		for ( int i = 1; i < Nx; ++i ) {
			cells[i][j].temp = alpha * dt * ( ( cells[i+1][j].temp - 2 * cells[i][j].temp + cells[i-1][j].temp ) / std::pow( dx, 2.0 ) + \
												( cells[i][j+1].temp - 2 * cells[i][j].temp + cells[i][j-1].temp ) / std::pow( dy, 2.0 ) ) + cells[i][j].temp_prev_ts;
		}
	}
}

double calcLNorm(
	double normPower
)
{
	double sumError = 0.0;

	for ( int j = 1; j < Ny; ++j ) {
		for ( int i = 1; i < Nx; ++i ) {
			auto & c( cells[i][j] );
			double err = (c.temp - ( c.x_loc + c.y_loc ) );
			sumError += std::abs( std::pow( err, normPower ) );
		}
	}

	return std::pow( sumError, 1 / normPower );
}
