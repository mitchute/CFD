// C++ headers
#include <iostream>
#include <fstream>

// My headers
#include "header.hh"

// Number of points in N x N domain
int const Nx = 20;
int const Ny = 20;

double const L = 1;
double const H = 1;

double const dx = L / Nx;
double const dy = H / Ny;

double const alpha = 1;

double const errorTolerance = 0.0001;

// Create cell array
cell cells[Nx+1][Ny+1];

int main()
{

	 explicitMethod();

	 implicitMethod();

}

void initDomain() {

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
				cells[i][j].temp_prev_iter = x + y;
				cells[i][j].temp_prev_ts = x + y;
			} else { // Everywhere else initialized to 0.0
				cells[i][j].temp = 0.0;
				cells[i][j].temp_prev_iter = 0.0;
				cells[i][j].temp_prev_ts = 0.0;
				//cells[i][j].temp = x + y; // for testing the contour plot only
			}
		}
	}

}

void reportResults( std::string const fileName ) {

	std::ofstream file( fileName, std::ofstream::out );

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

void explicitMethod() {

	initDomain();

	simulateExplicit();

	reportResults( "explicit.csv" );

}

void implicitMethod() {

	initDomain();

	simulateImplicit();

	reportResults( "implicit.csv" );

}

void simulateExplicit() {

	bool converged = false;

	int counter = 0;
	while ( !converged )
	{
		++ counter;
		fieldUpdateExplicit();

		converged = isConverged_ts();

		shiftTempsForNextTimeStep();

	}

}

void fieldUpdateExplicit() {

	double dt = std::pow( L / Nx, 2.0 ) * 0.25 * alpha;

	for ( int j = 1; j < Ny; ++j ) { // looping from i/j = 1 to N-1
		for ( int i = 1; i < Nx; ++i ) {
			cells[i][j].temp = alpha * dt * ( ( cells[i+1][j].temp_prev_ts - 2 * cells[i][j].temp_prev_ts + cells[i-1][j].temp_prev_ts ) / std::pow( dx, 2.0 ) \
											+ ( cells[i][j+1].temp_prev_ts - 2 * cells[i][j].temp_prev_ts + cells[i][j-1].temp_prev_ts ) / std::pow( dy, 2.0 ) ) + cells[i][j].temp_prev_ts;
		}
	}

}

bool isConverged_ts() {

	bool converged = true;

	for ( int j = 1; j < Ny; ++j ) {
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

	for ( int j = 1; j < Ny; ++j ) {
		for ( int i = 1; i < Nx; ++i ) {
			auto & thisCell( cells[i][j]);
			if ( abs( thisCell.temp - thisCell.temp_prev_iter ) > errorTolerance ) {
				converged = false;
				return converged;
			}
		}
	}

	return converged;

}

void shiftTempsForNextTimeStep() {

	for ( int j = 1; j < Ny; ++j ) { // looping from i/j = 1 to N-1
		for ( int i = 1; i < Nx; ++i ) {
			cells[i][j].temp_prev_ts = cells[i][j].temp;
		}
	}

}

void shiftTempsForNextIteration() {
	for ( int j = 1; j < Ny; ++j ) { // looping from i/j = 1 to N-1
		for ( int i = 1; i < Nx; ++i ) {
			cells[i][j].temp_prev_iter = cells[i][j].temp;
		}
	}
}

void simulateImplicit() {

	bool finalConverged = false;

	while ( !finalConverged ) {

		bool iterationConverged = false;

		while( !iterationConverged ) {
			fieldUpdateImplicit();

			iterationConverged = isConverged_iter();

			shiftTempsForNextIteration();
		}
	}

}

void fieldUpdateImplicit() {

}
