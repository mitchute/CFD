// C++ headers
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <sstream>
#include <string>
#include <time.h>

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
double const gamma = 1.5;

double const errorTolerance = 0.0001;

double LNorm = 0.0;

double t1= 0.05;
double t2 = 0.1;
double t3 = 0.15;
double t4 = 0.2;

clock_t t_init;
double t_final = 2.0;

// Create cell array
cell cells[Nx+1][Ny+1];

int main() { // Start of program
    /**
    Test a comment here.
    This is where the program starts
    */
	//explicitMethod();
	implicitMethod();

}

void initDomain() { // Initialize domain correctly for each run.

	t_init = clock();

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

void reportResults(
	std::string const solutionScheme,
	double simTime
)
{ // Report results to file

	// Convert simTime to String
	std::stringstream simTime_stream;
	simTime_stream << std::fixed << std::setprecision(2) << simTime;
	std::string simTime_str = simTime_stream.str();

	// Convert gamma to String
	std::stringstream gamma_stream;
	gamma_stream << std::fixed << std::setprecision(2) << gamma;
	std::string gamma_str = gamma_stream.str();

	std::ofstream file( solutionScheme + "-" + gamma_str + "-" + simTime_str + "-" + std::to_string( Nx ) + ".csv", std::ofstream::out );

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

void reportStatistics(
	std::string const solutionScheme,
	double simTime
)
{	// Report statistics from run

	double l1 = calcLNorm_exact( 1.0 );
	double l2 = calcLNorm_exact( 2.0 );

	// Convert simTime to String
	std::stringstream simTime_stream;
	simTime_stream << std::fixed << std::setprecision(2) << simTime;
	std::string simTime_str = simTime_stream.str();

	// Convert gamma to String
	std::stringstream gamma_stream;
	gamma_stream << std::fixed << std::setprecision(2) << gamma;
	std::string gamma_str = gamma_stream.str();

	std::ofstream file( solutionScheme + "-" + gamma_str + "-" + simTime_str + "-" + std::to_string( Nx ) + ".txt", std::ofstream::out );

	file << "L1 Norm: " << l1 << std::endl;
	file << "L2 Norm: " << l2 << std::endl;
	file << "Physical Time (s): " << simTime_str << std::endl;
	file << "Runtime (s): " << float( clock() - t_init ) / CLOCKS_PER_SEC << std::endl;

	file.close();
}

void explicitMethod() { // Simulate explicit method

	double simTime = 0.0;

	initDomain();
	simTime = simulateExplicit();
	reportResults( "explicit", simTime );
	reportStatistics( "explicit", simTime );
}

void implicitMethod() { // Simulate implicit method

	double simTime = 0.0;

	initDomain();
	simTime = simulateImplicit();
	reportResults( "implicit", simTime );
	reportStatistics( "implicit", simTime );
}

double simulateExplicit() {

	double dt = gamma * std::pow( dx, 2.0 ) / alpha;
	double t = 0.0;

	for ( t = 0; t < t_final; t = t + dt )
	{
		fieldUpdateExplicit( dt );
		shiftTempsForNewTimeStep();

		if ( t < t1 && t + dt > t1 ) {
			dt = t1 - t;
			t += dt;
			fieldUpdateExplicit( dt );
			shiftTempsForNewTimeStep();
			reportResults( "explicit", t );
			reportStatistics( "explicit", t );
		} else if ( t < t2 && t + dt > t2 ) {
			dt = t2 - t;
			t += dt;
			fieldUpdateExplicit( dt );
			shiftTempsForNewTimeStep();
			reportResults( "explicit", t );
			reportStatistics( "explicit", t );
		} else if ( t < t3 && t + dt > t3 ) {
			dt = t3 - t;
			t += dt;
			fieldUpdateExplicit( dt );
			shiftTempsForNewTimeStep();
			reportResults( "explicit", t );
			reportStatistics( "explicit", t );
		} else if ( t < t4 && t + dt > t4 ) {
			dt = t4 - t;
			t += dt;
			fieldUpdateExplicit( dt );
			shiftTempsForNewTimeStep();
			reportResults( "explicit", t );
			reportStatistics( "explicit", t );
		}

		dt = gamma * std::pow( dx, 2.0 ) / alpha;

	}

	return t;
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

double simulateImplicit() {

	double dt = ( gamma * std::pow( dx, 2.0 ) ) / alpha;
	double t = 0.0;

	for ( t = 0; t < t_final; t += dt ) {

		performTimeStepImplicit( dt );

		if ( t < t1 && t + dt > t1 ) {
			dt = t1 - t;
			t += dt;
			performTimeStepImplicit( dt );
			reportResults( "implicit", t );
			reportStatistics( "implicit", t );
		} else if ( t < t2 && t + dt > t2 ) {
			dt = t2 - t;
			t += dt;
			performTimeStepImplicit( dt );
			reportResults( "implicit", t );
			reportStatistics( "implicit", t );
		} else if ( t < t3 && t + dt > t3 ) {
			dt = t3 - t;
			t += dt;
			performTimeStepImplicit( dt );
			reportResults( "implicit", t );
			reportStatistics( "implicit", t );
		} else if ( t < t4 && t + dt > t4 ) {
			dt = t4 - t;
			t += dt;
			performTimeStepImplicit( dt );
			reportResults( "implicit", t );
			reportStatistics( "implicit", t );
		}

		dt = gamma * std::pow( dx, 2.0 ) / alpha;

	}

	return t;
}

void performTimeStepImplicit(
	double dt
)
{
	bool converged = false;

	while( !converged ) { // Iteration loop
		shiftTempsForNewIteration();
		fieldUpdateImplicit( dt );

		if ( calcLNorm_iter( 1.0 ) < errorTolerance) {
			converged = true;
		}
	}

	shiftTempsForNewTimeStep();

}

void fieldUpdateImplicit(
	double dt
)
{
	double delta = alpha * dt / dx;
	double eta = alpha * dt / dy;

	for ( int j = 1; j < Ny; ++j ) { // looping from i/j = 1 to N-1
		for ( int i = 1; i < Nx; ++i ) {
			cells[i][j].temp =  ( cells[i][j].temp_prev_ts + delta * ( cells[i+1][j].temp + cells[i-1][j].temp ) \
								+ eta * ( cells[i][j+1].temp + cells[i][j-1].temp ) ) / ( 1 + ( 2 * delta ) + ( 2 * eta ) );
		}
	}
}

double calcLNorm_exact(
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

double calcLNorm_iter(
	double normPower
)
{
	double sumError = 0.0;

	for ( int j = 1; j < Ny; ++j ) {
		for ( int i = 1; i < Nx; ++i ) {
			auto & c( cells[i][j] );
			double err = (c.temp - c.temp_iter );
			sumError += std::abs( std::pow( err, normPower ) );
		}
	}

	return std::pow( sumError, 1 / normPower );
}
