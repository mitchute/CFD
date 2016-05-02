// C++ headers
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <memory>
#include <stdio.h>
#include <sstream>
#include <string>
#include <time.h>
#include <vector>

// My headers
#include "header.hh"

double const u_lid = 1;
double const errorTolerance = 0.0001;

double const t1 = 0.2;
double const t2 = 0.5;
double const t3 = 1.0;
double const t4 = 3.0;

int main() { // Start of program
    /**
    Docs auto generated with DOxygen
    Program starts here
    */

	// Add runs here
	// BurgersClass( Pe, Re, gamma, numScheme, advDifferencingScheme )
	std::shared_ptr< BurgersClass > run1( new BurgersClass( 1.0, 40.0, 0.1, "Explicit", "FOU" ) );
	run1->initAndSim();

	std::shared_ptr< BurgersClass > run2( new BurgersClass( 1.0, 40.0, 0.1, "Explicit", "CD" ) );
	run2->initAndSim();

	std::shared_ptr< BurgersClass > run3( new BurgersClass( 1.0, 40.0, 0.1, "Implicit", "FOU" ) );
	run3->initAndSim();

	std::shared_ptr< BurgersClass > run4( new BurgersClass( 1.0, 40.0, 0.1, "Implicit", "CD" ) );
	run4->initAndSim();

	std::shared_ptr< BurgersClass > run5( new BurgersClass( 5.0, 40.0, 0.1, "Explicit", "FOU" ) );
	run5->initAndSim();

	std::shared_ptr< BurgersClass > run6( new BurgersClass( 5.0, 40.0, 0.1, "Explicit", "CD" ) );
	run6->initAndSim();

	std::shared_ptr< BurgersClass > run7( new BurgersClass( 5.0, 40.0, 0.1, "Implicit", "FOU" ) );
	run7->initAndSim();

	std::shared_ptr< BurgersClass > run8( new BurgersClass( 5.0, 40.0, 0.1, "Implicit", "CD" ) );
	run8->initAndSim();

	std::shared_ptr< BurgersClass > run9( new BurgersClass( 5.0, 40.0, 1.0, "Implicit", "FOU" ) );
	run9->initAndSim();

	std::shared_ptr< BurgersClass > run10( new BurgersClass( 5.0, 200.0, 1.0, "Implicit", "FOU" ) );
	run10->initAndSim();

}

void BurgersClass::initAndSim()
{
	initDomain();
	simulate();
	reportResults();
}

void BurgersClass::initDomain() // Initialize domain correctly for each run.
{
	// Initial time
	t_start = clock();

	// Set domain dimensions
	L = 2.0;
	H = 1.0;

	// Set mesh spacing
	dx = Pe * H / Re;
	dy = dx;

	// Set timestep
	dt = ( gamma * std::pow( Pe, 2.0 ) * H ) / ( u_lid * Re );

	// Set viscosity
	viscosity = u_lid * H / Re;

	// Calculate number of cells
	Nx = ceil( L / dx );
	Ny = Nx / 2;

	// Initialize all cells
	for ( int j = 0; j <= Ny; ++j ) {

		double y = double(j) * dy;

		for ( int i = 0; i <= Nx; ++i ) {

			std::shared_ptr< CellClass > thisCell( new CellClass );

			double x = double( i ) * dx;

			thisCell->x_index = i;
			thisCell->y_index = j;
			thisCell->x_loc = x;
			thisCell->y_loc = y;

			if ( j == Ny ) { // Top plate
				double initVal = u_lid;
				// u-velocity component
				thisCell->u = initVal;
				thisCell->u_prev_iter = initVal;
				thisCell->u_prev_ts = initVal;
				// v-velocity component
				thisCell->v = 0.0;
				thisCell->v_prev_iter = 0.0;
				thisCell->v_prev_ts = 0.0;
			} else { // Everywhere else initialized to 0.0
				thisCell->u = 0.0;
				thisCell->u_prev_iter = 0.0; // initialized big
				thisCell->u_prev_ts = 0.0;
			}

			// Store cells in vector
			cellVect.push_back( thisCell );

		}
	}

	// Store all cell neighbors on the cell instance for future reference
	for ( int j = 1; j < Ny; ++j ) {
		for ( int i = 0; i <= Nx; ++i ) {
			auto & thisCell( getCell( i, j) );

			thisCell->bottomCell = getCell( i, j - 1 );
			thisCell->topCell = getCell( i, j + 1 );

			if ( i == 0 ) { // Left side periodic boundary condition
				thisCell->leftCell = getCell( Nx, j );
				thisCell->rightCell = getCell( i + 1, j );
			} else if ( i == Nx ) { // Right side periodic boundary condition
				thisCell->leftCell = getCell( i - 1, j );
				thisCell->rightCell = getCell( 0 , j );
			} else { // All other cells
				thisCell->leftCell = getCell( i - 1, j );
				thisCell->rightCell = getCell( i + 1, j );
			}
		}
	}
}


void BurgersClass::reportResults()
{ // Report results to file

	// Convert Pe to string
	std::stringstream Pe_stream;
	Pe_stream << std::fixed << std::setprecision(2) << Pe;
	std::string Pe_str = Pe_stream.str();

	// Convert Re to string
	std::stringstream Re_stream;
	Re_stream << std::fixed << std::setprecision(2) << Re;
	std::string Re_str = Re_stream.str();

	// Convert gamma to string
	std::stringstream gamma_stream;
	gamma_stream << std::fixed << std::setprecision(2) << gamma;
	std::string gamma_str = gamma_stream.str();

	// Convert simTime to string
	std::stringstream simTime_stream;
	simTime_stream << std::fixed << std::setprecision(2) << t_curr;
	std::string simTime_str = simTime_stream.str();

	std::ofstream file( numScheme + "-" + advDifferencingScheme + "-" + Pe_str + "-" + Re_str + "-" + gamma_str + "-" + simTime_str + ".csv", std::ofstream::out );

	for ( int j = 0; j <= Ny; ++j ) {
		if ( j == 0 ) { // print x locations

			file << ","; // pad space to get correct alignment
			for ( int i = 0; i <= Nx; ++i ) {
				auto & thisCell( getCell( i, j ) );

				if ( i == Nx ) {
					file << thisCell->x_loc << std::endl;
				} else {
					file << thisCell->x_loc << ",";
				}
			}
		}

		for ( int i = 0; i <= Nx; ++i ) {
			auto & thisCell( getCell( i, j ) );
			if ( i == 0 ) { // print y locations
				file << thisCell->y_loc << ",";
			}

			if ( i == Nx ) {
				file << thisCell->u;
			} else {
				file << thisCell->u << ",";
			}
		}
		file << std::endl;
	}

	file.close();

	reportStatistics();
}

std::shared_ptr< CellClass >BaseDomainClass::getCell(
	int const i,
	int const j
)
{
	int vectorIndex = i + j * ( Nx + 1 );
	return cellVect[ vectorIndex ];
}

void BurgersClass::reportStatistics()
{	// Report statistics from run

	// Convert Pe to string
	std::stringstream Pe_stream;
	Pe_stream << std::fixed << std::setprecision(2) << Pe;
	std::string Pe_str = Pe_stream.str();

	// Convert Re to string
	std::stringstream Re_stream;
	Re_stream << std::fixed << std::setprecision(2) << Re;
	std::string Re_str = Re_stream.str();

	// Convert gamma to string
	std::stringstream gamma_stream;
	gamma_stream << std::fixed << std::setprecision(2) << gamma;
	std::string gamma_str = gamma_stream.str();

	// Convert simTime to string
	std::stringstream simTime_stream;
	simTime_stream << std::fixed << std::setprecision(2) << t_curr;
	std::string simTime_str = simTime_stream.str();

	std::ofstream file( numScheme + "-" + advDifferencingScheme + "-" + Pe_str + "-" + Re_str + "-" + gamma_str + "-" + simTime_str + ".txt", std::ofstream::out );

	file << "Physical Time (s): " << simTime_str << std::endl;
	file << "Runtime (s): " << float( clock() - t_start ) / CLOCKS_PER_SEC << std::endl;

	file.close();
}

void BurgersClass::simulate()
{

	bool iterateAgain = false;
	bool converged = false;

	while ( !converged || iterateAgain ) {

		if ( t_curr < t1 && t_curr + dt > t1 ) {
			performTimestep( t1 - t_curr );
			reportResults();
			iterateAgain = true;
		} else if ( t_curr < t2 && t_curr + dt > t2 ) {
			performTimestep( t2 - t_curr );
			reportResults();
			iterateAgain = true;
		} else if ( t_curr < t3 && t_curr + dt > t3 ) {
			performTimestep( t3 - t_curr );
			reportResults();
			iterateAgain = true;
		} else if ( t_curr < t4 && t_curr + dt > t4 ) {
			performTimestep( t4 - t_curr );
			reportResults();
			iterateAgain = true;
		} else if ( t_curr == t1  || t_curr == t2 || t_curr == t3 || t_curr == t4 ) {
			performTimestep( dt );
			iterateAgain = false;
			reportResults();
		} else {
			performTimestep( dt );
			iterateAgain = false;
		}

		// Check convergence
		converged = ( t_end < t_curr );

	}

}

void BurgersClass::performTimestep( double const timestep ) {

	shiftValsForNewTimestep();

	if ( numScheme == "Explicit" ) {

		update_u( timestep );
		update_v( timestep );

	} else if ( numScheme == "Implicit" ) {

		bool converged = false;

		while( !converged ) {

			shiftValsForNewIteration();

			update_u( timestep );
			update_v( timestep );

			converged = isConverged_iter();

		}

	}

	t_curr += timestep;

}

void BurgersClass::update_u( double const timestep ) {

	if ( numScheme == "Explicit" ) {

		for ( int j = 1; j < Ny; ++j ) {
			for ( int i = 0; i <= Nx; ++i ) {
				// Get ref to current cell
				auto & thisCell( getCell( i, j ) );

				double advTerm = 0.0;
				// Calculate advective term
				if ( advDifferencingScheme == "FOU" ) {
					advTerm = thisCell->u_prev_ts * ( thisCell->u_prev_ts - thisCell->leftCell->u_prev_ts ) / dx; \
								+ thisCell->v_prev_ts * ( thisCell->u_prev_ts - thisCell->bottomCell->u_prev_ts ) / dy;
				} else if ( advDifferencingScheme == "CD" ) {
					advTerm = thisCell->u_prev_ts * ( thisCell->rightCell->u_prev_ts - thisCell->leftCell->u_prev_ts ) / ( 2.0 * dx ) \
								+ thisCell->v_prev_ts * ( thisCell->topCell->u_prev_ts - thisCell->bottomCell->u_prev_ts ) / ( 2.0 * dy );
				}

				// Calculate diffusive term
				double diffusiveTerm_partialX = ( thisCell->leftCell->u_prev_ts + thisCell->rightCell->u_prev_ts - 2 * thisCell->u_prev_ts ) / std::pow( dx, 2.0 );
				double diffusiveTerm_partialY = ( thisCell->topCell->u_prev_ts + thisCell->bottomCell->u_prev_ts - 2 * thisCell->u_prev_ts ) / std::pow( dy, 2.0 );

				thisCell->u = thisCell->u_prev_ts + timestep * ( viscosity * ( diffusiveTerm_partialX + diffusiveTerm_partialY ) - advTerm );

			}
		}

	} else if ( numScheme == "Implicit" ) {

		bool converged = false;

		for ( int j = 1; j < Ny; ++j ) {
			for ( int i = 0; i <= Nx; ++i ) {
				// Get ref to current cell
				auto & thisCell( getCell( i, j ) );

				double U = timestep * thisCell->u / dx;
				double V = timestep * thisCell->v / dy;
				double gamma_x = timestep * viscosity / std::pow( dx, 2.0 );
				double gamma_y = timestep * viscosity / std::pow( dx, 2.0 );

				double numerator = 0.0;
				double denominator = 0.0;

				if ( advDifferencingScheme == "FOU" ) {
					numerator = thisCell->u_prev_ts + U * thisCell->leftCell->u + V * thisCell->bottomCell->u \
								+ gamma_x * ( thisCell->rightCell->u + thisCell->leftCell->u ) + gamma_y * ( thisCell->topCell->u + thisCell->bottomCell->u );
					denominator = 1 + U + V + 2 * gamma_x + 2 * gamma_y;
				} else if ( advDifferencingScheme == "CD" ) {
					numerator = thisCell->u_prev_ts - U * ( thisCell->rightCell->u - thisCell->leftCell->u ) - V * ( thisCell->topCell->u - thisCell->bottomCell->u ) \
								+ gamma_x * ( thisCell->rightCell->u + thisCell->leftCell->u ) + gamma_y * ( thisCell->topCell->u + thisCell->bottomCell->u );
					denominator = 1 + 2 * gamma_x + 2 * gamma_y;
				}

				thisCell->u = numerator / denominator;

			}
		}
	}

}

void BurgersClass::update_v( double const timestep ) {

	if ( numScheme == "Explicit" ) {

		for ( int j = 1; j < Ny; ++j ) {
			for ( int i = 0; i <= Nx; ++i ) {
				// Get ref to current cell
				auto & thisCell( getCell( i, j ) );

				double advTerm = 0.0;
				// Calculate advective term
				if ( advDifferencingScheme == "FOU" ) {
					advTerm = thisCell->u_prev_ts * ( thisCell->v_prev_ts - thisCell->leftCell->v_prev_ts ) / dx \
								+ thisCell->v_prev_ts * ( thisCell->v_prev_ts - thisCell->bottomCell->v_prev_ts ) / dy;
				} else if ( advDifferencingScheme == "CD" ) {
					advTerm = thisCell->u_prev_ts * ( thisCell->rightCell->v_prev_ts - thisCell->leftCell->v_prev_ts ) / ( 2.0 * dx ); \
								+ thisCell->v_prev_ts * ( thisCell->topCell->v_prev_ts - thisCell->bottomCell->v_prev_ts ) / ( 2.0 * dy );
				}

				// Calculate diffusive term
				double diffusiveTerm_partialX = ( thisCell->leftCell->v_prev_ts + thisCell->rightCell->v_prev_ts - 2 * thisCell->v_prev_ts ) / std::pow( dx, 2.0 );
				double diffusiveTerm_partialY = ( thisCell->topCell->v_prev_ts + thisCell->bottomCell->v_prev_ts - 2 * thisCell->v_prev_ts ) / std::pow( dy, 2.0 );

				thisCell->v = thisCell->v_prev_ts + timestep * ( viscosity * ( diffusiveTerm_partialX + diffusiveTerm_partialY ) - advTerm );

			}
		}

	} else if ( numScheme == "Implicit" ) {

		bool converged = false;

		for ( int j = 1; j < Ny; ++j ) {
			for ( int i = 0; i <= Nx; ++i ) {
				// Get ref to current cell
				auto & thisCell( getCell( i, j ) );

				double U = timestep * thisCell->u / dx;
				double V = timestep * thisCell->v / dy;
				double gamma_x = timestep * viscosity / std::pow( dx, 2.0 );
				double gamma_y = timestep * viscosity / std::pow( dx, 2.0 );

				double numerator = 0.0;
				double denominator = 0.0;

				if ( advDifferencingScheme == "FOU" ) {
					numerator = thisCell->v_prev_ts + U * thisCell->leftCell->v + V * thisCell->bottomCell->v \
								+ gamma_x * ( thisCell->rightCell->v + thisCell->leftCell->v ) + gamma_y * ( thisCell->topCell->v + thisCell->bottomCell->v );
					denominator = 1 + U + V + 2 * gamma_x + 2 * gamma_y;
				} else if ( advDifferencingScheme == "CD" ) {
					numerator = thisCell->v_prev_ts - U * ( thisCell->rightCell->v - thisCell->leftCell->v ) - V * ( thisCell->topCell->v - thisCell->bottomCell->v ) \
								+ gamma_x * ( thisCell->rightCell->v + thisCell->leftCell->v ) + gamma_y * ( thisCell->topCell->v + thisCell->bottomCell->v );
					denominator = 1 + 2 * gamma_x + 2 * gamma_y;
				}

				thisCell->v = numerator / denominator;

			}
		}
	}

}

bool BurgersClass::isConverged_ts() {

	if ( calcLNorm_ts( 1.0 ) < errorTolerance ) {
		return true;
	} else {
		return false;
	}

}

bool BurgersClass::isConverged_iter() {

	if ( calcLNorm_iter( 1.0 ) < errorTolerance ) {
		return true;
	} else {
		return false;
	}

}

double BurgersClass::calcLNorm_ts(
	double const power
)
{
	double sumError = 0.0;

	for ( int j = 1; j < Ny; ++j ) {
		for ( int i = 0; i <= Nx; ++i ) {
			auto & thisCell( getCell( i, j ) );
			double thisCellError = std::abs( thisCell->u - thisCell->u_prev_ts );
			thisCellError += std::abs( thisCell->v - thisCell->v_prev_ts );
			sumError += std::pow( thisCellError, power );
		}
	}

	return std::pow( sumError, 1 / power );

}

double BurgersClass::calcLNorm_iter(
	double power
)
{
	double sumError = 0.0;

	for ( int j = 1; j < Ny; ++j ) {
		for ( int i = 0; i <= Nx; ++i ) {
			auto & thisCell( getCell( i, j ) );
			double thisCellError = std::abs( thisCell->u - thisCell->u_prev_iter );
			thisCellError += std::abs( thisCell->v - thisCell->v_prev_iter );
			sumError += std::pow( thisCellError, power );
		}
	}

	return std::pow( sumError, 1 / power );

}

void BurgersClass::shiftValsForNewTimestep() {

	for ( int j = 1; j < Ny; ++j ) {
		for ( int i = 0; i <= Nx; ++i ) {
			auto & thisCell( getCell( i, j ) );

			thisCell->u_prev_ts = thisCell->u;
			thisCell->v_prev_ts = thisCell->v;
		}
	}

}

void BurgersClass::shiftValsForNewIteration() {

	for ( int j = 1; j < Ny; ++j ) {
		for ( int i = 0; i <= Nx; ++i ) {
			auto & thisCell( getCell( i, j ) );

			thisCell->u_prev_iter = thisCell->u;
			thisCell->v_prev_iter = thisCell->v;
		}
	}

}
