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

double const errorTolerance = 0.00001;

double const t1 = 12.0;
double const t2 = 24.0;
double const t3 = 36.0;
double const t4 = 48.0;

int main() { // Start of program
    /**
    Docs auto generated with DOxygen
    Program starts here
    */

	// Add runs here
	// ChannelFlowClass( Pe, Re, gamma, numScheme, advDifferencingScheme )
	std::shared_ptr< ChannelFlowClass > run1( new ChannelFlowClass( 2.5, 100.0, 0.2, 0.01, 1.0 ) );
	run1->initAndSim();

}

void ChannelFlowClass::initAndSim()
{
	initDomain();
	simulate();
	reportResults();
}

void ChannelFlowClass::initDomain() // Initialize domain correctly for each run.
{
	// Initial time
	t_start = clock();

	// Set domain dimensions
	L = 2.0; // [cm]
	H = 1.0; // [cm]

	// Set u_max
	u_max = Re * viscosity / H; // [cm/s]

	// Set mesh spacing
	dx = Pe * viscosity / u_max; // [cm]
	dy = dx; // [cm]

	// Set timestep
	dt = 0.001; //gamma * std::pow( dx, 2.0 ) / viscosity; // [s]

	// Calculate number of cells
	Nx = ceil( L / dx );
	Ny = ceil( H / dy );

	// Set p_inlet and p_outlet
	p_outlet = 100; // [g/cm-s-s]
	dp_total = u_max * 8.0 * viscosity * density / std::pow( H, 2.0 ) * L; // [g/cm-s-s]
	p_inlet = p_outlet + dp_total;

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

			// Initialize pressures
			if ( i == 0 ) { // Inlet
				thisCell->p = p_inlet;
			} else if ( i == Nx ) { // Outlet
				thisCell->p = p_outlet;
			} else { // All other cells -- init to get faster convergence
				thisCell->p = p_inlet - ( p_inlet - p_outlet ) * float( i ) / float( Nx );
			}

			// Store cells in vector
			cellVect.push_back( thisCell );

		}
	}

	// Store all cell neighbors on the cell instance for future reference
	for ( int j = 0; j <= Ny; ++j ) {
		for ( int i = 0; i <= Nx; ++i ) {
			auto & thisCell( getCell( i, j ) );

			// Left and right neighbors
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

			// Top and bottom neighbors
			if ( j == 0 ) { // Bottom plate
				thisCell->bottomCell = nullptr;
				thisCell->topCell = getCell( i, j + 1 );
			} else if ( j == Ny ) { // Top plate
				thisCell->bottomCell = getCell( i, j - 1 );
				thisCell->topCell = nullptr;
			} else { // All other cells
				thisCell->bottomCell = getCell( i, j - 1 );
				thisCell->topCell = getCell( i, j + 1 );
			}
		}
	}
}


void ChannelFlowClass::reportResults()
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

	std::ofstream file( Pe_str + "-" + Re_str + "-" + gamma_str + "-" + simTime_str + "-U" + ".csv", std::ofstream::out );

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
				file << std::setprecision(32) << thisCell->u;
			} else {
				file << std::setprecision(32) << thisCell->u << ",";
			}
		}
		file << std::endl;
	}

	file.close();

	reportPressureField();
	reportStatistics();
}

void ChannelFlowClass::reportPressureField()
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

	std::ofstream file( Pe_str + "-" + Re_str + "-" + gamma_str + "-" + simTime_str + "-P" + ".csv", std::ofstream::out );

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
				file << std::setprecision(32) << thisCell->p;
			} else {
				file << std::setprecision(32) << thisCell->p << ",";
			}
		}
		file << std::endl;
	}

	file.close();

}

std::shared_ptr< CellClass >BaseDomainClass::getCell(
	int const i,
	int const j
)
{
	int vectorIndex = i + j * ( Nx + 1 );
	return cellVect[ vectorIndex ];
}

void ChannelFlowClass::reportStatistics()
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

	std::ofstream file( Pe_str + "-" + Re_str + "-" + gamma_str + "-" + simTime_str + ".txt", std::ofstream::out );

	file << "Physical Time (s): " << simTime_str << std::endl;
	file << "Runtime (s): " << float( clock() - t_start ) / CLOCKS_PER_SEC << std::endl;

	file.close();
}

void ChannelFlowClass::simulate()
{

	bool iterateAgain = false;
	bool converged = false;

	reportResults();

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

void ChannelFlowClass::performTimestep( double const timestep ) {

	shiftValsForNewTimestep();

	updatePressureField();

	update_u( timestep );

	update_v( timestep );

	t_curr += timestep;

}

void ChannelFlowClass::update_u( double const timestep ) {

	for ( int j = 1; j < Ny; ++j ) {
		for ( int i = 0; i <= Nx; ++i ) {
			// Get ref to current cell
			auto & thisCell( getCell( i, j ) );

			double advTerm = 0.0;
			// Calculate advective term
			advTerm = thisCell->u_prev_ts * ( thisCell->u_prev_ts - thisCell->leftCell->u_prev_ts ) / dx; \
						+ thisCell->v_prev_ts * ( thisCell->u_prev_ts - thisCell->bottomCell->u_prev_ts ) / dy;

			// Calculate diffusive term
			double diffusiveTerm_partialX = ( thisCell->leftCell->u_prev_ts + thisCell->rightCell->u_prev_ts - 2 * thisCell->u_prev_ts ) / std::pow( dx, 2.0 );
			double diffusiveTerm_partialY = ( thisCell->topCell->u_prev_ts + thisCell->bottomCell->u_prev_ts - 2 * thisCell->u_prev_ts ) / std::pow( dy, 2.0 );

			double pressureTerm;
			// Calculate pressure term
			if ( i == 0 ) {
				double dp_local = thisCell->p - thisCell->rightCell->p;
				double p_leftCell = thisCell->p + dp_local;
				pressureTerm = ( thisCell->rightCell->p - p_leftCell ) / ( 2 * dx );
			} else if ( i == Nx ) {
				double dp_local = thisCell->leftCell->p - thisCell->p;
				double p_rightCell = thisCell->p - dp_local;
				pressureTerm = ( p_rightCell - thisCell->leftCell->p ) / ( 2 * dx );
			} else {
				pressureTerm = ( thisCell->rightCell->p - thisCell->leftCell->p ) / ( 2 * dx );
			}

			// u-momentum equation
			thisCell->u = thisCell->u_prev_ts + timestep * ( viscosity * ( diffusiveTerm_partialX + diffusiveTerm_partialY ) - advTerm - pressureTerm );

		}
	}

}

void ChannelFlowClass::update_v( double const timestep ) {

	for ( int j = 1; j < Ny; ++j ) {
		for ( int i = 0; i <= Nx; ++i ) {
			// Get ref to current cell
			auto & thisCell( getCell( i, j ) );

			double advTerm = 0.0;
			// Calculate advective term
			advTerm = thisCell->u_prev_ts * ( thisCell->v_prev_ts - thisCell->leftCell->v_prev_ts ) / dx \
						+ thisCell->v_prev_ts * ( thisCell->v_prev_ts - thisCell->bottomCell->v_prev_ts ) / dy;

			// Calculate diffusive term
			double diffusiveTerm_partialX = ( thisCell->leftCell->v_prev_ts + thisCell->rightCell->v_prev_ts - 2 * thisCell->v_prev_ts ) / std::pow( dx, 2.0 );
			double diffusiveTerm_partialY = ( thisCell->topCell->v_prev_ts + thisCell->bottomCell->v_prev_ts - 2 * thisCell->v_prev_ts ) / std::pow( dy, 2.0 );

			// Calculate pressure term
			double pressureTerm = ( thisCell->topCell->p - thisCell->bottomCell->p ) / ( 2 * dy );

			// v-momentum equation
			thisCell->v = thisCell->v_prev_ts + timestep * ( viscosity * ( diffusiveTerm_partialX + diffusiveTerm_partialY ) - advTerm - pressureTerm );

		}
	}

}

void ChannelFlowClass::updatePressureField() {

	// Update 'H's
	for ( int j = 1; j < Ny; ++j ) { // Upper and lower plate are not needed
		for ( int i = 1; i < Nx; ++i ) {
			auto & thisCell( getCell( i, j ) );

			double advectionTerm_Hx_E = -thisCell->u * ( thisCell->rightCell->u - thisCell->u ) / dx - thisCell->v * ( thisCell->topCell->u - thisCell->u ) / dy;
			double diffusiveTerm_Hx_E = ( thisCell->rightCell->u + thisCell->leftCell->u - 2 * thisCell->u ) / std::pow( dx, 2.0 )  \
									+ ( thisCell->topCell->u + thisCell->bottomCell->u - 2 * thisCell->u ) / std::pow( dy, 2.0 );

			double advectionTerm_Hx_W = -thisCell->u * ( thisCell->u - thisCell->leftCell->u ) / dx - thisCell->v * ( thisCell->topCell->u - thisCell->u ) / dy;
			double diffusiveTerm_Hx_W = ( thisCell->rightCell->u + thisCell->leftCell->u - 2 * thisCell->u ) / std::pow( dx, 2.0 )  \
									+ ( thisCell->topCell->u + thisCell->bottomCell->u - 2 * thisCell->u ) / std::pow( dy, 2.0 );

			double advectionTerm_Hy_N = -thisCell->u * ( thisCell->rightCell->v - thisCell->v ) / dx - thisCell->v * ( thisCell->topCell->v - thisCell->v ) / dy;
			double diffusiveTerm_Hy_N = ( thisCell->rightCell->v + thisCell->leftCell->v - 2 * thisCell->v ) / std::pow( dx, 2.0 )  \
									+ ( thisCell->topCell->v + thisCell->bottomCell->v - 2 * thisCell->v ) / std::pow( dy, 2.0 );

			double advectionTerm_Hy_S = -thisCell->u * ( thisCell->v - thisCell->leftCell->v ) / dx - thisCell->v * ( thisCell->v - thisCell->bottomCell->v ) / dy;
			double diffusiveTerm_Hy_S = ( thisCell->rightCell->v + thisCell->leftCell->v - 2 * thisCell->v ) / std::pow( dx, 2.0 )  \
									+ ( thisCell->topCell->v + thisCell->bottomCell->v - 2 * thisCell->v ) / std::pow( dy, 2.0 );

			thisCell->Hx_E = advectionTerm_Hx_E + diffusiveTerm_Hx_E / Re;
			thisCell->Hx_W = advectionTerm_Hx_W + diffusiveTerm_Hx_W / Re;
			thisCell->Hy_N = advectionTerm_Hy_N + diffusiveTerm_Hy_N / Re;
			thisCell->Hy_S = advectionTerm_Hy_S + diffusiveTerm_Hy_S / Re;

		}
	}

	int counter = 0;

	while ( !isConvergedPressure() ) {

		shiftPressForNewIter();

		// Update pressure field
		for ( int j = 0; j <= Ny; ++j ) { // Need upper and lower plates
			for ( int i = 1; i < Nx; ++i ) {
				auto & thisCell( getCell( i, j ) );

				if ( j == 0 ) { // Bottom plate
					thisCell->p = thisCell->topCell->p; // dp/dy = 0
				} else if ( j == Ny ) { // Top plate
					thisCell->p = thisCell->bottomCell->p; // dp/dy = 0
				} else { // All other cells
					thisCell->p = ( 1 / 4.0 ) * ( thisCell->rightCell->p + thisCell->leftCell->p + thisCell->topCell->p + thisCell->bottomCell->p ) \
									+ ( dx / 8.0 ) * ( -thisCell->Hx_E - thisCell->Hy_N + thisCell->Hy_S + thisCell->Hx_W );
				}
			}
		}

		++counter;

	}
}

bool  ChannelFlowClass::isConvergedPressure()
{
	double power = 1.0;
	double sumError = 0.0;

	for ( int j = 1; j < Ny; ++j ) {
		for ( int i = 0; i <= Nx; ++i ) {
			auto & thisCell( getCell( i, j ) );
			double thisCellError = std::abs( thisCell->p - thisCell->p_prev_iter );
			sumError += std::pow( thisCellError, power );
		}
	}

	if ( std::pow( sumError, 1 / power ) < errorTolerance ) {
		return true;
	} else {
		return false;
	}

}

void ChannelFlowClass::shiftPressForNewIter() {

	for ( int j = 0; j <= Ny; ++j ) {
		for ( int i = 0; i <= Nx; ++i ) {
			auto & thisCell( getCell( i, j ) );

			thisCell->p_prev_iter = thisCell->p;

		}
	}

}

void ChannelFlowClass::shiftValsForNewTimestep() {

	for ( int j = 1; j < Ny; ++j ) {
		for ( int i = 0; i <= Nx; ++i ) {
			auto & thisCell( getCell( i, j ) );

			thisCell->u_prev_ts = thisCell->u;
			thisCell->v_prev_ts = thisCell->v;

		}
	}

}
