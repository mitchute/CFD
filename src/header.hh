#ifndef header_hh_inluded
#define header_hh_inluded

// Classes

class cell // cell class
{
public:
	double x_loc;
	double y_loc;
	int x_index;
	int y_index;

	double temp;
	double temp_iter;
	double temp_prev_ts;
};

// General Functions
void initDomain();
void reportResults( std::string const fileName );
void reportStatistics( std::string const fileName );
void shiftTempsForNewTimeStep();
void shiftTempsForNewIteration();
bool isConverged_ts();
bool isConverged_iter();
double calcLNorm( double normPower );

// Explicit Functions
void explicitMethod();
void simulateExplicit();
void fieldUpdateExplicit( double dt );

// Implicit Functions
void implicitMethod();
void simulateImplicit();
void fieldUpdateImplicit( double dt );

#endif
