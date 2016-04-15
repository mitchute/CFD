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
void reportResults( std::string const solutionScheme, double simTime );
void reportStatistics( std::string const solutionScheme, double simTime );
void shiftTempsForNewTimeStep();
void shiftTempsForNewIteration();
bool isConverged_ts();
bool isConverged_iter();
double calcLNorm_exact( double normPower );
double calcLNorm_iter( double normPower );

// Explicit Functions
void explicitMethod();
double simulateExplicit();
void fieldUpdateExplicit( double dt );

// Implicit Functions
void implicitMethod();
double simulateImplicit();
void performTimeStepImplicit( double dt );
void fieldUpdateImplicit( double dt );

#endif
