#ifndef header_hh_inluded
#define header_hh_inluded

// Classes

class cell
{
public:
	double x_loc;
	double y_loc;
	int x_index;
	int y_index;

	double temp = 0.0;
	double temp_prev_iter = 0.0;
	double temp_prev_ts = 0.0;
};

// General Functions
void initDomain();
void reportResults( std::string const fileName );
void shiftTempsForNextTimeStep();
void shiftTempsForNextIteration();
bool isConverged_ts();
bool isConverged_iter();

// Explicit Functions
void explicitMethod();
void simulateExplicit();
void fieldUpdateExplicit();

// Implicit Functions
void implicitMethod();
void simulateImplicit();
void fieldUpdateImplicit();

#endif
