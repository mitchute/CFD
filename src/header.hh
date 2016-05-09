#ifndef header_hh_inluded
#define header_hh_inluded

// Classes

class CellClass // cell class
{
public:
	// Location and index
	double x_loc;
	double y_loc;
	int x_index;
	int y_index;

	// u-velocity component
	double u;
	double u_prev_iter;
	double u_prev_ts;

	// v-velocity component
	double v;
	double v_prev_iter;
	double v_prev_ts;

	// Pressure
	double p;
	double p_prev_iter;
	double p_prev_ts;

	// H terms
	double Hx_E;
	double Hx_W;
	double Hy_N;
	double Hy_S;

	// Cell neighbors
	std::shared_ptr< CellClass > leftCell;
	std::shared_ptr< CellClass > rightCell;
	std::shared_ptr< CellClass > topCell;
	std::shared_ptr< CellClass > bottomCell;

	// Constructor
	CellClass() :
		x_loc( 0.0 ),
		y_loc( 0.0 ),
		x_index( 0 ),
		y_index( 0 ),

		u( 0.0 ),
		u_prev_iter( 0.0 ),
		u_prev_ts( 0.0 ),

		v( 0.0 ),
		v_prev_iter( 0.0 ),
		v_prev_ts( 0.0 ),

		p( 0.0 ),
		p_prev_iter( 0.0 ),
		p_prev_ts( 0.0 ),

		leftCell( nullptr ),
		rightCell( nullptr ),
		bottomCell( nullptr ),
		topCell( nullptr )

	{};

	// Destructor
	~CellClass(){};

};

class BaseDomainClass
{
public:

	// BaseDomain Members
	double u_max;
	double p_inlet;
	double dp_total;
	double p_outlet;
	double L;
	double H;
	double dx;
	double dy;
	double dt;
	int Nx;
	int Ny;
	clock_t t_start;
	double t_curr;
	double t_end;
	std::vector< std::shared_ptr< CellClass > > cellVect;

	// Default constructor
	BaseDomainClass() :
		t_curr( 0.0 ),
		t_end( 60.0 )
	{};

	// Destructor
	~BaseDomainClass(){};

	// Pure virtual functions
	virtual void initAndSim()=0;
	virtual void initDomain()=0;
	virtual void simulate()=0;
	virtual void reportResults()=0;
	virtual void reportStatistics()=0;
	virtual bool isConverged_ts()=0;
	virtual bool isConverged_iter()=0;
	virtual double calcLNorm_ts( double const power )=0;
	virtual double calcLNorm_iter( double const power )=0;
	virtual void shiftValsForNewTimestep()=0;
	virtual void shiftValsForNewIteration()=0;

	// BaseDomain member functions
	std::shared_ptr< CellClass > getCell( int const i, int const j );

};

class ChannelFlowClass : public BaseDomainClass
{
public:
	double Pe;
	double Re;
	double gamma;
	double viscosity;
	double density;

	// Default constructor
	ChannelFlowClass() :
		Pe( 0.0 ),
		Re( 0.0 ),
		gamma( 0.0 ),
		viscosity( 0.0 ),
		density( 0.0 )
	{}

	// Member constructor
	ChannelFlowClass(
		double const Pe_,
		double const Re_,
		double const gamma_,
		double const viscosity_,
		double const density_
	) :
		Pe( Pe_ ),
		Re( Re_ ),
		gamma( gamma_ ),
		viscosity( viscosity_ ),
		density( density_ )
	{}

	// Destructor
	~ChannelFlowClass(){};

	// Virtual functions
	void initAndSim();
	void initDomain();
	void simulate();
	void reportResults();
	void reportStatistics();
	bool isConverged_ts();
	bool isConverged_iter();
	double calcLNorm_ts( double const power );
	double calcLNorm_iter( double const power );
	void shiftValsForNewTimestep();
	void shiftValsForNewIteration();

	// Member Functions
	void updatePressureField();
	void update_u( double dt );
	void update_v( double dt );
	void performTimestep( double dt );
	void reportPressureField();
	bool isConvergedPressure();
	void shiftPressForNewIter();
};

#endif
