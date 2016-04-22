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
	double u_iter;
	double u_prev_iter;
	double u_prev_ts;

	// v-velocity component
	double v;
	double v_iter;
	double v_prev_iter;
	double v_prev_ts;

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
		u_iter( 0.0 ),
		u_prev_iter( 0.0 ),
		u_prev_ts( 10.0 ),

		v( 0.0 ),
		v_iter( 0.0 ),
		v_prev_iter( 0.0 ),
		v_prev_ts( 10.0 ),

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
	double L;
	double H;
	double dx;
	double dy;
	double dt;
	int Nx;
	int Ny;
	clock_t t_start;
	double t_curr;
	std::vector< std::shared_ptr< CellClass > > cellVect;

	// Default constructor
	BaseDomainClass() :
		t_curr( 0.0 )
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

class BurgersClass : public BaseDomainClass
{
public:
	double Pe;
	double Re;
	double gamma;
	double viscosity;
	std::string numScheme; // Explicit or implicit
	std::string advDifferencingScheme; // FOU or CD

	// Default constructor
	BurgersClass() :
		viscosity( 0.0 ),
		Pe( 0.0 ),
		Re( 0.0 ),
		gamma( 0.0 ),
		numScheme( "default" ),
		advDifferencingScheme( "NYBD" )
	{}

	// Member constructor
	BurgersClass(
		double const Pe_,
		double const Re_,
		double const gamma_,
		std::string const numScheme_,
		std::string const advDifferencingScheme_
	) :
		Pe( Pe_ ),
		Re( Re_ ),
		gamma( gamma_ ),
		numScheme( numScheme_ ),
		advDifferencingScheme( advDifferencingScheme_ )
	{}

	// Destructor
	~BurgersClass(){};

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
	void update_u( double dt );
	void update_v( double dt );

};

#endif
