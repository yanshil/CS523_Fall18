Class FluidQuantity
{
	bool isScalar;
	double * Phi;	// Matrix
	double * Phi_new;	// Matrix
	
	FluidQuantity();
	~FluidQuantity();
	
	/* Linear Interpolate */
	double linter();
	/* Compute Velocity */
	double computeVelocity();
	void applyAdvection();
	
	/*  */
	void advection();
	void projection();
	void update();
}