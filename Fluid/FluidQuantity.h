#include <nova/Tools/Grids/Grid.h>

using namespace Nova;

enum
{
	d = 2
};
typedef float T;
typedef Vector<T, d> TV;
typedef Vector<int, d> T_INDEX;

/* Auxiliary Function: Get exact offset for Column-based 1D array */
int offset(const T_INDEX &index, const T_INDEX &counts)
{

	int os = index[1] * counts[0] + index[0];
	if (d == 3)
		os += index[2] * counts[0] * counts[1];
	return os;

	// return (z * xSize * ySize) + (y * xSize) + x;
}

int offset_xyz(int x, int y, int z, const T_INDEX &counts)
{

	int os = y * counts[0] + x;
	if (d == 3)
		os += z * counts[0] * counts[1];
	return os;

	// return (z * xSize * ySize) + (y * xSize) + x;
}

class FluidQuantity
{
  public:
	float *Phi;		// Array of num_cell
	float *Phi_new; // Array of num_cell
	int faceIndicator;	// -1 for Scalar
	Grid<T, d> *grid;

	FluidQuantity();
	FluidQuantity(Grid<T, d> grid, int faceIndicator);
	~FluidQuantity();

	/* Linear Interpolate */
	float linter(float a, float b, float x)
	{
		return (1.0 - x) * a + x * b;
	}

	double at(T_INDEX &index)
	{
		return Phi[offset(index, (*grid).counts)];
	}

	double at_(int x, int y, int z)
	{
		return Phi[offset_xyz(x, y, z, (*grid).counts)];
	}

	/* Compute Velocity */
	TV computeVelocity(T_INDEX &index, FluidQuantity *FluidVelocity[d]);
	void applyAdvection()
	{
		std::swap(Phi, Phi_new);
	}

	/*  */
	void advection();
	void projection();
	void update();
};