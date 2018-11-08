#include <nova/Tools/Grids/Grid.h>

using namespace Nova;

enum
{
    d = 2
};

typedef float T;
typedef Vector<T, d> TV;
typedef Vector<int, d> T_INDEX;

// /*!
//  * Auxiliary Function: Get exact offset for Column-based 1D array
//  * return (z * xSize * ySize) + (y * xSize) + x;
//  */
// int index2offset(const T_INDEX &index, const T_INDEX &counts);

class FluidQuantity
{
  public:
    float *Phi;     // Array of num_cell
    float *Phi_new; // Array of num_cell
    int axis;       // -1 for Scalar. 0, 1, 2 to indicate faces axis
    int number_of_ghost_cells;
    Grid<T, d> *grid;

    FluidQuantity();
    FluidQuantity(Grid<T, d> &grid, int axis);
    ~FluidQuantity();

    // ********************************************

    /*!
     * 1D Linear Interpolate etween a and b for x in (0, 1)
     */
    float linter(float a, float b, float x)
    {
        return (1.0 - x) * a + x * b;
    }

    /*!
 * Auxiliary Function: Get exact offset for Column-based 1D array
 * return (z * xSize * ySize) + (y * xSize) + x;
 */
    int index2offset(const T_INDEX &index)
    {

        int os = index[1] * (*grid).counts[0] + index[0];
        if (d == 3)
            os += index[2] * (*grid).counts[0] * (*grid).counts[1];
        return os;
    }

    /*!
     * Linear Interpolator for TV{i, j, k} on grid
     * Coordinates will be clamped to lie in simulation domain
     */
    float linter(TV &location);

    /*!
     * Read access in the quantity field
     */
    float at(const T_INDEX &index)
    {
        return Phi[index2offset(index)];
    }

    /*!
     * Read & write access in the quantity field
     */
    float &modify_at(const T_INDEX &index)
    {
        return Phi[index2offset(index)];
    }

    // a.modify_at(1) = ...;

    /* Compute Velocity */
    TV computeVelocity(T_INDEX &index, FluidQuantity *FluidVelocity[d]);
    TV traceBack(TV &location, float timestep, TV &velocity);

    void flip()
    {
        std::swap(Phi, Phi_new);
    }

    /*  */
    void advection(float timestep, FluidQuantity *FluidVelocity[d]);
    void projection();
    void update();
};
