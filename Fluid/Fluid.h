#include <nova/Tools/Grids/Grid.h>
#include <nova/Tools/Utilities/Range_Iterator.h>

using namespace Nova;

enum
{
    d = 2
};

typedef double T;
typedef Vector<T, d> TV;
typedef Vector<int, d> T_INDEX;

//----------Auxiliary Function--------------------
T_INDEX Next_Cell(const int axis, const T_INDEX &index);
T_INDEX Previous_Cell(const int axis, const T_INDEX &index);

//------------------------------------------------

class FluidQuantity
{
  public:
    double *Phi;     // Array of num_cell
    double *Phi_new; // Array of num_cell
    int axis;        // -1 for Scalar. 0, 1, 2 to indicate faces axis
    int number_of_ghost_cells;
    Grid<T, d> *grid;

    FluidQuantity()
    {
        grid = nullptr;
        std::cout << "Wrong Constructor!" << std::endl;
    }
    FluidQuantity(Grid<T, d> &grid, int axis);
    ~FluidQuantity();

    // ********************************************

    /*!
     * 1D Linear Interpolate etween a and b for x in (0, 1)
     */
    double linter(double a, double b, double x)
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
    double linter(TV &location);

    /*!
     * Read and Write access in the quantity field for Phi
     */
    double &at(const T_INDEX &index)
    {
        return Phi[index2offset(index)];
    }

    /*!
     * Read & Write access in the quantity field for Phi_new
     */
    double &new_at(const T_INDEX &index)
    {
        return Phi_new[index2offset(index)];
    }

    /* Compute Velocity */
    TV computeVelocity(const T_INDEX &index, FluidQuantity *velocityField[d]);
    TV Clamp_To_Domain(const TV &location);
    TV traceBack(const TV &location, double timestep, TV &velocity);

    void flip()
    {
        std::swap(Phi, Phi_new);
    }

    /* Advection */
    void advect(const T_INDEX &index, double timestep, FluidQuantity *velocityField[d]);
};

class FluidSolver
{
    FluidQuantity *velocityField[d];
    FluidQuantity *density_field;

    // FluidQuantity pressure;
    Grid<T, d> *grid;

    double density;

    double *rhs;
    double *pressure_solution;

  public:
    FluidSolver();
    FluidSolver(Grid<T, d> &grid, double density)
    {
        std::cout << "Constructor of FluidSolver" << std::endl;

        for (int axis = 0; axis < d; axis++)
            velocityField[axis] = new FluidQuantity(grid, axis);

        // TODO
        this->density_field = new FluidQuantity(grid, -1);
        this->density = density;
        this->grid = &grid;

        // density = FluidQuantity(grid, -1);
        // pressure = FluidQuantity(grid, -1);

        rhs = new double[grid.counts.Product()];
        pressure_solution = new double[grid.counts.Product()];
    }

    ~FluidSolver()
    {
        // for (int axis = 0; axis < d; axis++)
        //     delete[] velocityField[axis];

        delete[] rhs;
        delete[] pressure_solution;

        // TODO: delete density_field;
    }

    /* Advection */
    void advection(double timestep);

    /* Calculate the RHS of Poisson Equation */
    void calculateRHS();

    /* Projection with CG */
    void projection();

    /* Update velocity with pressure */
    void updateVelocity(T_INDEX &index, double timestep);
    void updateVelocity(double timestep);
    /* Make result visulizable */
    void flip();

    //-----------------------------------------------

    void initialize();
    /* UPDATE */
    void addInflow(const T_INDEX &index, const double density, const TV &velocity);
    void update(double timestep);

    //-----------------------------------------------

    // TODO:
    int index2offset(const T_INDEX &index)
    {

        int os = index[1] * (*grid).counts[0] + index[0];
        if (d == 3)
            os += index[2] * (*grid).counts[0] * (*grid).counts[1];
        return os;
    }

    double nablapOnI(T_INDEX &index, int axis)
    {
        return (pressure_solution[index2offset(index)] -
                pressure_solution[index2offset(Previous_Cell(axis, index))]) *
               (*grid).one_over_dX(axis);
    }
};
