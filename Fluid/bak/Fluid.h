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

//------------------------------------------------

// T_INDEX Next_Cell(const int axis, const T_INDEX &index, const int number_of_ghost_cells = 0);
// T_INDEX Previous_Cell(const int axis, const T_INDEX &index, const int number_of_ghost_cells = 0);

class FluidQuantity
{
  public:
    double *Phi;     // Array of num_cell
    double *Phi_new; // Array of num_cell
    int axis;        // -1 for Scalar. 0, 1, 2 to indicate faces axis
    int number_of_ghost_cells;
    Grid<T, d> *grid;

    T_INDEX storing_counts;

    FluidQuantity()
    {
        grid = nullptr;
        std::cout << "Wrong Constructor!" << std::endl;
    }
    FluidQuantity(Grid<T, d> &grid, int axis, int number_of_ghost_cells);
    ~FluidQuantity();

    void fill(double content);

    // ********************************************

    /*!
     * 1D Linear Interpolate etween a and b for x in (0, 1)
     */
    double linter(double a, double b, double x)
    {
        return (1.0 - x) * a + x * b;
    }

    bool Inside_Domain(const T_INDEX &index)
    {
        T_INDEX diff = storing_counts - index;

        for (int i = 0; i < d; i++)
        {
            if (diff(i) < 0)
            {
                return false;
            }
        }
        return true;
    }

    bool Inside_Domain(const TV &location)
    {
        TV tmp_location = TV(location);
        // Clamp to domain if not in domain
        T_INDEX tmp_index = (*grid).Cell(tmp_location, number_of_ghost_cells);

        T_INDEX diff = storing_counts - tmp_index;

        for (int i = 0; i < d; i++)
        {
            if (diff(i) < 0)
            {
                return false;
            }
        }
        return true;
    }

    /*!
    * Auxiliary Function: Get exact offset for Column-based 1D array
    * return (z * xSize * ySize) + (y * xSize) + x;
    */
    // TODO:
    int index2offset(const T_INDEX &index)
    {
        // Becuase index in the grid start from (1,1)...
        T_INDEX tmp_index = index - T_INDEX(1);

        int os = tmp_index[1] * storing_counts[0] + tmp_index[0];
        if (d == 3)
            os += tmp_index[2] * storing_counts[0] * storing_counts[1];
        return os;
    }

    T_INDEX offset2index(const int os)
    {
        // 3D: os = z * m * n + y * m + x
        // 2D: os = y * m + x
        T_INDEX tmp_index = T_INDEX();

        // x <- os mod m
        tmp_index[0] = os % storing_counts[0];

        if (d == 2)
            // y <- (os - x) / m
            tmp_index[1] = (os - tmp_index[0]) / storing_counts[1];
        else
        {
            // y <- (os - x) mod n
            tmp_index[1] = (os - tmp_index[0]) % storing_counts[1];

            // z <- (os - x - y * m) / (m*n)
            tmp_index[2] = (os - tmp_index[0] - tmp_index[1] * storing_counts[0]) / storing_counts[0] / storing_counts[1];
        }

        // Becuase index in the grid start from (1,1)...
        tmp_index += T_INDEX(1);
        return tmp_index;
    }

    /*!
     * Linear Interpolator for TV{i, j, k} on grid
     * Coordinates will be clamped to lie in simulation domain
    */
    double linter(TV &location);

    //----------Auxiliary Function--------------------
    void printPhi()
    {
        std::cout << "Phi____: ";

        for (int i = 0; i < storing_counts.Product(); i++)
        {
            if (i % storing_counts[0] == 0)
            {
                std::cout << "\n";
            }
            std::cout << Phi[i] << ", ";
        }

        std::cout << "\n================================" << std::endl;
    }

    void printPhi_new()
    {
        std::cout << "Phi_new: ";

        for (int i = 0; i < storing_counts.Product(); i++)
        {
            if (i % storing_counts[0] == 0)
            {
                std::cout << "\n";
            }
            std::cout << Phi_new[i] << ", ";
        }

        std::cout << "\n================================" << std::endl;
    }

    T_INDEX Next_Cell(const int axis, const T_INDEX &index)
    {
        T_INDEX shifted_index(index);
        shifted_index(axis) += 1;

        return shifted_index;
    }
    // TODO
    T_INDEX Previous_Cell(const int axis, const T_INDEX &index)
    {
        T_INDEX shifted_index(index);
        shifted_index(axis) -= 1;

        return shifted_index;
    }

    /*!
     * Read and Write access in the quantity field for Phi
     */
    double &modify_at(const T_INDEX &index)
    {
        if (!Inside_Domain(index))
        {
            std::cout << "index = " << index << std::endl;
            // Raise exception
            throw std::runtime_error("Try to write at an Out_of_domain area");
        }
        // if (!(*grid).Inside_Domain(index))
        // {
        //     std::cout << "index = " << index << std::endl;
        //     // Raise exception
        //     throw std::runtime_error("Try to write at an Out_of_domain area");
        // }
        return Phi[index2offset(index)];
    }

    /*!
     * Read access for the quantity field of Phi
     */
    double at(const T_INDEX &index)
    {
        // e.g. For n = 4
        // Density: 4 * 4 * 4;      velocityU: 5 * 4 * 4;   velocity V: 4 * 5 * 4;  velocity W: 4 * 4 * 5
        // this->Inside_Domain are checking if is in (1,1,1) to (5, 4, 4) [For U]
        // (*grid).Inside_Domain(index) is checking if is in (1,1,1) to (4,4,4)
        // (*grid).Inside_Domain(index, num_of_ghost_cell = 1) is checking if is in (0,0,0) to (5,5,5)

        //----------------- Velocity /Scalar Field Domain (d/u/v domain) ------------------
        if (!Inside_Domain(index))
        {
            TV location = (*grid).Center(index);
            T_INDEX clamped_index = (*grid).Clamp_To_Cell(location, number_of_ghost_cells);
            return Phi[index2offset(clamped_index)];
        }

        //------------------ Grid Domain (grid domain (ghost_cell = 0)) ----------------------------
        // if (!(*grid).Inside_Domain(index))
        // {
        //     TV location = (*grid).Center(index);
        //     T_INDEX clamped_index = (*grid).Clamp_To_Cell(location, number_of_ghost_cells);

        //     return Phi[index2offset(clamped_index)];
        // }

        return Phi[index2offset(index)];
    }

    /*!
     * RGB color range 0-1 access in the
     */
    double rgb_at(const T_INDEX &index)
    {
        return std::max(std::min(1.0 - at(index), 1.0), 0.0);
    }

    /*!
     * Read & Write access in the quantity field for Phi_new
     */
    double &new_at(const T_INDEX &index)
    {
        if (!Inside_Domain(index))
        {
            // Raise exception
            throw std::runtime_error("Try to write at an Out_of_domain area");
        }
        // if (!(*grid).Inside_Domain(index))
        // {
        //     // Raise exception
        //     throw std::runtime_error("Try to write at an Out_of_domain area");
        // }
        return Phi_new[index2offset(index)];
    }

    /* Compute Velocity */
    TV computeVelocity(const T_INDEX &index, FluidQuantity *velocityField[d]);
    TV Clamp_To_Domain(const TV &location);
    // TV traceBack(const TV &location, double timestep, TV &velocity);

    void flip()
    {
        std::swap(Phi, Phi_new);
    }

    /* Advection */
    void advect(const T_INDEX &index, double timestep, FluidQuantity *velocityField[d]);
    void advect(double timestep, FluidQuantity *velocityField[d]);

    /* */
    // void updateVelocity(const T_INDEX &index, double timestep, double np)
    // {
    //     new_at(index) -= np * timestep;
    // }
};

class FluidSolver
{
    FluidQuantity *velocityField[d];
    FluidQuantity *density_field;

    // FluidQuantity pressure;
    Grid<T, d> *grid;

    double density;

    double *divG;
    double *pressure_solution;
    int number_of_ghost_cells;
    int size;
    T_INDEX storing_counts;

  public:
    FluidSolver();
    FluidSolver(Grid<T, d> &grid, double density, int number_of_ghost_cells)
    {
        std::cout << "Constructor of FluidSolver" << std::endl;

        for (int axis = 0; axis < d; axis++)
            velocityField[axis] = new FluidQuantity(grid, axis, number_of_ghost_cells);

        // TODO
        this->density_field = new FluidQuantity(grid, -1, number_of_ghost_cells);
        this->density = density;
        this->grid = &grid;
        this->number_of_ghost_cells = number_of_ghost_cells;
        this->size = grid.counts.Product();

        this->storing_counts = grid.counts;

        // density = FluidQuantity(grid, -1);
        // pressure = FluidQuantity(grid, -1);

        divG = new double[size];
        pressure_solution = new double[size];

        for (int i = 0; i < size; i++)
        {
            divG[i] = 0;
            pressure_solution[i] = 0;
        }
    }

    ~FluidSolver()
    {
        for (int axis = 0; axis < d; axis++)
            delete velocityField[axis];

        delete density_field;

        delete divG;
        delete pressure_solution;
    }

    double getRGBcolorDensity(T_INDEX &index);

    /* Advection */
    void advection(double timestep);

    /* Calculate the RHS of Poisson Equation */
    void calculateDivergence();

    /* Projection with CG */
    void pressure_solution_Jacobi();
    void projection(int limit);
    void projection(int limit, double timestep = 0.12);
    void Project(int limit);
    /* ================== */
    void SetDivBoundary();
    void SetPressureBoundary();
    void SetVelocityBoundary();
    void SetVelocityBoundary_as0();

    /* Update velocity with pressure */
    void updateVelocity(double timestep);
    /* Make result visulizable */
    void flip();

    //-----------------------------------------------

    void initialize();
    /* UPDATE */
    void addInflow(const T_INDEX &index, const double density, const TV &velocity);

    void update(double timestep);

    //-----------------------------------------------

    //----------Auxiliary Function--------------------
    // TODO
    T_INDEX Next_Cell(const int axis, const T_INDEX &index)
    {
        T_INDEX shifted_index(index);
        shifted_index(axis) += 1;

        return shifted_index;
    }
    // TODO
    T_INDEX Previous_Cell(const int axis, const T_INDEX &index)
    {
        T_INDEX shifted_index(index);
        shifted_index(axis) -= 1;

        return shifted_index;
    }

    // TODO:
    int index2offset(const T_INDEX &index)
    {
        // Becuase index in the grid start from (1,1)...
        T_INDEX tmp_index = index - T_INDEX(1);

        int os = tmp_index[1] * storing_counts[0] + tmp_index[0];
        if (d == 3)
            os += tmp_index[2] * storing_counts[0] * storing_counts[1];
        return os;
    }

    T_INDEX offset2index(const int os)
    {
        // 3D: os = z * m * n + y * m + x
        // 2D: os = y * m + x
        T_INDEX tmp_index = T_INDEX();

        // x <- os mod m
        tmp_index[0] = os % storing_counts[0];

        if (d == 2)
            // y <- (os - x) / m
            tmp_index[1] = (os - tmp_index[0]) / storing_counts[1];
        else
        {
            // y <- (os - x) mod n
            tmp_index[1] = (os - tmp_index[0]) % storing_counts[1];

            // z <- (os - x - y * m) / (m*n)
            tmp_index[2] = (os - tmp_index[0] - tmp_index[1] * storing_counts[0]) / storing_counts[0] / storing_counts[1];
        }

        // Becuase index in the grid start from (1,1)...
        tmp_index += T_INDEX(1);
        return tmp_index;
    }

    void projection_INCSAMPLE(int limit, double timestep);
};
