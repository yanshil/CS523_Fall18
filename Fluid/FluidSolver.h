//!#####################################################################
//! \file FluidSolver.h
// https://github.com/OrionQuest/Nova_Examples/blob/master/embedded_deformables/Embedded_Deformables_Example.h
//!#####################################################################
// Class FluidSolver
//######################################################################

#ifndef __FluidSolver__
#define __FluidSolver__

#include "FluidQuantity.h"

namespace Nova
{
template <typename T, int d>
class FluidSolver
{
    using T_INDEX = Vector<int, d>;
    using TV = Vector<T, d>;

  public:
    FluidQuantity<T,d> *velocityField[d];
    FluidQuantity<T,d> *density_field;

    // FluidQuantity pressure;
    Grid<T, d> *grid;

    T density;

    T *divG;
    T *pressure_solution;
    int number_of_ghost_cells;
    int size;
    T_INDEX storing_counts;

    FluidSolver();
    FluidSolver(Grid<T, d> &grid, T density, int number_of_ghost_cells)
    {
        std::cout << "Constructor of FluidSolver" << std::endl;

        for (int axis = 0; axis < d; axis++)
            velocityField[axis] = new FluidQuantity<T,d>(grid, axis, number_of_ghost_cells);

        // TODO
        this->density_field = new FluidQuantity<T,d>(grid, -1, number_of_ghost_cells);
        this->density = density;
        this->grid = &grid;
        this->number_of_ghost_cells = number_of_ghost_cells;
        this->size = grid.counts.Product();

        this->storing_counts = grid.counts;

        // density = FluidQuantity(grid, -1);
        // pressure = FluidQuantity(grid, -1);

        divG = new T[size];
        pressure_solution = new T[size];

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

    T getRGBcolorDensity(T_INDEX &index);

    /* Advection */
    void advection(T timestep);

    /* Calculate the RHS of Poisson Equation */
    void calculateDivergence();

    /* Projection with CG */
    void pressure_solution_Jacobi();
    void projection(int limit);
    void projection(int limit, T timestep = 0.12);
    void Project(int limit);
    /* ================== */
    void SetDivBoundary();
    void SetPressureBoundary();
    void SetVelocityBoundary();
    void SetVelocityBoundary_as0();

    /* Update velocity with pressure */
    void updateVelocity(T timestep);
    /* Make result visulizable */
    void flip();

    //-----------------------------------------------

    void initialize();
    /* UPDATE */
    void addInflow(const T_INDEX &index, const T density, const TV &velocity);

    void update(T timestep);

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

    void projection_INCSAMPLE(int limit, T timestep);
};

} // namespace Nova

#endif