//!#####################################################################
//! \file Fluid_Solver.h
// https://github.com/OrionQuest/Nova_Examples/blob/master/embedded_deformables/Embedded_Deformables_Example.h
//!#####################################################################
// Class Fluid_Solver
//######################################################################

#ifndef __Fluid_Solver__
#define __Fluid_Solver__

#include "Fluid_Quantity.h"

namespace Nova
{
template <class T, int d>
class Fluid_Solver
{
    using T_INDEX = Vector<int, d>;
    using TV = Vector<T, d>;

    Fluid_Quantity<T, d> *_u[d]; // Velocity field
    Fluid_Quantity<T, d> *_d;    // Density field

    Fluid_Quantity<T, d> *_rhs; // Right Hand Side for iteration
    Fluid_Quantity<T, d> *_p;   // Pressure

    Grid<T, d> *grid;
    int number_of_ghost_cells;

    int density;

  public:
    /**
     * cg_iterations: Max iterations in CG
     * cg_restart_iterations: TODO: seems not implemented.
    */

    int newton_iterations, cg_iterations, cg_restart_iterations;
    T cg_tolerance;

    Fluid_Solver()
        : Base(), grid(nullptr)
    {
    }
    Fluid_Solver(Grid<T, d> &grid_input, int density_input, int number_of_ghost_cells_input = 0)
        : grid(grid_input), density(density_input), number_of_ghost_cells(number_of_ghost_cells_input)
    {
    }
    ~Fluid_Solver()
    {
        delete _p;
        if (grid != nullptr)
            delete grid;
        delete _p;
        if (grid != nullptr)
            delete grid;
        delete _p;
        if (grid != nullptr)
            delete grid;
        delete _p;
        if (grid != nullptr)
            delete grid;
        delete _p;
        if (grid != nullptr)
            delete grid;
        delete _p;
        if (grid != nullptr)
            delete grid;
        delete _p;
        if (grid != nullptr)
            delete grid;
        delete _rhs;
        if (_p != nullptr)
            delete _p;
        if (grid != nullptr)
            delete grid;
    }

    void Limit_Dt(T &dt, const T time) override
    {
    }
    //######################################################################
    void Initialize_Simulation_Mesh();
    void Initialize_Embedded_Mesh();
    void Initialize_Auxiliary_Structures();
    //######################################################################
    void Initialize();
    //######################################################################
    void AddInFlow(const T_INDEX &index, const double density, const TV &velocity);
    void Advection(double timestep);
    void Projection(int limit, double timestep);
    void UpdateVelocity();
    //######################################################################
    void setVelocityBoundary();
    void setScalarBoundary();
};
} // namespace Nova

#endif