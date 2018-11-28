//!#####################################################################
//! \file FluidQuantity.h
//!#####################################################################
// Class FluidQuantity
//######################################################################

#ifndef __FluidQuantity__
#define __FluidQuantity__

#include "FluidSimulator_Grid.h"
#include <nova/Tools/Utilities/Range_Iterator.h>

namespace Nova
{
template <typename T, int d>
class FluidQuantity
{
    using T_INDEX = Vector<int, d>;
    using TV = Vector<T, d>;

  public:
    T *Phi;     // Array of num_cell
    T *Phi_new; // Array of num_cell
    int axis;   // -1 for Scalar. 0, 1, 2 to indicate faces axis
    int number_of_ghost_cells;

    int size_whole_domain;
    int size_interior_domain;
    int size_simulation_domain;
    FluidSimulator_Grid<T, d> *grid;

    T_INDEX interior_domain;
    T_INDEX whole_domain;
    T_INDEX simulation_domain;

    FluidQuantity();
    FluidQuantity(FluidSimulator_Grid<T, d> &grid, int axis, int number_of_ghost_cells);
    ~FluidQuantity();

    /////////////////////////////////////////////////
    /// Fundamental Read & Write
    /////////////////////////////////////////////////
    void fill(T content);
    T at(const T_INDEX &index);
    T &modify_at(const T_INDEX &index);
    T &new_at(const T_INDEX &index);
    T rgb_at(const T_INDEX &index);
    void flip()
    {
        std::swap(Phi, Phi_new);
    }
    /////////////////////////////////////////////////
    /// Fluid Simulate Steps
    /////////////////////////////////////////////////
    /* Compute Velocity */
    TV computeVelocity(const TV &location, FluidQuantity *velocityField[d]);
    /* Advection */
    void advect(const T_INDEX &index, T timestep, FluidQuantity *velocityField[d]);
    void advect(T timestep, FluidQuantity *velocityField[d]);
    /////////////////////////////////////////////////
    /// Auxiliary Function
    /////////////////////////////////////////////////
    TV Clamp_To_Domain(const TV &location);
    T linter(T a, T b, T x);
    T linter(const TV &location);
    bool Inside_Domain(const T_INDEX &index);

    /////////////////////////////////////////////////
    /// Printing Function
    /////////////////////////////////////////////////
    void printPhi()
    {
        std::cout << "Phi____: ";

        for (int i = 0; i < simulation_domain.Product(); i++)
        {
            if (i % simulation_domain[0] == 0)
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

        for (int i = 0; i < simulation_domain.Product(); i++)
        {
            if (i % simulation_domain[0] == 0)
            {
                std::cout << "\n";
            }
            std::cout << Phi_new[i] << ", ";
        }

        std::cout << "\n================================" << std::endl;
    }
};
} // namespace Nova

#endif