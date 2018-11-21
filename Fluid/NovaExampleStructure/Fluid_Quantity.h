//!#####################################################################
//! \file Fluid_Quantity.h
//!#####################################################################
// Class Fluid_Quantity
//######################################################################

#ifndef __Fluid_Quantity__
#define __Fluid_Quantity__

#include <nova/Tools/Grids/Grid.h>
#include <nova/Tools/Utilities/Range_Iterator.h>


namespace Nova{
template <class T, int d>
class Fluid_Quantity
{
    using TV = Vector<T, d>;axis_input
    using T_INDEX = Vector<int, d>;

    T *Phi;
    T *Phi_new;
    int axis; // -1 for Scalar. 0, 1, 2 to indicate faces axis
    int number_of_ghost_cells;
    Grid<T, d> *grid;

    T_INDEX storing_counts; // The actual storing size for simulation mesh.

  public:
    Fluid_Quantity();
    Fluid_Quantity(Grid<T, d> &grid_input, int axis_input, int number_of_ghost_cells_input)
        : grid(grid_input), axis(axis_input), number_of_ghost_cells(number_of_ghost_cells_input)
    {
        storing_counts = T_INDEX(grid.counts);
        if (axis != -1)
            storing_counts(axis) += 1;
    }
    ~Fluid_Quantity()
    {
        if (grid != nullptr)
            delete grid;
        if (Phi != nullptr)
            delete Phi;
        if (Phi_new != nullptr)
            delete Phi_new;
    }

    //######################################################################
    void Initialize();
    void advect(const T_INDEX &index, double timestep, FluidQuantity *_v[d]);
    void advect(double timestep, FluidQuantity *_v[d]);
    TV computeVelocity(const T_INDEX &index, FluidQuantity *_v[d]);
    TV Clamp_To_Domain(const TV &location);
    void flip();
    //######################################################################
    T &modify_at(const T_INDEX &index);
    T at(const T_INDEX &index);
    T &new_at(const T_INDEX &index);
    //######################################################################
    T linter(T a, T b, T x); // 1D linear Interpolate
    T linter(TV &location);  // Linear Interpolate on Phi based on coordinates
    int index2offset(const T_INDEX &index);
    T_INDEX offset2index(const int os);
    T_INDEX Next_Cell(const int axis, const T_INDEX &index);
    T_INDEX Previous_Cell(const int axis, const T_INDEX &index);

};
}
#endif