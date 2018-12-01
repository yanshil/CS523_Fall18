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
  private:
    using TV = Vector<T, d>;
    using T_INDEX = Vector<int, d>;

    T *_Phi;
    T *_Phi_new;

    FluidSimulator_Grid<T, d> *grid; // with grid size

    T hx;

    int axis; // indicate face offset

    T_INDEX simulation_counts;
    Range<T, d> simulation_domain;
    int size;

    T linp(T a, T b, T delta) const
    {
        return a * (1 - delta) + b * delta;
    }

  public:
    FluidQuantity()
    {
        std::cout << "??" << std::endl;
    }
    FluidQuantity(FluidSimulator_Grid<T, d> &grid, int axis)
        : grid(&grid), axis(axis)
    {
        // Simulation Domain
        simulation_counts = T_INDEX(grid.counts);
        simulation_domain = Range<T, d>(grid.domain);
        if (axis != -1)
        {
            simulation_counts(axis) += 1;
            simulation_domain.max_corner(axis) += grid.dX(axis);
        }

        size = simulation_counts.Product();
        _Phi = new T[size];
        _Phi_new = new T[size];

        memset(_Phi, 0, size * sizeof(T));
    }
    ~FluidQuantity()
    {
        delete[] _Phi;
        delete[] _Phi_new;
    }

    int index2offset(const T_INDEX &index) const
    {
        return grid->index2offset(index, simulation_counts);
    }
    T_INDEX offset2index(const int os) const
    {
        return grid->offset2index(os, simulation_counts);
    }

    void flip()
    {
        std::swap(_Phi, _Phi_new);
    }

    const T *Phi() const
    {
        return _Phi;
    }

    T at(const T_INDEX &index) const
    {
        return _Phi[index2offset(index)];
    }
    T &at(const T_INDEX &index)
    {
        return _Phi[index2offset(index)];
    }

    T linp(TV location) const;
    void advect(T timestep, FluidQuantity *_v[d]);
    void addInflow(const T_INDEX &min_corner, const T_INDEX &max_corner, const T input);

    void printPhi()
    {
        std::cout << "axis = : " << axis;

        for (int i = 0; i < size; i++)
        {
            if (i % simulation_counts[0] == 0)
            {
                std::cout << "\n";
            }
            std::cout << _Phi[i] << ", ";
        }

        std::cout << "\n================================" << std::endl;
    }
};

} // namespace Nova

#endif