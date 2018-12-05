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

    // Grid with Inertia domain m, n, k
    FluidSimulator_Grid<T, d> *grid;

    // Face Indiicator
    int axis;
    // Grid Cell size
    T hx;

    // Simulation Domain: For velocity[i], simulation_counts(i) = counts(i) + 1
    T_INDEX simulation_counts;
    Range<T, d> simulation_domain;

    // Simulation Domain size.
    int size;

    // 1D Linear Interpolate
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
            // TODO: Fix
            simulation_domain.max_corner(axis) += grid.dX(axis) - grid.dX(axis) * 1e-03;
        }

        hx = grid.hx;

        size = simulation_counts.Product();
        _Phi = new T[size];
        _Phi_new = new T[size];

        // std::cout << "simulation count: " << simulation_counts << std::endl;
        // printf("axis = %d, size = %d\n", axis, size);
        // printf("\n");

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
    void addInflow(const T_INDEX &index, const T input);

    void printPhi()
    {
        // std::cout << "axis = : " << axis << std::endl;

        printf("=========================: m= %d, size = %d\n", simulation_counts[0], size);

        for (int i = 0; i < size; i++)
        {
            printf("%.1f", _Phi[i]);

            if ((i + 1) % simulation_counts[0] == 0)
            {
                printf("\n");
            }
            else
            {
                printf(",");
            }
        }
    }
};

} // namespace Nova

#endif