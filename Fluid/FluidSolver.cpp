//!#####################################################################
//! \file FluidSolver.cpp
//!#####################################################################
#include "FluidSolver.h"
#include <nova/Tools/Utilities/Range_Iterator.h>

using namespace Nova;
/**
 * calculateRHS
 */
template <typename T, int d>
void FluidSolver<T, d>::calculateRHS()
{
    memset(_rhs, 0, size * sizeof(T));
    for (int idx = 0; idx < size; idx++)
    {
        T_INDEX index = offset2index(idx);

        for (int axis = 0; axis < d; axis++)
        {
            _rhs[idx] -= (_v[axis]->at(grid->Next_Cell(axis, index)) -
                          _v[axis]->at(index)) /
                         hx;
        }
    }
}

/**
 * setBoundaryCondition
 * 
 */
// TODO 3D
template <typename T, int d>
void FluidSolver<T, d>::setBoundaryCondition()
{
    int m = grid->counts[0];
    int n = grid->counts[1];

    for (int ix = 0; ix < m; ix++)
    {
        _v[1]->at(T_INDEX{ix, 0}) = 0.0;
        _v[1]->at(T_INDEX{ix, n}) = 0.0;
    }
    for (int iy = 0; iy < n; iy++)
    {
        _v[0]->at(T_INDEX{m, iy}) = 0.0;
        _v[0]->at(T_INDEX{0, iy}) = 0.0;
    }
}

/**
 * project
 */
template <typename T, int d>
void FluidSolver<T, d>::project_GS(int limit, T timestep, bool output)
{
    T scale = timestep / rho / hx / hx;
    T maxDelta;

    for (int iteration = 0; iteration < limit; iteration++)
    {
        maxDelta = __DBL_MIN__;

        for (int idx = 0; idx < size; idx++)
        {
            T_INDEX index = offset2index(idx);

            T Aii = 0;
            T sum = 0;

            for (int axis = 0; axis < d; axis++)
            {
                T_INDEX p_index = grid->Previous_Cell(axis, index);
                if (grid->Inside_Domain(p_index))
                {
                    Aii += scale;
                    sum -= scale * _p[index2offset(p_index)];
                }

                T_INDEX n_index = grid->Next_Cell(axis, index);
                if (grid->Inside_Domain(n_index))
                {
                    Aii += scale;
                    sum -= scale * _p[index2offset(n_index)];
                }
            }

            T newP = (_rhs[idx] - sum) / Aii;
            maxDelta = std::max(maxDelta, std::fabs(_p[idx] - newP));
            _p[idx] = newP;
        }

        if (maxDelta < 1e-5)
        {
            if (output)
                printf("Converge with %d iteration, with Norm1 = %f\n", iteration, maxDelta);
            return;
        }
    }

    if (output)
        printf("Exceed Limit of %d, with Norm1 = %f\n", limit, maxDelta);
}

/**
 * 
 * ApplyPressure
 */
template <typename T, int d>
void FluidSolver<T, d>::applyPressure(T timestep)
{
    T scale = timestep / rho / hx;

    for (int idx = 0; idx < size; idx++)
    {
        T_INDEX index = offset2index(idx);

        for (int axis = 0; axis < d; axis++)
        {
            _v[axis]->at(index) -= _p[idx] * scale;
            _v[axis]->at(grid->Next_Cell(axis, index)) += _p[idx] * scale;
        }
    }

    setBoundaryCondition();
}

template <typename T, int d>
void FluidSolver<T, d>::advection(T timestep)
{
    _d->advect(timestep, _v);

    for (int axis = 0; axis < d; axis++)
        _v[axis]->advect(timestep, _v);
}

template <typename T, int d>
void FluidSolver<T, d>::flip()
{
    _d->flip();

    for (int axis = 0; axis < d; axis++)
        _v[axis]->flip();
}
/**
 * Update
 */
template <typename T, int d>
void FluidSolver<T, d>::update(T timestep)
{
    calculateRHS();
    project_GS(1000, timestep);
    applyPressure(timestep);

    advection(timestep);
    // _v[1]->printPhi();
    // _v[0]->printPhi();
    flip();
}

/**
 * Addinflow
 */
template <typename T, int d>
void FluidSolver<T, d>::addInflow(const T_INDEX &min_corner, const T_INDEX &max_corner,
                                  int axis, T input_value)
{
    T_INDEX currIndex;
    for (Range_Iterator<d> iterator(Range<int, d>(min_corner, max_corner)); iterator.Valid(); iterator.Next())
    {
        currIndex = T_INDEX() + iterator.Index();

        if (axis == -1)
            _d->addInflow(currIndex, input_value);
        else
            _v[axis]->addInflow(currIndex, input_value);
    }
}
//######################################################################
template class Nova::FluidSolver<float, 2>;
template class Nova::FluidSolver<float, 3>;
#ifdef COMPILE_WITH_DOUBLE_SUPPORT
template class Nova::FluidSolver<double, 2>;
template class Nova::FluidSolver<double, 3>;
#endif