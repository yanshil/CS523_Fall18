#include "FluidQuantity.h"
#include <nova/Tools/Utilities/Range_Iterator.h>

namespace Nova
{
template <class T, int d>
class FluidSolver
{
    using TV = Vector<T, d>;
    using T_INDEX = Vector<int, d>;

    // Density Field
    FluidQuantity<T, d> *_d;
    // Velocity Field
    FluidQuantity<T, d> *_v[d];

    FluidSimulator_Grid<T, d> *grid;

    int size;
    T hx;
    T rho;

    T *_rhs;
    T *_p;

    int index2offset(const T_INDEX &index)
    {
        return grid->index2offset(index, grid->counts);
    }

    T_INDEX offset2index(const int os)
    {
        return grid->offset2index(os, grid->counts);
    }

  public:
    FluidSolver(FluidSimulator_Grid<T, d> &grid, T rho) : grid(&grid), rho(rho)
    {
        this->hx = grid.hx;
        this->size = grid.counts.Product();

        _d = new FluidQuantity<T, d>(grid, -1);

        for (int axis = 0; axis < d; axis++)
            _v[axis] = new FluidQuantity<T, d>(grid, axis);

        _rhs = new T[size];
        _p = new T[size];
    }

    ~FluidSolver()
    {
        delete _d;
        for (int axis = 0; axis < d; axis++)
            delete _v[axis];

        delete[] _rhs;
        delete[] _p;
    }

    const T *rhs() const
    {
        return _rhs;
    }

    void calculateRHS()
    {
        memset(_rhs, 0, size * sizeof(T));
        for (int idx = 0; idx < size; idx++)
        {
            T_INDEX index = offset2index(idx);

            for (int axis = 0; axis < d; axis++)
                _rhs[idx] -= (_v[axis]->at(grid->Next_Cell(axis, index)) - _v[axis]->at(index)) / hx;
        }
    }

    void project_GS(int limit, T timestep, bool output = true)
    {
        T scale = timestep / rho / hx / hx;

        T maxDelta;
        for (int iteration = 0; iteration < limit; iteration++)
        {
            maxDelta = __DBL_MIN__;

            for (int idx = 0; idx < size; idx++)
            {
                T_INDEX index = offset2index(idx);
                T Aii = 0, sum = 0;

                for (int axis = 0; axis < d; axis++)
                {
                    // Previous u
                    T_INDEX p_index = grid->Previous_Cell(axis, index);

                    if (grid->Inside_Domain(p_index))
                    {
                        Aii += scale;
                        sum -= scale * _p[index2offset(p_index)];
                    }

                    // Next u
                    T_INDEX n_index = grid->Next_Cell(axis, index);

                    if (grid->Inside_Domain(n_index))
                    {
                        Aii += scale;
                        sum -= scale * _p[index2offset(n_index)];
                    }
                }

                T newP = (_rhs[idx] - sum) / Aii;

                maxDelta = std::max(maxDelta, fabs(_p[idx] - newP));

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

    void setBoundaryCondition()
    {
        for (int axis = 0; axis < d; axis++)
        {
            _v[axis]->setBoundaryValue();
        }
    }

    void applyPressure(T timestep)
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

    void initialize()
    {
        for (int axis = 0; axis < d; axis++)
        {
            _v[axis]->calculateBCFlag();
        }
    }

    void update(T timestep)
    {
        // Projection
        calculateRHS();

        project_GS(1000, timestep, false);
        applyPressure(timestep);

        //Advection
        _d->advect(timestep, _v);
        
        for(int axis = 0; axis < d; axis++)
            _v[axis]->advect(timestep, _v);
        

        // Flip
        _d->flip();

        for(int axis = 0; axis < d; axis++)
            _v[axis]->flip();
        

    }

    void addInflow(int ix0, int iy0, int ix1, int iy1, int axis, T value)
    {
        for (int x = std::max(ix0, 0); x < std::min(ix1, grid->counts[0]); x++)
            for (int y = std::max(iy0, 0); y < std::min(iy1, grid->counts[1]); y++)
            {
                int idx = y * grid->counts[0] + x;
                T_INDEX index = offset2index(idx);
                std::cout << "index = " << index << std::endl;

                if (axis == -1)
                    _d->addInflow(index, value);
                else
                {
                    _v[axis]->addInflow(index, value);
                }
            }
    }

    // ?????? Why Not Work ???
    void addInflow(const T_INDEX &min_corner, const T_INDEX &max_corner, int axis, T value)
    {
        T_INDEX index;
        for (Range_Iterator<d> iterator(Range<int, d>(min_corner, max_corner)); iterator.Valid(); iterator.Next())
        {
            index = T_INDEX() + iterator.Index();
            std::cout << "index = " << index << std::endl;

            if (axis == -1)
                _d->addInflow(index, value);
            else
            {
                _v[axis]->addInflow(index, value);
            }
        }
    }

    /* Convert fluid density to RGB color scaled in 0 - 1 */
    T toRGB(const T_INDEX &index)
    {
        return std::max(std::min(1.0 - _d->at(index), 1.0), 0.0);
    }

    //===========================================================

    // Projection output as a tent-like plot
    void test_projection_tent(T timestep)
    {
        // Projection
        // calculateRHS();
        memset(_rhs, 0, size * sizeof(T));
        _rhs[7740] = 1;
        project_GS(1000, timestep, false);
        applyPressure(timestep);
    }
};
} // namespace Nova