//!#####################################################################
//! \file FluidSolver.cpp
//!#####################################################################
#include "FluidSolver.h"
using namespace Nova;
//######################################################################
template <typename T, int d>
FluidSolver<T, d>::FluidSolver(Grid<T, d> &grid, T density, int number_of_ghost_cells)
{
    std::cout << "Constructor of FluidSolver" << std::endl;

    for (int axis = 0; axis < d; axis++)
        velocityField[axis] = new FluidQuantity<T, d>(grid, axis, number_of_ghost_cells);

    // TODO
    this->density_field = new FluidQuantity<T, d>(grid, -1, number_of_ghost_cells);
    this->density = density;
    this->grid = &grid;
    this->number_of_ghost_cells = number_of_ghost_cells;
    this->size = grid.counts.Product();

    this->storing_counts = grid.counts;

    // density = FluidQuantity(grid, -1);
    // pressure = FluidQuantity(grid, -1);

    _rhs = new T[size];
    _p = new T[size];

    for (int i = 0; i < size; i++)
    {
        _rhs[i] = 0;
        _p[i] = 0;
    }
}

template <typename T, int d>
FluidSolver<T, d>::~FluidSolver()
{
    for (int axis = 0; axis < d; axis++)
        delete velocityField[axis];

    delete density_field;

    delete _rhs;
    delete _p;
}
//######################################################################
template <typename T, int d>
T FluidSolver<T, d>::getRGBcolorDensity(T_INDEX &index)
{
    return (*density_field).rgb_at(index);
}

//######################################################################
// TODO
//######################################################################
template <typename T, int d>
void FluidSolver<T, d>::initialize()
{

    // for (int i = 0; i < d; i++)
    // {
    //     (*velocityField[i]).fill(0.5);
    // }
    (*velocityField[0]).fill(0);
    (*velocityField[1]).fill(0.1);
    // (*velocityField[2]).fill(0.5);

    (*density_field).fill(0);
}
//######################################################################
// TODO
//######################################################################
template <typename T, int d>
void FluidSolver<T, d>::advection(T timestep)
{
    (*density_field).advect(timestep, velocityField);

    for (int i = 0; i < d; i++)
    {
        (*velocityField[i]).advect(timestep, velocityField);
    }
}

//----------Auxiliary Function--------------------
// TODO
template <typename T, int d>
Vector<int, d> FluidSolver<T, d>::Next_Cell(const int axis, const T_INDEX &index)
{
    T_INDEX shifted_index(index);
    shifted_index(axis) += 1;

    return shifted_index;
}
// TODO
template <typename T, int d>
Vector<int, d> FluidSolver<T, d>::Previous_Cell(const int axis, const T_INDEX &index)
{
    T_INDEX shifted_index(index);
    shifted_index(axis) -= 1;

    return shifted_index;
}

// TODO:
template <typename T, int d>
int FluidSolver<T, d>::index2offset(const T_INDEX &index)
{
    // Becuase index in the grid start from (1,1)...
    T_INDEX tmp_index = index - T_INDEX(1);

    int os = tmp_index[1] * storing_counts[0] + tmp_index[0];
    if (d == 3)
        os += tmp_index[2] * storing_counts[0] * storing_counts[1];
    return os;
}

template <typename T, int d>
Vector<int, d> FluidSolver<T, d>::offset2index(const int os)
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

// TODO
template <typename T, int d>
void copy_(T src[], T dst[], int size)
{
    for (int i = 0; i < size; i++)
    {
        dst[i] = src[i];
    }
}

// TODO
template <typename T, int d>
void printArray(T arr[], int size)
{
    std::cout << "----\n";
    for (int i = 0; i < size; i++)
    {

        if (i % 16 == 0)
        {
            std::cout << "\n";
        }

        std::cout << arr[i] << ", ";
    }
    std::cout << std::endl;
}

template <typename T, int d>
void FluidSolver<T, d>::calculateDivergence()
{

    T_INDEX currIndex;
    for (Range_Iterator<d> iterator(Range<int, d>(T_INDEX(1), storing_counts)); iterator.Valid(); iterator.Next())
    {
        currIndex = T_INDEX() + iterator.Index();
        int idx = index2offset(currIndex);
        _rhs[idx] = 0;

        for (int axis = 0; axis < d; axis++)
        {
            T_INDEX n_index = Next_Cell(axis, currIndex);
            if ((*grid).Inside_Domain(n_index))
            {
                _rhs[idx] += ((*velocityField[axis]).at(n_index) - (*velocityField[axis]).at(currIndex)) * (*grid).one_over_dX[axis];
                // _rhs[idx] += ((*velocityField[axis]).new_at(n_index) - (*velocityField[axis]).new_at(currIndex)) * (*grid).one_over_dX[axis];
            }
        }
    }
}
template <typename T, int d>
T calculate1Norm(T a[], T b[], int size)
{
    T max = __DBL_MIN__;
    for (int i = 0; i < size; i++)
        max = std::max(std::abs(a[i] - b[i]), max);

    return max;
}
template <typename T, int d>
void FluidSolver<T, d>::Project(int limit)
{
    T scale = (*grid).one_over_dX[0] * (*grid).one_over_dX[0];
    double maxDelta;

    for (int iteration = 0; iteration < limit; iteration++)
    {
        maxDelta = __DBL_MIN__;
        T_INDEX index;
        for (Range_Iterator<d> iterator(Range<int, d>(T_INDEX(1), (*grid).Number_Of_Cells())); iterator.Valid(); iterator.Next())
        {
            index = T_INDEX() + iterator.Index();

            T Aii = 0;
            T AijSum = 0;
            int os = index2offset(index);

            for (int axis = 0; axis < d; axis++)
            {
                T_INDEX p_index = Previous_Cell(axis, index);
                T_INDEX n_index = Next_Cell(axis, index);
                if ((*grid).Inside_Domain(p_index))
                {
                    Aii -= scale;
                    AijSum += scale * _p[index2offset(p_index)];
                }
                if ((*grid).Inside_Domain(n_index))
                {
                    Aii -= scale;
                    AijSum += scale * _p[index2offset(n_index)];
                }
            }

            T newP = (_rhs[os] - AijSum) / Aii;

            maxDelta = std::max(maxDelta, fabs(_p[os] - newP));

            _p[os] = newP;
        }

        if (maxDelta < 1e-5)
        {
            printf("Exiting solver after %d iterations, maximum change is %f\n", iteration, maxDelta);
            return;
        }
    }
    printf("Exceeded budget of %d iterations, maximum change was %f\n", limit, maxDelta);
}
template <typename T, int d>
void FluidSolver<T, d>::updateVelocity(T timestep)
{
    // Update Velocity
    T_INDEX index;
    for (Range_Iterator<d> iterator(Range<int, d>(T_INDEX(1), (*grid).Number_Of_Cells())); iterator.Valid(); iterator.Next())
    {
        index = T_INDEX() + iterator.Index();
        // _p[index2offset(index)] = size * _rhs[index2offset(index)];

        for (int axis = 0; axis < d; axis++)
        {
            if ((*grid).Inside_Domain(Previous_Cell(axis, index)))
            {
                (*velocityField[axis]).new_at(index) -= timestep * (_p[index2offset(index)] - _p[index2offset(Previous_Cell(axis, index))]) *
                                                        (*grid).one_over_dX(axis);
            }
        }
    }
}
template <typename T, int d>
void FluidSolver<T, d>::flip()
{
    (*density_field).flip();

    for (int i = 0; i < d; i++)
    {
        (*velocityField[i]).flip();
    }
}
template <typename T, int d>
void FluidSolver<T, d>::addInflow(const T_INDEX &index, const T density, const TV &velocity)
{
    (*density_field).modify_at(index) = density;

    for (int i = 0; i < d; i++)
    {
        (*velocityField[i]).modify_at(index) = velocity[i];
    }
}
template <typename T, int d>
void FluidSolver<T, d>::SetPressureBoundary()
{
    //1.m -> m,n
    Range<int, d> vTop(T_INDEX{1, (*grid).counts[0]}, (*grid).counts);
    // 1,1 -> 1, m
    Range<int, d> uLeft(T_INDEX(1), T_INDEX{1, (*grid).counts[0]});
    // n,1 -> m,n
    Range<int, d> uRight(T_INDEX{(*grid).counts[1], 1}, (*grid).counts);

    T_INDEX currIndex;
    for (Range_Iterator<d> iterator(vTop); iterator.Valid(); iterator.Next())
    {
        currIndex = T_INDEX() + iterator.Index();
        _p[index2offset(currIndex)] = _p[index2offset(Previous_Cell(1, currIndex))];
    }

    for (Range_Iterator<d> iterator(uLeft); iterator.Valid(); iterator.Next())
    {
        currIndex = T_INDEX() + iterator.Index();
        _p[index2offset(currIndex)] = _p[index2offset(Next_Cell(0, currIndex))];
    }

    for (Range_Iterator<d> iterator(uRight); iterator.Valid(); iterator.Next())
    {
        currIndex = T_INDEX() + iterator.Index();
        _p[index2offset(currIndex)] = _p[index2offset(Previous_Cell(0, currIndex))];
    }
}
template <typename T, int d>
void FluidSolver<T, d>::SetVelocityBoundary()
{
    //1.m -> m,n
    Range<int, d> vTop(T_INDEX{1, (*grid).counts[0]}, (*grid).counts);
    // 1,1 -> 1, m
    Range<int, d> uLeft(T_INDEX(1), T_INDEX{1, (*grid).counts[0]});
    // n,1 -> m,n
    Range<int, d> uRight(T_INDEX{(*grid).counts[1], 1}, (*grid).counts);

    T_INDEX currIndex;
    for (Range_Iterator<d> iterator(vTop); iterator.Valid(); iterator.Next())
    {
        currIndex = T_INDEX() + iterator.Index();
        (*velocityField[0]).new_at(currIndex) = (*velocityField[0]).new_at(Previous_Cell(1, currIndex));
    }

    for (Range_Iterator<d> iterator(uLeft); iterator.Valid(); iterator.Next())
    {
        currIndex = T_INDEX() + iterator.Index();
        (*velocityField[1]).new_at(currIndex) = (*velocityField[1]).new_at(Next_Cell(0, currIndex));
    }

    for (Range_Iterator<d> iterator(uRight); iterator.Valid(); iterator.Next())
    {
        currIndex = T_INDEX() + iterator.Index();
        (*velocityField[1]).new_at(currIndex) = (*velocityField[1]).new_at(Previous_Cell(0, currIndex));
    }
}
template <typename T, int d>
void FluidSolver<T, d>::SetVelocityBoundary_as0()
{
    //1.m -> m,n
    Range<int, d> vTop(T_INDEX{1, (*grid).counts[0]}, (*grid).counts);
    // 1,1 -> 1, m
    Range<int, d> uLeft(T_INDEX(1), T_INDEX{1, (*grid).counts[0]});
    // n,1 -> m,n
    Range<int, d> uRight(T_INDEX{(*grid).counts[1], 1}, (*grid).counts);

    T_INDEX currIndex;
    for (Range_Iterator<d> iterator(vTop); iterator.Valid(); iterator.Next())
    {
        currIndex = T_INDEX() + iterator.Index();
        (*velocityField[1]).new_at(currIndex) = 0;
    }

    for (Range_Iterator<d> iterator(uLeft); iterator.Valid(); iterator.Next())
    {
        currIndex = T_INDEX() + iterator.Index();
        (*velocityField[0]).new_at(currIndex) = 0;
    }

    for (Range_Iterator<d> iterator(uRight); iterator.Valid(); iterator.Next())
    {
        currIndex = T_INDEX() + iterator.Index();
        (*velocityField[0]).new_at(currIndex) = 0;
    }
}
template <typename T, int d>
void FluidSolver<T, d>::SetDivBoundary()
{
    //1.m -> m,n
    Range<int, d> vTop(T_INDEX{1, (*grid).counts[0]}, (*grid).counts);
    // 1,1 -> 1, m
    Range<int, d> uLeft(T_INDEX(1), T_INDEX{1, (*grid).counts[0]});
    // n,1 -> m,n
    Range<int, d> uRight(T_INDEX{(*grid).counts[1], 1}, (*grid).counts);

    T_INDEX currIndex;
    for (Range_Iterator<d> iterator(vTop); iterator.Valid(); iterator.Next())
    {
        currIndex = T_INDEX() + iterator.Index();
        _rhs[index2offset(currIndex)] = _rhs[index2offset(Previous_Cell(1, currIndex))];
    }

    for (Range_Iterator<d> iterator(uLeft); iterator.Valid(); iterator.Next())
    {
        currIndex = T_INDEX() + iterator.Index();
        _rhs[index2offset(currIndex)] = _rhs[index2offset(Next_Cell(0, currIndex))];
    }

    for (Range_Iterator<d> iterator(uRight); iterator.Valid(); iterator.Next())
    {
        currIndex = T_INDEX() + iterator.Index();
        _rhs[index2offset(currIndex)] = _rhs[index2offset(Previous_Cell(0, currIndex))];
    }
}
template <typename T, int d>
void FluidSolver<T, d>::update(T timestep)
{
    //  std::cout<< "Running" <<std::endl;
    calculateDivergence();
    advection(timestep);
    // SetVelocityBoundary_as0();

    Project(1000);
    updateVelocity(timestep);

    flip();
    // SetVelocityBoundary_as0();

    for (int i = 0; i < size; i++)
    {
        _p[i] = 0;
    }
}
//######################################################################
template class Nova::FluidSolver<float, 2>;
template class Nova::FluidSolver<float, 3>;
#ifdef COMPILE_WITH_DOUBLE_SUPPORT
template class Nova::FluidSolver<double, 2>;
template class Nova::FluidSolver<double, 3>;
#endif