//!#####################################################################
//! \file FluidSolver.cpp
//!#####################################################################
#include "FluidSolver.h"
using namespace Nova;
//######################################################################
template <typename T, int d>
FluidSolver<T, d>::FluidSolver(FluidSimulator_Grid<T, d> &grid, T density, int number_of_ghost_cells)
    : grid(&grid), density(density), number_of_ghost_cells(number_of_ghost_cells)
{
    std::cout << "Constructor of FluidSolver" << std::endl;

    for (int axis = 0; axis < d; axis++)
        _v[axis] = new FluidQuantity<T, d>(grid, axis, number_of_ghost_cells);

    this->_d = new FluidQuantity<T, d>(grid, -1, number_of_ghost_cells);

    this->interior_domain = grid.counts;
    this->whole_domain = grid.counts + T_INDEX(number_of_ghost_cells * 2);

    this->size_interior_domain = interior_domain.Product();
    this->size_whole_domain = whole_domain.Product();

    _rhs = new T[size_interior_domain];
    _p = new T[size_interior_domain];

    // _BCflag = new Boundary_Condition[size_whole_domain];

    _BCflag = new int[size_whole_domain];

    for (int i = 0; i < size_interior_domain; i++)
    {
        _rhs[i] = 0;
        _p[i] = 0;
    }
}

template <typename T, int d>
FluidSolver<T, d>::~FluidSolver()
{
    for (int axis = 0; axis < d; axis++)
        delete _v[axis];

    delete _d;

    delete _rhs;
    delete _p;
    delete _BCflag;
}
//######################################################################
template <typename T, int d>
T FluidSolver<T, d>::getRGBcolorDensity(T_INDEX &index)
{
    return (*_d).rgb_at(index);
}

//######################################################################
// TODO
//######################################################################
template <typename T, int d>
void FluidSolver<T, d>::initialize()
{

    (*_d).fill(0);

    // for (int i = 0; i < d; i++)
    // {
    //     (*_v[i]).fill(0.5);
    // }
    (*_v[0]).fill(0);
    (*_v[1]).fill(0.1);

    // // Circle Test
    // T_INDEX index;
    // for (Range_Iterator<d> iterator(Range<int, d>(T_INDEX(1), interior_domain)); iterator.Valid(); iterator.Next())
    // {
    //     index = T_INDEX() + iterator.Index();
    //     int idx = grid->index2offset(index);

    //     for (int axis = 0; axis < d; axis++)
    //     {
    //         // std::cout << "index = " << index << std::endl;
    //         double sum = (index(0) - 16) * (index(0) - 16) + (index(1) - 16) * (index(1) - 16);
    //         // double sum = 1.0;
    //         (*_v[axis]).modify_at(index) = (index(1 - axis) - 16) / sum;
    //     }
    // }
}
//######################################################################
// TODO
//######################################################################
template <typename T, int d>
void FluidSolver<T, d>::advection(T timestep)
{
    (*_d).advect(timestep, _v);

    for (int i = 0; i < d; i++)
    {
        (*_v[i]).advect(timestep, _v);
    }
}

//----------Auxiliary Function--------------------

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
    for (Range_Iterator<d> iterator(Range<int, d>(T_INDEX(1), interior_domain)); iterator.Valid(); iterator.Next())
    {
        currIndex = T_INDEX() + iterator.Index();
        int idx = grid->index2offset(currIndex);
        _rhs[idx] = 0;

        for (int axis = 0; axis < d; axis++)
        {
            T_INDEX n_index = grid->Next_Cell(axis, currIndex);
            // Next Cell Inside whole_domain
            if (grid->Inside_Domain(n_index, number_of_ghost_cells))
            {
                _rhs[idx] += ((*_v[axis]).at(n_index) - (*_v[axis]).at(currIndex)) * (*grid).one_over_dX[axis];
                // _rhs[idx] += ((*_v[axis]).new_at(n_index) - (*_v[axis]).new_at(currIndex)) * (*grid).one_over_dX[axis];
            }
            else
            {
                // Next Cell not in Inertia domain
                std::cout << "should not happen?" << std::endl;
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
    // TODO: Modify 
    calculateDivergence();

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
            int os = grid->index2offset(index);

            for (int axis = 0; axis < d; axis++)
            {
                T_INDEX p_index = grid->Previous_Cell(axis, index);
                T_INDEX n_index = grid->Next_Cell(axis, index);
                if ((*grid).Inside_Domain(p_index))
                {
                    Aii -= scale;
                    AijSum += scale * _p[grid->index2offset(p_index)];
                }
                if ((*grid).Inside_Domain(n_index))
                {
                    Aii -= scale;
                    AijSum += scale * _p[grid->index2offset(n_index)];
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
        // _p[grid->index2offset(index)] = size_whole_domain * _rhs[grid->index2offset(index)];

        for (int axis = 0; axis < d; axis++)
        {
            if ((*grid).Inside_Domain(grid->Previous_Cell(axis, index)))
            {
                (*_v[axis]).new_at(index) -= timestep * (_p[grid->index2offset(index)] - _p[grid->index2offset(grid->Previous_Cell(axis, index))]) *
                                             (*grid).one_over_dX(axis);
            }
        }
    }
}

template <typename T, int d>
void FluidSolver<T, d>::flip()
{
    (*_d).flip();

    for (int i = 0; i < d; i++)
    {
        (*_v[i]).flip();
    }
}
template <typename T, int d>
void FluidSolver<T, d>::addInflow(const T_INDEX &index, const T density)
{
    (*_d).modify_at(index) = density;
}

template <typename T, int d>
void FluidSolver<T, d>::setBoundaryCondition()
{
    T_INDEX max_corner = grid->Cell_Indices(number_of_ghost_cells).max_corner;
    int m = max_corner(0);
    int n = max_corner(1);

    // 0, 33 -> 33, 33
    Range<int, d> vTop(T_INDEX{0, n}, max_corner);
    // 0,0 -> 33, 0
    Range<int, d> uLeft(T_INDEX(0), T_INDEX{0, n});
    // 33,0 -> 33, 33
    Range<int, d> uRight(T_INDEX{m, 0}, max_corner);
    // 0,0 -> m,n
    Range<int, d> vBottom(T_INDEX(0), T_INDEX{m, 0});

    T_INDEX currIndex;
    // On Top: u -> ghost; v <- 0
    for (Range_Iterator<d> iterator(vTop); iterator.Valid(); iterator.Next())
    {
        currIndex = T_INDEX() + iterator.Index();
        (*_v[0]).new_at(currIndex) = (*_v[0]).new_at(grid->Previous_Cell(1, currIndex));
        _d->new_at(currIndex) = _d->new_at(grid->Previous_Cell(1, currIndex));
        (*_v[1]).new_at(currIndex) = 0;
    }

    // On Left: u <- 0; v -> spread
    for (Range_Iterator<d> iterator(uLeft); iterator.Valid(); iterator.Next())
    {
        currIndex = T_INDEX() + iterator.Index();
        (*_v[1]).new_at(currIndex) = (*_v[1]).new_at(grid->Next_Cell(0, currIndex));
        _d->new_at(currIndex) = _d->new_at(grid->Next_Cell(0, currIndex));
        (*_v[0]).new_at(currIndex) = 0;
    }

    for (Range_Iterator<d> iterator(uRight); iterator.Valid(); iterator.Next())
    {
        currIndex = T_INDEX() + iterator.Index();
        (*_v[1]).new_at(currIndex) = (*_v[1]).new_at(grid->Previous_Cell(0, currIndex));
        _d->new_at(currIndex) = _d->new_at(grid->Previous_Cell(0, currIndex));
        (*_v[0]).new_at(currIndex) = 0;
    }
    for (Range_Iterator<d> iterator(vBottom); iterator.Valid(); iterator.Next())
    {
        currIndex = T_INDEX() + iterator.Index();
        (*_v[0]).new_at(currIndex) = (*_v[0]).new_at(grid->Next_Cell(1, currIndex));
        _d->new_at(currIndex) = _d->new_at(grid->Next_Cell(1, currIndex));
        (*_v[1]).new_at(currIndex) = 0;
    }
}

template <typename T, int d>
void FluidSolver<T, d>::spread2Ghost()
{
    T_INDEX max_corner = grid->Cell_Indices(number_of_ghost_cells).max_corner;
    int m = max_corner(0);
    int n = max_corner(1);

    // 0, 33 -> 33, 33
    Range<int, d> vTop(T_INDEX{0, n}, max_corner);
    // 0,0 -> 33, 0
    Range<int, d> uLeft(T_INDEX(0), T_INDEX{0, n});
    // 33,0 -> 33, 33
    Range<int, d> uRight(T_INDEX{m, 0}, max_corner);
    // 0,0 -> m,n
    Range<int, d> vBottom(T_INDEX(0), T_INDEX{m, 0});

    T_INDEX currIndex;
    // On Top: u -> ghost; v <- 0
    for (Range_Iterator<d> iterator(vTop); iterator.Valid(); iterator.Next())
    {
        currIndex = T_INDEX() + iterator.Index();
        (*_v[0]).new_at(currIndex) = (*_v[0]).new_at(grid->Previous_Cell(1, currIndex));
        _d->new_at(currIndex) = _d->new_at(grid->Previous_Cell(1, currIndex));
    }

    // On Left: u <- 0; v -> spread
    for (Range_Iterator<d> iterator(uLeft); iterator.Valid(); iterator.Next())
    {
        currIndex = T_INDEX() + iterator.Index();
        (*_v[1]).new_at(currIndex) = (*_v[1]).new_at(grid->Next_Cell(0, currIndex));
        _d->new_at(currIndex) = _d->new_at(grid->Next_Cell(0, currIndex));
    }

    for (Range_Iterator<d> iterator(uRight); iterator.Valid(); iterator.Next())
    {
        currIndex = T_INDEX() + iterator.Index();
        (*_v[1]).new_at(currIndex) = (*_v[1]).new_at(grid->Previous_Cell(0, currIndex));
        _d->new_at(currIndex) = _d->new_at(grid->Previous_Cell(0, currIndex));
    }

    for (Range_Iterator<d> iterator(vBottom); iterator.Valid(); iterator.Next())
    {
        currIndex = T_INDEX() + iterator.Index();
        (*_v[0]).new_at(currIndex) = (*_v[0]).new_at(grid->Next_Cell(1, currIndex));
        _d->new_at(currIndex) = _d->new_at(grid->Next_Cell(1, currIndex));
    }
}

template <typename T, int d>
void FluidSolver<T, d>::update(T timestep)
{
    //  std::cout<< "Running" <<std::endl;
    spread2Ghost();
    advection(timestep);
    // _d->printPhi();
    // _v[0]->printPhi();
    // _v[1]->printPhi();
    setBoundaryCondition();

    // Project(1000);
    // updateVelocity(timestep);

    flip();

    for (int i = 0; i < size_interior_domain; i++)
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