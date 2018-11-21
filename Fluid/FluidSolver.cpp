//!#####################################################################
//! \file FluidSolver.cpp
//!#####################################################################
#include "FluidSolver.h"
using namespace Nova;
//######################################################################
// TODO
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
    (*velocityField[1]).fill(0);
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
        divG[idx] = 0;

        for (int axis = 0; axis < d; axis++)
        {
            T_INDEX n_index = Next_Cell(axis, currIndex);
            if ((*grid).Inside_Domain(n_index))
            {
                divG[idx] += ((*velocityField[axis]).at(n_index) - (*velocityField[axis]).at(currIndex)) * (*grid).one_over_dX[axis];
                // divG[idx] += ((*velocityField[axis]).new_at(n_index) - (*velocityField[axis]).new_at(currIndex)) * (*grid).one_over_dX[axis];
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
                    AijSum += scale * pressure_solution[index2offset(p_index)];
                }
                if ((*grid).Inside_Domain(n_index))
                {
                    Aii -= scale;
                    AijSum += scale * pressure_solution[index2offset(n_index)];
                }
            }

            T newP = (divG[os] - AijSum) / Aii;

            maxDelta = std::max(maxDelta, fabs(pressure_solution[os] - newP));

            pressure_solution[os] = newP;
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
        // pressure_solution[index2offset(index)] = size * divG[index2offset(index)];

        for (int axis = 0; axis < d; axis++)
        {
            if ((*grid).Inside_Domain(Previous_Cell(axis, index)))
            {
                (*velocityField[axis]).new_at(index) -= timestep * (pressure_solution[index2offset(index)] - pressure_solution[index2offset(Previous_Cell(axis, index))]) *
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
        pressure_solution[index2offset(currIndex)] = pressure_solution[index2offset(Previous_Cell(1, currIndex))];
    }

    for (Range_Iterator<d> iterator(uLeft); iterator.Valid(); iterator.Next())
    {
        currIndex = T_INDEX() + iterator.Index();
        pressure_solution[index2offset(currIndex)] = pressure_solution[index2offset(Next_Cell(0, currIndex))];
    }

    for (Range_Iterator<d> iterator(uRight); iterator.Valid(); iterator.Next())
    {
        currIndex = T_INDEX() + iterator.Index();
        pressure_solution[index2offset(currIndex)] = pressure_solution[index2offset(Previous_Cell(0, currIndex))];
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
        divG[index2offset(currIndex)] = divG[index2offset(Previous_Cell(1, currIndex))];
    }

    for (Range_Iterator<d> iterator(uLeft); iterator.Valid(); iterator.Next())
    {
        currIndex = T_INDEX() + iterator.Index();
        divG[index2offset(currIndex)] = divG[index2offset(Next_Cell(0, currIndex))];
    }

    for (Range_Iterator<d> iterator(uRight); iterator.Valid(); iterator.Next())
    {
        currIndex = T_INDEX() + iterator.Index();
        divG[index2offset(currIndex)] = divG[index2offset(Previous_Cell(0, currIndex))];
    }
}
template <typename T, int d>
void FluidSolver<T, d>::update(T timestep)
{
    calculateDivergence();
    advection(timestep);

    Project(1000);
    updateVelocity(timestep);
    SetVelocityBoundary_as0();

    flip();

    for (int i = 0; i < size; i++)
    {
        pressure_solution[i] = 0;
    }

}
//######################################################################
template class Nova::FluidSolver<float,2>;
template class Nova::FluidSolver<float,3>;
#ifdef COMPILE_WITH_DOUBLE_SUPPORT
template class Nova::FluidSolver<double,2>;
template class Nova::FluidSolver<double,3>;
#endif