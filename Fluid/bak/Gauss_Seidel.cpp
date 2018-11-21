#include "Fluid.h"

using namespace Nova;

main(int argc, char const *argv[])
{

    return 0;
}

void FluidSolver::projection(int limit, double timestep)
{
    double maxDelta;
    int limit;
    double c = timestep * (*grid).one_over_dX.Product();

    for (int iteration = 0; iteration < limit; iteration++)
    {
        T_INDEX currIndex;
        for (Range_Iterator<d> iterator(Range<int, d>(T_INDEX(1), (*grid).counts)); iterator.Valid(); iterator.Next())
        {
            currIndex = T_INDEX() + iterator.Index();
            int os = index2offset(currIndex);

            double Aijsum = 0, Aii = 0;

            for (int axis = 0; axis < d; axis++)
            {
                if (currIndex[axis] > 0)
                {
                    Aii += c;
                    Aijsum -= c * pressure_solution[index2offset(Previous_Cell(axis, currIndex))];
                }
                if (currIndex[axis] < (*grid).counts[axis])
                {
                    Aii += c;
                    Aijsum -= c * pressure_solution[index2offset(Next_Cell(axis, currIndex))];
                }

                double newP = (divG[os] - Aijsum) / Aii;

                maxDelta = std::max(maxDelta, fabs(pressure_solution[os] - newP));

                pressure_solution[os] = newP;
            }
        }

        if (maxDelta < 1e-05)
        {
            printf("Exiting solver after %d iterations, maximum change is %f\n", iteration, maxDelta);
            return;
        }

    }
    printf("Exceeded budget of %d iterations, maximum change was %f\n", limit, maxDelta);
}
