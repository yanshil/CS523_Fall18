#include "FluidQuantity.h"

FluidQuantity::FluidQuantity(Grid<T, d> grid, int faceIndicator)
{

    Phi = new float[grid.counts.Product];
    Phi_new = new float[grid.counts.Product];
    this->grid = &grid;

    for (int i = 0; i < d; i++)
    {
        this->faceIndicator = faceIndicator;
    }
}

FluidQuantity::~FluidQuantity()
{
    delete[] Phi;
    delete[] Phi_new;
}

TV FluidQuantity::computeVelocity(T_INDEX &index, FluidQuantity *FluidVelocity[d])
{
    TV velocity;
    // *Phi + offset(index, (*grid).counts);
    int ix = index[0], iy = index[1], iz = 0;
    if (d == 3)
        iz = index[2];

    // For Scalar
    if (faceIndicator == -1)
    {
        velocity[0] = 0.5 * (FluidVelocity[0]->at_(ix, iy, iz) + FluidVelocity[0]->at_(ix + 1, iy, iz));
        velocity[1] = 0.5 * (FluidVelocity[1]->at_(ix, iy, iz) + FluidVelocity[1]->at_(ix, iy + 1, iz));
        if (d == 3)
            velocity[2] = 0.5 * (FluidVelocity[2]->at_(ix, iy, iz) + FluidVelocity[2]->at_(ix, iy, iz + 1));

        return velocity;
    }

    // For Face Velocity
    for (int i = 0; i < d; i++)
    {
        // Copy value for that exact face
        if (faceIndicator == i)
        {
            velocity[i] = FluidVelocity[i]->at_(ix, iy, iz);
        }
        else
        {
            if(d==3) int rest_index = 3 - faceIndicator - i;
            
            // TODO
            int ox, oy, oz;

            double re = 0;
            re += FluidVelocity[i]->at_(ix+ox, iy+oy, iz+oz);
            re *= 0.25;
            velocity[i] = re;

        }

        return velocity;
    }
}
