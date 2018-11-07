# include <nova/Tools/Grids/Grid.h>
# include "FluidQuantity.h"

using namespace Nova;

enum
{
    d = 3
};
typedef float T;
typedef Vector<T, d> TV;
typedef Vector<int, d> T_INDEX;

/* Auxiliary Function: Get exact offset for Column-based 1D array */
int offset(const T_INDEX &index, const T_INDEX &counts)
{

    int os = index[1] * counts[0] + index[0];
    if (d == 3)
        os += index[2] * counts[0] * counts[1];
    return os;

    // return (z * xSize * ySize) + (y * xSize) + x;
}

main(int argc, char const *argv[])
{

    
    if (true) {
        FluidQuantity u;
    }
    

    if (false)
    {
        T_INDEX v = T_INDEX(1);
        std::cout << v << std::endl;
        T_INDEX counts = T_INDEX(5);
        std::cout << offset(v, counts) << std::endl;
    }

    if (false)
    {
        Vector<int, 2> v_FI = Vector<int, 2>{0, -1}; // For Face Indicator == 0
        Vector<int, 2> v_UV = Vector<int, 2>{0, 1};  // For The update velocity
        Vector<int, 1> v3 = Vector<int, 1>{0};       // For the rest

        int faceIndicator = 0;
        int uv = 1;
        int rest_index;
        if (d == 3)
            rest_index = 3 - uv - faceIndicator;

        TV offset = TV(0);

        for (int i = 0; i < 2; i++)
        {

            for (int j = 0; j < 2; j++)
            {
                offset[faceIndicator] = v_FI[i];
                offset[uv] = v_UV[j];
                if (d == 3)
                {
                    offset[rest_index] = 0;
                }

                std::cout << offset << std::endl;
            }
        }
    }

    return 0;
}
