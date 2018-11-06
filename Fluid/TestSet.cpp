#include <nova/Tools/Grids/Grid.h>

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
    // T_INDEX v = T_INDEX(1);
    // std::cout<<v<<std::endl;
    // T_INDEX counts = T_INDEX(5);
    // std::cout<<offset(v, counts)<<std::endl;

    Vector<int, 2> v1 = Vector<int, 2>{0, -1};  // For Face Indicator == 0
    Vector<int, 2> v2 = Vector<int, 2>{0, 1};   // For The update velocity
    Vector<int, 1> v3 = Vector<int, 1>{0};  // For the rest

    int faceIndicator = 0;
    int _ii_ = 1;
    int rest_index = 3 - _ii_ -faceIndicator;

    for (int i = 0; i < v1.Size(); i++)
    {

        for (int j = 0; j < v2.Size(); j++)
        {

            for (int k = 0; k < v3.Size(); k++)
            {
                std::cout << v1(i) << ", " << v2(j) << ", " << v3(k) << std::endl;
            }
        }
    }

    return 0;
}
