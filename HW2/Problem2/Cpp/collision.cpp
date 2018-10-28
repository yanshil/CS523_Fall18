#include <Eigen/Geometry>
#include <iostream>

using namespace Eigen;

double sdf_box(Vector3d input, Vector3d origin, Vector3d width)
{
    //   vec3 d = abs(p) - b;
    //   return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));

    Vector3d dTmp;

    for (int i = 0; i < dTmp.size; i++)
    {
        dTmp[i] = std::abs(input[i]-origin[i]) - width[i];
    }

    // If in the box, this will be the SDF.
    double dTmpMax = dTmp.maxCoeff;

    double dTmpLength = 0;
    for (int i = 0; i < dTmp.size; i++)
    {
        if (dTmp(i) < 0)
            dTmp(i) = 0;
        dTmpLength += std::pow(dTmp(i), 2.0);
    }

    dTmpLength = std::sqrt(dTmpLength);

    return std::min(dTmpMax, 0.0) + dTmpLength;
}

