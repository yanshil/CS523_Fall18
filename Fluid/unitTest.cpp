#include <algorithm>
#include <math.h>
#include <iostream>

using namespace std;

double at(int p1, int p2)
{
    return p1 + p2;
}

/* Linear intERPolate between a and b for x ranging from 0 to 1 */
double lerp(double a, double b, double x)
{
    return a * (1.0 - x) + b * x;
}

/* Linear intERPolate on grid at coordinates (x, y).
     * Coordinates will be clamped to lie in simulation domain
     */
double lerp(double x, double y)
{
    double _ox = 0.5, _oy = 0.5; // center
    int _w = 128, _h = 128;

    x = min(max(x - _ox, 0.0), _w - 1.001);
    y = min(max(y - _oy, 0.0), _h - 1.001);
    int ix = (int)x;
    int iy = (int)y;
    x -= ix;
    y -= iy;

    cout<<"ox = "<<_ox<<",\toy = "<< _oy <<endl;

    cout<<"x = "<<x<<",\ty = "<< y<<endl;
    cout<<"ix = "<<ix<<",\tiy = "<< iy<<endl;
    
    double x00 = at(ix + 0, iy + 0), x10 = at(ix + 1, iy + 0);
    double x01 = at(ix + 0, iy + 1), x11 = at(ix + 1, iy + 1);

    return lerp(lerp(x00, x10, x), lerp(x01, x11, x), y);
}

main(int argc, char const *argv[])
{
    cout<<lerp(4, 4.5)<<endl;
    return 0;
}
