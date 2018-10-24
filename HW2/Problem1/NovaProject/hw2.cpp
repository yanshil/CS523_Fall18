#define EPSILON 0.0001
#define ELLIPSE_A 4
#define ELLIPSE_B 2

#include <nova/Tools/Grids/Grid.h>
#include <nova/Tools/Parsing/Parse_Args.h>
#include <nova/Tools/Utilities/File_Utilities.h>
#include <iostream>

#include <float.h> /* DBL_MAX */
#include <math.h>  /* sqrt */

using namespace Nova;

/*-------- Global typedef & varaible------------*/

enum
{
    d = 2
};
typedef double T;
using T_INDEX = Vector<int, d>;

typedef Vector<double, 2> Vector2d;

/* Global variabble for shape selection */
int shape = 0;

/*-------------------------------------------*/

class CELL
{
  public:
    Vector2d bottemleft;
    Vector2d bottemright;
    Vector2d topleft;
    Vector2d topright;

    double dx;
    double dy;
};

void menu_GetShape()
{
    while (shape == 0)
    {
        std::cout << "Please select a shape:" << std::endl;
        std::cout << "1. Circle [default with origin (0,0) ]" << std::endl;
        std::cout << "2. Ellipse [default with origin (0,0) ]" << std::endl;
        std::cout << "3. Rectangular [default with origin (0,0) ]" << std::endl;

        std::cin >> shape;

        switch (shape)
        {
        case 1:
            std::cout << "Circle is selected" << std::endl;
            break;
        case 2:
            std::cout << "Ellipse is selected" << std::endl;
            break;
        case 3:
            std::cout << "Rectangular is selected" << std::endl;
            break;
        default:
            shape = 0;
            std::cout << "Unvalid selection. Please enter '1', '2' or '3'.---------------\n"
                      << std::endl;
            break;
        }
    }
}

/* input: point in space; origin: center of box; width: "radius" of the box */
double sdf_box(Vector2d input, Vector2d origin, Vector2d width)
{
    //   vec3 d = abs(p) - b;
    //   return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));

    Vector2d dTmp;
    dTmp = (input - origin).Abs() - width;

    // If in the box, this will be the SDF.
    double dTmpMax = dTmp.Max();

    double dTmpLength = 0;
    for (int i = 0; i < dTmp.Size(); i++)
    {
        if (dTmp(i) < 0)
            dTmp(i) = 0;
        dTmpLength += std::pow(dTmp(i), 2.0);
    }

    dTmpLength = std::sqrt(dTmpLength);

    return std::min(dTmpMax, 0.0) + dTmpLength;
}

/* Given the implicit function, 
get value of phi and stored themgrids's nodes */

double getPhi_analytic(Vector2d input)
{
    double phi = 0;

    /* Circle */
    if (shape == 1)
        phi = std::sqrt(std::pow(input(0), 2.0) + std::pow(input(1), 2.0)) - 1;
    /* Ellipse */
    // Ref: http://impact.byu.edu/Image%20Processing%20Seminar/W02_SignedDistanceFunction.pdf
    if (shape == 2)
        phi = std::sqrt(std::pow(input(0), 2.0) / ELLIPSE_A + std::pow(input(1), 2.0) / ELLIPSE_B) - 1;
    /* Retangular */
    // Ref: https://stackoverflow.com/questions/35326366/glsl-cube-signed-distance-field-implementation-explanation
    // Ref: https://www.alanzucconi.com/2016/07/01/signed-distance-functions/#part3
    if (shape == 3)
    {
        phi = sdf_box(input, Vector2d({0.0, 0.0}), Vector2d({1.0, 1.0}));
    }

    // std::cout<< "phi("  << input << ") = "<< phi <<std::endl;

    return phi;
}

/* Given a point in space, return its stored cell */

CELL getCell(Vector2d point, Grid<T, d> grid)
{
    CELL cell;

    cell.dx = grid.dX(0);
    cell.dy = grid.dX(1);

    cell.bottemleft(0) = floor((point(0) - grid.domain.min_corner(0)) / cell.dx) * cell.dx;
    cell.bottemleft(1) = floor((point(1) - grid.domain.min_corner(1)) / cell.dy) * cell.dy;

    cell.bottemright(0) = floor(1 + (point(0) - grid.domain.min_corner(0)) / cell.dx) * cell.dx;
    cell.bottemright(1) = floor((point(1) - grid.domain.min_corner(1)) / cell.dy) * cell.dy;

    cell.topleft(0) = floor((point(0) - grid.domain.min_corner(0)) / cell.dx) * cell.dx;
    cell.topleft(1) = floor(1 + (point(1) - grid.domain.min_corner(1)) / cell.dy) * cell.dy;

    cell.topright(0) = floor(1 + (point(0) - grid.domain.min_corner(0)) / cell.dx) * cell.dx;
    cell.topright(1) = floor(1 + (point(1) - grid.domain.min_corner(1)) / cell.dy) * cell.dy;

    // std::cout<<"cell.bottomleft = ("<<cell.bottemleft<<")"<<std::endl;
    // std::cout<<"cell.bottemright = ("<<cell.bottemright<<")"<<std::endl;
    // std::cout<<"cell.topleft = ("<<cell.topleft<<")"<<std::endl;
    // std::cout<<"cell.topright = ("<<cell.topright<<")"<<std::endl;

    return cell;
}

/* Given a point in space, return its phi value */
double getPhi(Vector2d point, Grid<T, d> grid)
{
    /* Main Logic */
    double phi_p, phi00, phi01, phi10, phi11, phi_x0, phi_x1;
    double a, b;
    CELL cell;

    cell = getCell(point, grid);

    /* residual */
    a = (point(0) - cell.bottemleft(0)) / cell.dx;
    b = (point(1) - cell.bottemleft(1)) / cell.dy;

    /* If the point is on the cell vertex */
    if ((std::abs(a) < EPSILON) && (std::abs(b) < EPSILON))
    {
        phi_p = getPhi_analytic(point);
        std::cout << "Getting Phi analytically =" << phi_p << std::endl;
        return phi_p;
    }

    /* Get phi's analytic result for points on cell boundaries */
    phi00 = getPhi_analytic(cell.bottemleft);
    phi01 = getPhi_analytic(cell.bottemright);
    phi10 = getPhi_analytic(cell.topleft);
    phi11 = getPhi_analytic(cell.topright);

    /* Interpolate temparory value */
    phi_x0 = a * phi01 + (1 - a) * phi00;
    phi_x1 = a * phi11 + (1 - a) * phi10;

    /* result */
    phi_p = b * phi_x1 + (1 - b) * phi_x0;

    return phi_p;
}

/* Get nabla phi */
Vector2d getNablaPhi(Vector2d point, Grid<T, d> grid)
{
    Vector2d np;
    double phi00, phi01, phi10;

    CELL cell = getCell(point, grid);
    phi00 = getPhi_analytic(cell.bottemleft);
    phi01 = getPhi_analytic(cell.bottemright);
    phi10 = getPhi_analytic(cell.topleft);

    np(0) = (phi01 - phi00) / cell.dx;
    np(1) = (phi10 - phi00) / cell.dy;

    // std::cout << "phi00(" << cell.bottemleft << ") = " << phi00 << std::endl;
    // std::cout << "phi01(" << cell.bottemright << ") = " << phi01 << std::endl;
    // std::cout << "phi10(" << cell.topleft << ") = " << phi10 << std::endl;

    // std::cout << "(phi01 - phi00)" << (phi01 - phi00) << std::endl;
    // std::cout << "(phi10 - phi00)" << (phi10 - phi00) << std::endl;
    // std::cout << "np = " << np << std::endl;

    return np;
}

void movePointToBoundary(Vector2d point, Grid<T, d> grid)
{
    double currentDistance = DBL_MAX; // Max Double from float.h
    Vector2d nabla_phi;
    int iteration = 1;

    /* Display Out/In/On information before iteration */
    double phi_p = getPhi(point, grid);

    std::cout << "---------------" << std::endl;

    if (phi_p == 0)
        std::cout << "* Given point is on the boundary" << std::endl;
    if (phi_p < 0)
        std::cout << "* Given point is in the boundary" << std::endl;
    if (phi_p > 0)
        std::cout << "* Given point is out of boundary" << std::endl;

    // Stop until phi is close to 0
    while (std::abs(currentDistance) > EPSILON)
    {
        std::cout << "---------------\n"
                  << "Iteration " << iteration << std::endl;
        // std::cout << "Iteration " << iteration << std::endl;

        currentDistance = getPhi(point, grid);
        nabla_phi = getNablaPhi(point, grid);
        std::cout << "* Current point location: (" << point << ")" << std::endl;
        // std::cout << "phi= " << currentDistance << std::endl;
        // std::cout << "Estimated nable_phi = " << nabla_phi << std::endl;

        point(0) = point(0) - currentDistance * nabla_phi(0);
        point(1) = point(1) - currentDistance * nabla_phi(1);
        iteration++;

        if (iteration >= 50)
        {
            std::cout << "Too many iterations!" << std::endl;
            break;
        }
    }

    std::cout << "---------------\nFinal point location: \n"
              << point << std::endl;
}

int main(int argc, char **argv)
{

    Parse_Args parse_args;
    if (d == 2)
        parse_args.Add_Vector_2D_Argument("-size", Vector<double, 2>(500), "", "Grid resolution");
    else
        parse_args.Add_Vector_3D_Argument("-size", Vector<double, 3>(500), "", "Grid resolution");
    parse_args.Add_String_Argument("-o", ".", "", "Output directory");
    parse_args.Parse(argc, argv);

    T_INDEX counts;
    if (d == 2)
    {
        auto counts_2d = parse_args.Get_Vector_2D_Value("-size");
        for (int v = 0; v < d; ++v)
            counts(v) = counts_2d(v);
    }
    else
    {
        auto counts_3d = parse_args.Get_Vector_3D_Value("-size");
        for (int v = 0; v < d; ++v)
            counts(v) = counts_3d(v);
    }
    std::string output_directory = parse_args.Get_String_Value("-o");

    // Log::cout << "Counts: " << counts << std::endl;

    Grid<T, d> grid(counts, Range<T, d>::Unit_Box());
    // File_Utilities::Write_To_File(output_directory+"/grid.grid",grid);
    std::cout << "Cell size: " << grid.dX(0) << std::endl;

    menu_GetShape();

    std::cout << "Please enter a coordinate to evaluate: " << std::endl;

    Vector2d point;
    for (int i = 0; i < d; i++)
        std::cin >> point(i);

    movePointToBoundary(point, grid);
    //0.7071067

    return 0;
}
