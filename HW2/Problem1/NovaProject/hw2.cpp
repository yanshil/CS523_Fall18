#define EPSILON 0.0001
#define ELLIPSE_A 4
#define ELLIPSE_B 2
#define ELLIPSE_C 3

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

// typedef Vector<double, 2> Vector2d;

/* Global variabble for shape selection */
int shape = 0;

/*-------------------------------------------*/

class CELL
{
  public:
    Vector<double, d> c000;
    Vector<double, d> c010;
    Vector<double, d> c100;
    Vector<double, d> c111;

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

double sdf_circle(Vector<double, d> inputn)
{
    double tmp = 0;

    for (int i = 0; i < d; i++)
    {
        tmp += std::pow(input(i), 2.0)
    }

    return std::sqrt(tmp) - 1;
    // 2d: phi = std::sqrt(std::pow(input(0)-origin(0), 2.0) + std::pow(input(1)-origin(0), 2.0)) - 1;
}

double sdf_ellipse(Vector<double, d> input, Vector<double, d> ellipse_para)
{
    double tmp = 0;

    for (int i = 0; i < d; i++)
    {
        tmp += std::pow(input(i), 2.0) / ellipse_para(i);
    }
    return std::sqrt(tmp) - 1;
    // 2d: phi = std::sqrt(std::pow(input(0), 2.0) / ELLIPSE_A + std::pow(input(1), 2.0) / ELLIPSE_B) - 1;
}

/* input: point in space; origin: center of box; width: "radius" of the box */
double sdf_box(Vector<double, d> input, Vector<double, d> width)
{
    Vector<double, d> dTmp;
    dTmp = (input).Abs() - width;

    // If in the box, this will be the SDF.
    double dTmpMax = dTmp.Max();

    double dTmpLength = 0;
    for (int i = 0; i < d; i++)
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

double getPhi_analytic(Vector<double, d> input)
{
    double phi = 0;

    /* Circle */
    if (shape == 1)
        phi = sdf_circle(input);
    /* Ellipse */
    // Ref: http://impact.byu.edu/Image%20Processing%20Seminar/W02_SignedDistanceFunction.pdf
    if (shape == 2)
    {
        Vector<double, d> ellipse_para;

        ellipse_para(0) = ELLIPSE_A;
        ellipse_para(1) = ELLIPSE_B;
        if(d==3) ellipse_para(2) = ELLIPSE_C;

        phi = sdf_ellipse(input, ellipse_para);
    }
    // phi = std::sqrt(std::pow(input(0), 2.0) / ELLIPSE_A + std::pow(input(1), 2.0) / ELLIPSE_B) - 1;
    /* Retangular */
    // Ref: https://stackoverflow.com/questions/35326366/glsl-cube-signed-distance-field-implementation-explanation
    // Ref: https://www.alanzucconi.com/2016/07/01/signed-distance-functions/#part3
    if (shape == 3)
    {
        phi = sdf_box(input, Vector<double, d>({0.0, 0.0}), Vector<double, d>({1.0, 1.0}));
    }

    // std::cout<< "phi("  << input << ") = "<< phi <<std::endl;

    return phi;
}

/* Given a point in space, return its stored cell */

CELL getCell(Vector<double, d> point, Grid<T, d> grid)
{
    CELL cell;

    cell.c000(0) = floor((point(0) - grid.domain.min_corner(0)) / grid.dX(0)) * grid.dX(0);
    cell.c000(1) = floor((point(1) - grid.domain.min_corner(1)) / grid.dX(1)) * grid.dX(1);

    cell.c010(0) = floor(1 + (point(0) - grid.domain.min_corner(0)) / grid.dX(0)) * grid.dX(0);
    cell.c010(1) = floor((point(1) - grid.domain.min_corner(1)) / grid.dX(1)) * grid.dX(1);

    cell.c100(0) = floor((point(0) - grid.domain.min_corner(0)) / grid.dX(0)) * grid.dX(0);
    cell.c100(1) = floor(1 + (point(1) - grid.domain.min_corner(1)) / grid.dX(1)) * grid.dX(1);

    cell.c111(0) = floor(1 + (point(0) - grid.domain.min_corner(0)) / grid.dX(0)) * grid.dX(0);
    cell.c111(1) = floor(1 + (point(1) - grid.domain.min_corner(1)) / grid.dX(1)) * grid.dX(1);

    return cell;
}

/* Given a point in space, return its phi value */
double getPhi(Vector<double, d> point, Grid<T, d> grid)
{
    /* Main Logic */
    double phi_p, phi00, phi01, phi10, phi11, phi_x0, phi_x1;
    double a, b;
    CELL cell;

    cell = getCell(point, grid);

    /* residual */
    a = (point(0) - cell.c000(0)) / grid.dX(0);
    b = (point(1) - cell.c000(1)) / grid.dX(1);
    //TODO 3D

    /* If the point is on the cell vertex */
    if ((std::abs(a) < EPSILON) && (std::abs(b) < EPSILON))
    {
        phi_p = getPhi_analytic(point);
        std::cout << "Getting Phi analytically =" << phi_p << std::endl;
        return phi_p;
    }

    /* Get phi's analytic result for points on cell boundaries */
    phi00 = getPhi_analytic(cell.c000);
    phi01 = getPhi_analytic(cell.c010);
    phi10 = getPhi_analytic(cell.c100);
    phi11 = getPhi_analytic(cell.c111);

    //TODO: 3D

    /* Interpolate temparory value */
    phi_x0 = a * phi01 + (1 - a) * phi00;
    phi_x1 = a * phi11 + (1 - a) * phi10;

    /* result */
    phi_p = b * phi_x1 + (1 - b) * phi_x0;

    return phi_p;
}

/* Get nabla phi */
Vector<double, d> getNablaPhi(Vector<double, d> point, Grid<T, d> grid)
{
    Vector<double, d> np;
    double phi00, phi01, phi10;

    CELL cell = getCell(point, grid);
    phi00 = getPhi_analytic(cell.c000);
    phi01 = getPhi_analytic(cell.c010);
    phi10 = getPhi_analytic(cell.c100);

    np(0) = (phi01 - phi00) / grid.dX(0);
    np(1) = (phi10 - phi00) / grid.dX(1);
    // TODO: 3D

    // std::cout << "phi00(" << cell.c000 << ") = " << phi00 << std::endl;
    // std::cout << "phi01(" << cell.c010 << ") = " << phi01 << std::endl;
    // std::cout << "phi10(" << cell.c100 << ") = " << phi10 << std::endl;

    // std::cout << "(phi01 - phi00)" << (phi01 - phi00) << std::endl;
    // std::cout << "(phi10 - phi00)" << (phi10 - phi00) << std::endl;
    // std::cout << "np = " << np << std::endl;

    return np;
}

void movePointToBoundary(Vector<double, d> point, Grid<T, d> grid)
{
    double currentDistance = DBL_MAX; // Max Double from float.h
    Vector<double, d> nabla_phi;
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
        //TODO: 3D
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

    Vector<double, d> point;
    for (int i = 0; i < d; i++)
        std::cin >> point(i);

    movePointToBoundary(point, grid);
    //0.7071067

    return 0;
}
