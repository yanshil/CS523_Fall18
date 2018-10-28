#define EPSILON 0.005
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
    d = 3
};
typedef double T;
using T_INDEX = Vector<int, d>;

/* Global variabble for shape selection */
int shape = 0;

/*-------------------------------------------*/

class CELL
{

    /* 
    Verteices for 3D Grid Cell
    Tagged with c:xyz

    For 2D Grid Cell, only use c: xy0
    */
  public:
    Vector<double, d> c000;
    Vector<double, d> c001; // z
    Vector<double, d> c010; // y
    Vector<double, d> c011;
    Vector<double, d> c100; // x
    Vector<double, d> c101;
    Vector<double, d> c110;
    Vector<double, d> c111;

    CELL()
    {
        c000 = Vector<double, d>(0.0);
        c001 = Vector<double, d>(0.0);
        c010 = Vector<double, d>(0.0);
        c011 = Vector<double, d>(0.0);
        c100 = Vector<double, d>(0.0);
        c101 = Vector<double, d>(0.0);
        c110 = Vector<double, d>(0.0);
        c111 = Vector<double, d>(0.0);
    }
    // using point = Vector<double, d>;

    // Vector<double, d> cellVertex[8];
};

void menu_GetShape()
{
    while (shape == 0)
    {
        if (d == 2)
        {
            std::cout << "Please select a shape: [default with origin (0,0)]" << std::endl;
            std::cout << "1. Circle (r = 1)" << std::endl;
            std::cout << "2. Ellipse (A = 4, B = 2)" << std::endl;
            std::cout << "3. Square (width = 1)" << std::endl;
        }
        else
        {
            std::cout << "Please select a shape: [default with origin (0,0,0)]" << std::endl;
            std::cout << "1. Shpere (r = 1)" << std::endl;
            std::cout << "2. Ellipsoid (A = 4, B = 2, C = 3)" << std::endl;
            std::cout << "3. Cubic (width = 1)" << std::endl;
        }

        std::cin >> shape;

        switch (shape)
        {
        case 1:
            if (d == 2)
                std::cout << "Circle is selected" << std::endl;
            else
                std::cout << "Shpere is selected" << std::endl;
            break;
        case 2:
            if (d == 2)
                std::cout << "Ellipse is selected" << std::endl;
            else
                std::cout << "Ellipsoid is selected" << std::endl;
            break;
        case 3:
            if (d == 2)
                std::cout << "Square is selected" << std::endl;
            else
                std::cout << "Cubic is selected" << std::endl;
            break;
        default:
            shape = 0;
            std::cout << "Unvalid selection. Please enter '1', '2' or '3'.---------------\n"
                      << std::endl;
            break;
        }
    }
}

double sdf_circle(Vector<double, d> input)
{
    double tmp = 0;

    for (int i = 0; i < d; i++)
    {
        tmp += std::pow(input(i), 2.0);
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
        if (d == 3)
            ellipse_para(2) = ELLIPSE_C;

        phi = sdf_ellipse(input, ellipse_para);
    }
    // phi = std::sqrt(std::pow(input(0), 2.0) / ELLIPSE_A + std::pow(input(1), 2.0) / ELLIPSE_B) - 1;
    /* Retangular */
    // Ref: https://stackoverflow.com/questions/35326366/glsl-cube-signed-distance-field-implementation-explanation
    // Ref: https://www.alanzucconi.com/2016/07/01/signed-distance-functions/#part3
    if (shape == 3)
    {
        phi = sdf_box(input, Vector<double, d>(1.0));
    }

    return phi;
}

/* Given a point in space, return its stored cell */

CELL getCell(Vector<double, d> point, Grid<T, d> grid)
{
    CELL cell;

    for (int i = 0; i < d; i++)
    {
        cell.c000(i) = floor((point(i) - grid.domain.min_corner(i)) * grid.one_over_dX(i)) * grid.dX(i);

        cell.c001(i) = cell.c000(i);
        cell.c010(i) = cell.c000(i);
        cell.c011(i) = cell.c000(i);
        cell.c100(i) = cell.c000(i);
        cell.c101(i) = cell.c000(i);
        cell.c110(i) = cell.c000(i);
        cell.c111(i) = cell.c000(i);
    }

    cell.c100(0) += grid.dX(0);
    cell.c010(1) += grid.dX(1);

    cell.c110(0) += grid.dX(0);
    cell.c110(1) += grid.dX(1);

    if (d == 3)
    {
        cell.c001(2) += grid.dX(2);

        cell.c011(1) += grid.dX(1);
        cell.c011(2) += grid.dX(2);

        cell.c101(0) += grid.dX(0);
        cell.c101(2) += grid.dX(2);

        cell.c111(0) += grid.dX(0);
        cell.c111(1) += grid.dX(1);
        cell.c111(2) += grid.dX(2);
    }

    // std::cout << cell.c000 << std::endl;
    // std::cout << cell.c100 << std::endl;
    // std::cout << cell.c010 << std::endl;
    // std::cout << cell.c001 << std::endl;

    return cell;
}

/* Given a point in space, return its phi value */
double getPhi(Vector<double, d> point, Grid<T, d> grid)
{
    /* Main Logic */
    double phi_p;
    double a, b, c;
    CELL cell;

    cell = getCell(point, grid);

    /* residual */
    a = (point(0) - cell.c000(0)) / grid.dX(0);
    b = (point(1) - cell.c000(1)) / grid.dX(1);
    if (d == 3)
        c = (point(2) - cell.c000(2)) / grid.dX(2);

    /* If the point is on the cell vertex */
    if ((std::abs(a) < EPSILON) && (std::abs(b) < EPSILON))
    {
        phi_p = getPhi_analytic(point);
        // std::cout << "(On Cell Boundary!) Getting Phi analytically =" << phi_p << std::endl;
        return phi_p;
    }

    /* Get phi's analytic result for points on cell boundaries */
    double phi000, phi001, phi010, phi011;
    phi000 = getPhi_analytic(cell.c000);
    phi001 = getPhi_analytic(cell.c001);
    phi010 = getPhi_analytic(cell.c010);
    phi011 = getPhi_analytic(cell.c011);

    if (d == 2)
    {
        double phi_x0, phi_x1;
        phi_x0 = a * phi001 + (1 - a) * phi000;
        phi_x1 = a * phi011 + (1 - a) * phi010;

        return b * phi_x1 + (1 - b) * phi_x0;
    }
    else
    {
        double phi100, phi101, phi110, phi111;
        phi100 = getPhi_analytic(cell.c100);
        phi101 = getPhi_analytic(cell.c101);
        phi110 = getPhi_analytic(cell.c110);
        phi111 = getPhi_analytic(cell.c111);

        double phi_x00, phi_x10, phi_x01, phi_x11, phi_0y0, phi_0y1;
        // double phi_000, phi_001, phi_010, phi_011, phi_100, phi_101, phi_110, phi_111;
        phi_x00 = a * phi100 + (1 - a) * phi000;
        phi_x10 = a * phi110 + (1 - a) * phi010;
        phi_x01 = a * phi101 + (1 - a) * phi001;
        phi_x11 = a * phi111 + (1 - a) * phi011;

        phi_0y0 = b * phi_x10 + (1 - b) * phi_x00;
        phi_0y1 = b * phi_x11 + (1 - b) * phi_x01;

        /* result */
        return c * phi_0y0 + (1 - c) * phi_0y1;
    }
}

/* Get nabla phi */
Vector<double, d> getNablaPhi(Vector<double, d> point, Grid<T, d> grid)
{
    Vector<double, d> np;
    double phi000, phi100, phi010, phi001;

    CELL cell = getCell(point, grid);
    // phi00 = getPhi_analytic(cell.c000);
    // phi01 = getPhi_analytic(cell.c010);
    // phi10 = getPhi_analytic(cell.c100);

    phi000 = getPhi_analytic(cell.c000);
    if (d == 2)
    {
        phi100 = getPhi_analytic(cell.c100);
        phi010 = getPhi_analytic(cell.c010);
    }
    else
    {
        phi100 = getPhi_analytic(cell.c100);
        phi010 = getPhi_analytic(cell.c010);
        phi001 = getPhi_analytic(cell.c001);

        np(2) = (phi001 - phi000) * grid.one_over_dX(2);
    }

    np(0) = (phi100 - phi000) * grid.one_over_dX(0);
    np(1) = (phi010 - phi000) * grid.one_over_dX(1);

    return np;
}

void movePointToBoundary(Vector<double, d> point, Grid<T, d> grid)
{
    Vector<double, d> originalLocation = point;
    double currentDistance = DBL_MAX; // Max Double from float.h
    Vector<double, d> nabla_phi;
    int iteration = 0;

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
        std::cout << "==> Iterate to (" << point << ")" << std::endl;

        // std::cout << "phi= " << currentDistance << std::endl;
        // std::cout << "Estimated nable_phi = " << nabla_phi << std::endl;

        for (int i = 0; i < d; i++)
        {
            point(i) = point(i) - currentDistance * nabla_phi(i);
        }

        iteration++;

        if (iteration >= 200)
        {
            std::cout << "Too many iterations!" << std::endl;
            double tmp;
            if (d == 2)
            {
                if (shape == 1)
                    tmp = std::pow(point(0), 2.0) + std::pow(point(1), 2.0);
                if (shape == 2)
                    tmp = std::pow(point(0), 2.0) / ELLIPSE_A + std::pow(point(1), 2.0) / ELLIPSE_B;
                std::cout << "Radius(final position)" << tmp << std::endl;
            }
            else
            {
                if (shape == 1)
                    tmp = std::pow(point(0), 2.0) + std::pow(point(1), 2.0) + std::pow(point(2), 2.0);
                if (shape == 2)
                    tmp = std::pow(point(0), 2.0) / ELLIPSE_A + std::pow(point(1), 2.0) / ELLIPSE_B + std::pow(point(2), 2.0) / ELLIPSE_C;

                std::cout << "Radius(final position)" << tmp << std::endl;
            }

            break;
        }
    }

    std::cout << "Final point location: \n"
              << point << std::endl;
}

int main(int argc, char **argv)
{
    std::cout << "Current Dimension = " << d << std::endl;
    std::cout << "\tNote: 2D / 3D is modified from #20 of code. (Setting enum{d=2} or enum{d=3})." << std::endl;

    /* Customize Grid Size with Parsing Arguments */
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
    std::cout << "Cell width: " << grid.dX(0) << "\n---------------" << std::endl;

    menu_GetShape();

    /* Input arbitrary point */
    std::cout << "Please enter a coordinate to evaluate: " << std::endl;

    Vector<double, d> point;
    for (int i = 0; i < d; i++)
        std::cin >> point(i);

    std::cout << "Input Point: (" << point << ")" << std::endl;

    /* Main Logic for moving points */
    movePointToBoundary(point, grid);

    return 0;
}
