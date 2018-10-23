# define NUM_CELLS 100
# define EPSILON 0.0001

#include <nova/Tools/Grids/Grid.h>
#include <nova/Tools/Parsing/Parse_Args.h>
#include <nova/Tools/Utilities/File_Utilities.h>
#include <iostream>

# include <Eigen/Geometry>

# include <float.h>      /* DBL_MAX */
# include <math.h>       /* sqrt */

using namespace Nova;

enum {d=2};
typedef float T;
using T_INDEX       = Vector<int,d>;

class CELL
{
    public:
    Eigen::Vector2d bottemleft;
    Eigen::Vector2d bottemright;
    Eigen::Vector2d topleft;
    Eigen::Vector2d topright;

    double dx;
    double dy;

};

/* Given the implicit function, 
get value of phi and stored themgrids's nodes */

double getPhi_analytic(Eigen::Vector2d input)
{
    double phi;

    /* Circle */
    phi = std::sqrt(std::pow(input(0), 2.0) + std::pow(input(1), 2.0)) - 1;

    /* Ellipse */


    /* Rectangular */


    return phi;
}

/* Given a point in space, return its stored cell */

CELL getCell(Eigen::Vector2d point, Grid<T,d> grid)
{
    CELL cell;
    
    cell.dx = grid.dX(0);
    cell.dy = grid.dX(1);

    cell.bottemleft(0) = floor((point(0) - grid.domain.min_corner(0)) / cell.dx) * cell.dx;
    cell.bottemleft(1) = floor((point(1) - grid.domain.min_corner(1)) / cell.dy) * cell.dy;

    cell.bottemright(0) = floor((point(0) - grid.domain.min_corner(0)) / cell.dx) * cell.dx;
    cell.bottemright(1) = floor(1 + (point(1) - grid.domain.min_corner(1)) / cell.dy) * cell.dy;

    cell.topleft(0) = floor(1 + (point(0) - grid.domain.min_corner(0)) / cell.dx) * cell.dx;
    cell.topleft(1) = floor((point(1) - grid.domain.min_corner(1)) / cell.dy) * cell.dy;

    cell.topright(0) = floor(1 + (point(0) - grid.domain.min_corner(0)) / cell.dx) * cell.dx;
    cell.topright(1) = floor(1 + (point(1) - grid.domain.min_corner(1)) / cell.dy) * cell.dy;

    // std::cout<<"Bottom Left = \n"<<cell.bottemleft<<std::endl;
    // std::cout<<"Bottom Right = \n"<<cell.bottemright<<std::endl;
    // std::cout<<"Top Left = \n"<<cell.topleft<<std::endl;
    // std::cout<<"Top Right = \n"<<cell.topright<<std::endl;

    return cell;
}


/* Given a point in space, return its phi value */
double getPhi(Eigen::Vector2d point, Grid<T,d> grid)
{
    /* Main Logic */
    double phi_p, phi00, phi01, phi10, phi11, phi_x0, phi_x1;
    double a, b;
    CELL cell;

    cell = getCell(point, grid);

    /* residual */
    a = point(0) - cell.bottemleft(0);
    b = point(1) - cell.bottemleft(1);

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
    phi_x0 = a * phi10 + (cell.dx - a) * phi00;
    phi_x1 = a * phi11 + (cell.dx - a) * phi01;

    /* result */
    phi_p = b * phi_x1 + (cell.dx - b) * phi_x0;

    return phi_p;
}


/* For points on cell, get nabla phi analytically */
Eigen::Vector2d getNablaPhi_analytic(Eigen::Vector2d point)
{
    Eigen::Vector2d np;

    /* Circle */
    np(0) = point(0) / std::sqrt(std::pow(point(0), 2.0) + std::pow(point(1), 2.0));
    np(1) = point(1) / std::sqrt(std::pow(point(0), 2.0) + std::pow(point(1), 2.0));

    return np;

}


/* Get nabla phi */
Eigen::Vector2d getNablaPhi(Eigen::Vector2d point, Grid<T,d> grid)
{
    Eigen::Vector2d np;
    double a, b, phi00;

    CELL cell = getCell(point, grid);
    phi00 = getPhi_analytic(cell.bottemleft);

    /* residual */
    a = point(0) - cell.bottemleft(0);
    b = point(1) - cell.bottemleft(1);

    /* If the point is on the cell vertex */
    if ((std::abs(a) < EPSILON) && (std::abs(b) < EPSILON))
    {
        np = getNablaPhi_analytic(point);
        std::cout << "Getting NablaPhi analytically =" << np << std::endl;
        return np;
    }

    np(0) = (getPhi(point, grid) - phi00) / a; // TODO
    np(1) = (getPhi(point, grid) - phi00) / b;
    
    return np;
}

void movePointToBoundary(Eigen::Vector2d point, Grid<T,d> grid)
{
    double currentDistance = DBL_MAX;   // Max Double from float.h
    Eigen::Vector2d nabla_phi;

    // Stop until phi is close to 0
    while(std::abs(currentDistance) > EPSILON)
    {
        std::cout<<"Current point location: \n"<< point <<std::endl;
        
        currentDistance = getPhi(point, grid);
        std::cout<<"Current Phi: \n"<< currentDistance <<std::endl;

        nabla_phi = getNablaPhi(point, grid);
        std::cout << "nable_phi = \n" << nabla_phi << std::endl;


        point(0) = point(0) - currentDistance * nabla_phi(0); //TODo
        point(1) = point(0) - currentDistance * nabla_phi(1);
    }

    std::cout<<"Final point location: \n"<< point <<std::endl;

}

int main(int argc,char** argv)
{
    // enum {d=2};
    // typedef float T;
    // using T_INDEX       = Vector<int,d>;

    Parse_Args parse_args;
    if(d==2) parse_args.Add_Vector_2D_Argument("-size",Vector<double,2>(64),"","Grid resolution");
    else parse_args.Add_Vector_3D_Argument("-size",Vector<double,3>(64),"","Grid resolution");
    parse_args.Add_String_Argument("-o",".","","Output directory");
    parse_args.Parse(argc,argv);

    T_INDEX counts;
    if(d==2){auto counts_2d=parse_args.Get_Vector_2D_Value("-size");for(int v=0;v<d;++v) counts(v)=counts_2d(v);}
    else{auto counts_3d=parse_args.Get_Vector_3D_Value("-size");for(int v=0;v<d;++v) counts(v)=counts_3d(v);}
    std::string output_directory=parse_args.Get_String_Value("-o");

    Log::cout<<"Counts: "<<counts<<std::endl;

    using TV = Vector<T,d>;
    // Grid<T,d> grid(counts,Range<T,d>(TV(-10.0),TV(10.0)));
    Grid<T,d> grid(counts,Range<T,d>::Unit_Box());
    // File_Utilities::Write_To_File(output_directory+"/grid.grid",grid);

    Eigen::Vector2d point;
    point << 2, 2;

    double phi_p = getPhi(point, grid);
    if(phi_p == 0) std::cout << "Given point is on the boundary" << std::endl;
    if(phi_p < 0) std::cout << "Given point is in the boundary" << std::endl;
    if(phi_p > 0) std::cout << "Given point is out of boundary" << std::endl;

    movePointToBoundary(point, grid);
    // 0.071067


    return 0;
}

// g++ test.cpp -I ./include/ -std=gnu++11 -lboost_system
