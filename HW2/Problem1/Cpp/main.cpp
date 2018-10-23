# define NUM_CELLS 100
# define EPSILON 0.0001

# include <iostream>
# include <Eigen/Geometry>

# include <float.h>      /* DBL_MAX */
# include <math.h>       /* sqrt */


/*   Uniform Grid    */

class GRID
{
    public:

    Eigen::Vector2d origin;
    Eigen::Vector2d count;
    double dx;

    GRID()
    {
        origin << 0,0;
        count << NUM_CELLS, NUM_CELLS;
        dx = 0.001;
    }
};

GRID grid;

class CELL
{
    public:
    Eigen::Vector2d bottemleft;
    Eigen::Vector2d bottemright;
    Eigen::Vector2d topleft;
    Eigen::Vector2d topright;

    double dx;

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

CELL getCell(Eigen::Vector2d point)
{
    CELL cell;
    
    cell.dx = grid.dx;

    cell.bottemleft(0) = floor((point(0) - grid.origin(0)) / cell.dx) * cell.dx;
    cell.bottemleft(1) = floor((point(1) - grid.origin(1)) / cell.dx) * cell.dx;

    cell.bottemright(0) = floor((point(0) - grid.origin(0)) / cell.dx) * cell.dx;
    cell.bottemright(1) = floor(1 + (point(1) - grid.origin(1)) / cell.dx) * cell.dx;

    cell.topleft(0) = floor(1 + (point(0) - grid.origin(0)) / cell.dx) * cell.dx;
    cell.topleft(1) = floor((point(1) - grid.origin(1)) / cell.dx) * cell.dx;

    cell.topright(0) = floor(1 + (point(0) - grid.origin(0)) / cell.dx) * cell.dx;
    cell.topright(1) = floor(1 + (point(1) - grid.origin(1)) / cell.dx) * cell.dx;

    // std::cout<<"Bottom Left = \n"<<cell.bottemleft<<std::endl;
    // std::cout<<"Bottom Right = \n"<<cell.bottemright<<std::endl;
    // std::cout<<"Top Left = \n"<<cell.topleft<<std::endl;
    // std::cout<<"Top Right = \n"<<cell.topright<<std::endl;

    return cell;
}


/* Given a point in space, return its phi value */
double getPhi(Eigen::Vector2d point)
{
    /* Main Logic */
    double phi_p, phi00, phi01, phi10, phi11, phi_x0, phi_x1;
    double a, b;
    CELL cell;

    cell = getCell(point);

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
Eigen::Vector2d getNablaPhi(Eigen::Vector2d point)
{
    Eigen::Vector2d np;
    double a, b, phi00;

    CELL cell = getCell(point);
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

    np(0) = (getPhi(point) - phi00) / a; // TODO
    np(1) = (getPhi(point) - phi00) / b;
    
    return np;
}

void movePointToBoundary(Eigen::Vector2d point)
{
    double currentDistance = DBL_MAX;   // Max Double from float.h
    Eigen::Vector2d nabla_phi;

    // Stop until phi is close to 0
    while(std::abs(currentDistance) > EPSILON)
    {
        std::cout<<"Current point location: \n"<< point <<std::endl;
        
        currentDistance = getPhi(point);
        std::cout<<"Current Phi: \n"<< currentDistance <<std::endl;

        nabla_phi = getNablaPhi(point);
        std::cout << "nable_phi = \n" << nabla_phi << std::endl;


        point(0) = point(0) - currentDistance * nabla_phi(0); //TODo
        point(1) = point(0) - currentDistance * nabla_phi(1);
    }

    std::cout<<"Final point location: \n"<< point <<std::endl;

}

int main(int argc, char const *argv[])
{
    Eigen::Vector2d point;
    point << 2, 2;

    double phi_p = getPhi(point);
    if(phi_p == 0) std::cout << "Given point is on the boundary" << std::endl;
    if(phi_p < 0) std::cout << "Given point is in the boundary" << std::endl;
    if(phi_p > 0) std::cout << "Given point is out of boundary" << std::endl;

    movePointToBoundary(point);
    // 0.071067
    return 0;
}
