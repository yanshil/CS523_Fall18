#include <nova/Tools/Grids/Grid.h>
#include <nova/Tools/Utilities/Range_Iterator.h>

using namespace Nova;

enum
{
    d = 2
};

typedef float T;
typedef Vector<T, d> TV;
typedef Vector<int, d> T_INDEX;

void unitTest_Ghost(const Grid<T, d> &grid, const TV &location, const int ngc)
{
    std::cout << "================== " << "\n";

    T_INDEX clamp1 = grid.Clamp_To_Cell(location, ngc);
    T_INDEX cell1 = grid.Cell(location, ngc);
    bool inside_domain = grid.Inside_Domain(cell1, ngc);

    std::cout << "TV = " << location << "\t";
    std::cout << "number_of_ghost_cell = " << ngc << "\t";

    std::cout << "Clamp to Index = " << clamp1 << "\t";
    std::cout << "Cell Index = " << cell1 << "\t";

    std::cout << "Inside Domain = " << inside_domain << std::endl;

    std::cout << "Cell_Indices Range =";
    std::cout << grid.Cell_Indices(ngc).min_corner << " --- ";
    std::cout << grid.Cell_Indices(ngc).max_corner << "\t";

    // 0,0 - 6,6
    std::cout << "Node_Indices Range = ";
    std::cout << grid.Node_Indices(ngc).min_corner << " --- ";
    std::cout << grid.Node_Indices(ngc).max_corner << std::endl;
}

int main(int argc, char **argv)
{

    T_INDEX counts = T_INDEX(4);
    Grid<T, d> grid(counts, Range<T, d>::Unit_Box());

    // Expected 4,4
    std::cout << "grid.Number_Of_Cells()";
    std::cout << grid.Number_Of_Cells() << std::endl;

    // Expected 5,5
    std::cout << "grid.Number_Of_Nodes()";
    std::cout << grid.Number_Of_Nodes() << std::endl;

    T_INDEX index = T_INDEX{2, 2};

    // Expected [index]Cell.min_corder
    std::cout << "grid.Node(index)";
    std::cout << grid.Node(index) << std::endl;

    // Expected center location for cell.
    std::cout << "grid.Center(index)";
    std::cout << grid.Center(index) << std::endl;

    // Expected Face location for given index and axis
    // 0.375, 0.375
    std::cout << "grid.Face(1, index)";
    std::cout << grid.Face(1, index) << std::endl;

    // Expected Face location for given index and axis
    // 0.25, 0.375
    std::cout << "grid.Face(0, index)";
    std::cout << grid.Face(0, index) << std::endl;

    // 0.25, 0.25 - 0.5, 0.5
    std::cout << "grid.Cell_Domain(index)";
    std::cout << grid.Cell_Domain(index).min_corner << std::endl;
    std::cout << grid.Cell_Domain(index).max_corner << std::endl;

    // 0,0 - 5,5
    std::cout << "grid.Cell_Indices(ngc)";
    std::cout << grid.Cell_Indices(1).min_corner << std::endl;
    std::cout << grid.Cell_Indices(1).max_corner << std::endl;

    // 0,0 - 6,6
    std::cout << "grid.Node_Indices(ngc)";
    std::cout << grid.Node_Indices(1).min_corner << std::endl;
    std::cout << grid.Node_Indices(1).max_corner << std::endl;

    // --------------- Num_of_ghost_cell -----------------------

    unitTest_Ghost(grid, TV(3), 1);
    unitTest_Ghost(grid, TV(3), 2);
    unitTest_Ghost(grid, TV{1.2, 1}, 1);
    unitTest_Ghost(grid, TV{1.2, 1}, 0);
    unitTest_Ghost(grid, TV{0.05, 0.1}, 1);
    unitTest_Ghost(grid, TV{0.05, 0.1}, 2);
    unitTest_Ghost(grid, TV{-1, -3}, 2);

    return 0;
}

/*
 * 1. Cell Domain: From (1,1) To (n,n)
 * 2. If num_ghost_cell != 0, thickened...
 * 3. Clamp into domain (including ghost_cells domain)
 * 
 */
