//!#####################################################################
//! \file FluidQuantity.h
//!#####################################################################
// Class FluidQuantity
//######################################################################

#ifndef __FluidQuantity__
#define __FluidQuantity__

#include <nova/Tools/Grids/Grid.h>
#include <nova/Tools/Utilities/Range_Iterator.h>

namespace Nova{
template <typename T, int d>
class FluidQuantity
{
    using T_INDEX = Vector<int, d>;
    using TV = Vector<T, d>;

  public:
    T *Phi;     // Array of num_cell
    T *Phi_new; // Array of num_cell
    int axis;        // -1 for Scalar. 0, 1, 2 to indicate faces axis
    int number_of_ghost_cells;
    Grid<T, d> *grid;

    T_INDEX storing_counts;

    FluidQuantity()
    {
        grid = nullptr;
        std::cout << "Wrong Constructor!" << std::endl;
    }
    FluidQuantity(Grid<T, d> &grid, int axis, int number_of_ghost_cells);
    ~FluidQuantity();

    void fill(T content);

    // ********************************************

    /*!
     * 1D Linear Interpolate etween a and b for x in (0, 1)
     */
    T linter(T a, T b, T x)
    {
        return (1.0 - x) * a + x * b;
    }

    bool Inside_Domain(const T_INDEX &index)
    {
        T_INDEX diff = storing_counts - index;

        for (int i = 0; i < d; i++)
        {
            if (diff(i) < 0)
            {
                return false;
            }
        }
        return true;
    }

    bool Inside_Domain(const TV &location)
    {
        TV tmp_location = TV(location);
        // Clamp to domain if not in domain
        T_INDEX tmp_index = (*grid).Cell(tmp_location, number_of_ghost_cells);

        T_INDEX diff = storing_counts - tmp_index;

        for (int i = 0; i < d; i++)
        {
            if (diff(i) < 0)
            {
                return false;
            }
        }
        return true;
    }

    /*!
    * Auxiliary Function: Get exact offset for Column-based 1D array
    * return (z * xSize * ySize) + (y * xSize) + x;
    */
    // TODO:
    int index2offset(const T_INDEX &index)
    {
        // Becuase index in the grid start from (1,1)...
        T_INDEX tmp_index = index - T_INDEX(1);

        int os = tmp_index[1] * storing_counts[0] + tmp_index[0];
        if (d == 3)
            os += tmp_index[2] * storing_counts[0] * storing_counts[1];
        return os;
    }

    T_INDEX offset2index(const int os)
    {
        // 3D: os = z * m * n + y * m + x
        // 2D: os = y * m + x
        T_INDEX tmp_index = T_INDEX();

        // x <- os mod m
        tmp_index[0] = os % storing_counts[0];

        if (d == 2)
            // y <- (os - x) / m
            tmp_index[1] = (os - tmp_index[0]) / storing_counts[1];
        else
        {
            // y <- (os - x) mod n
            tmp_index[1] = (os - tmp_index[0]) % storing_counts[1];

            // z <- (os - x - y * m) / (m*n)
            tmp_index[2] = (os - tmp_index[0] - tmp_index[1] * storing_counts[0]) / storing_counts[0] / storing_counts[1];
        }

        // Becuase index in the grid start from (1,1)...
        tmp_index += T_INDEX(1);
        return tmp_index;
    }

    /*!
     * Linear Interpolator for TV{i, j, k} on grid
     * Coordinates will be clamped to lie in simulation domain
    */
    T linter(TV &location);

    //----------Auxiliary Function--------------------
    void printPhi()
    {
        std::cout << "Phi____: ";

        for (int i = 0; i < storing_counts.Product(); i++)
        {
            if (i % storing_counts[0] == 0)
            {
                std::cout << "\n";
            }
            std::cout << Phi[i] << ", ";
        }

        std::cout << "\n================================" << std::endl;
    }

    void printPhi_new()
    {
        std::cout << "Phi_new: ";

        for (int i = 0; i < storing_counts.Product(); i++)
        {
            if (i % storing_counts[0] == 0)
            {
                std::cout << "\n";
            }
            std::cout << Phi_new[i] << ", ";
        }

        std::cout << "\n================================" << std::endl;
    }

    T_INDEX Next_Cell(const int axis, const T_INDEX &index)
    {
        T_INDEX shifted_index(index);
        shifted_index(axis) += 1;

        return shifted_index;
    }
    // TODO
    T_INDEX Previous_Cell(const int axis, const T_INDEX &index)
    {
        T_INDEX shifted_index(index);
        shifted_index(axis) -= 1;

        return shifted_index;
    }

    /*!
     * Read and Write access in the quantity field for Phi
     */
    T &modify_at(const T_INDEX &index)
    {
        if (!Inside_Domain(index))
        {
            std::cout << "index = " << index << std::endl;
            // Raise exception
            throw std::runtime_error("Try to write at an Out_of_domain area");
        }
        // if (!(*grid).Inside_Domain(index))
        // {
        //     std::cout << "index = " << index << std::endl;
        //     // Raise exception
        //     throw std::runtime_error("Try to write at an Out_of_domain area");
        // }
        return Phi[index2offset(index)];
    }

    /*!
     * Read access for the quantity field of Phi
     */
    T at(const T_INDEX &index)
    {
        // e.g. For n = 4
        // Density: 4 * 4 * 4;      velocityU: 5 * 4 * 4;   velocity V: 4 * 5 * 4;  velocity W: 4 * 4 * 5
        // this->Inside_Domain are checking if is in (1,1,1) to (5, 4, 4) [For U]
        // (*grid).Inside_Domain(index) is checking if is in (1,1,1) to (4,4,4)
        // (*grid).Inside_Domain(index, num_of_ghost_cell = 1) is checking if is in (0,0,0) to (5,5,5)

        //----------------- Velocity /Scalar Field Domain (d/u/v domain) ------------------
        if (!Inside_Domain(index))
        {
            TV location = (*grid).Center(index);
            T_INDEX clamped_index = (*grid).Clamp_To_Cell(location, number_of_ghost_cells);
            return Phi[index2offset(clamped_index)];
        }

        //------------------ Grid Domain (grid domain (ghost_cell = 0)) ----------------------------
        // if (!(*grid).Inside_Domain(index))
        // {
        //     TV location = (*grid).Center(index);
        //     T_INDEX clamped_index = (*grid).Clamp_To_Cell(location, number_of_ghost_cells);

        //     return Phi[index2offset(clamped_index)];
        // }

        return Phi[index2offset(index)];
    }

    /*!
     * RGB color range 0-1 access in the
     */
    T rgb_at(const T_INDEX &index)
    {
        return std::max(std::min(1.0 - at(index), 1.0), 0.0);
    }

    /*!
     * Read & Write access in the quantity field for Phi_new
     */
    T &new_at(const T_INDEX &index)
    {
        if (!Inside_Domain(index))
        {
            // Raise exception
            throw std::runtime_error("Try to write at an Out_of_domain area");
        }
        // if (!(*grid).Inside_Domain(index))
        // {
        //     // Raise exception
        //     throw std::runtime_error("Try to write at an Out_of_domain area");
        // }
        return Phi_new[index2offset(index)];
    }

    /* Compute Velocity */
    TV computeVelocity(const T_INDEX &index, FluidQuantity *velocityField[d]);
    TV Clamp_To_Domain(const TV &location);
    // TV traceBack(const TV &location, T timestep, TV &velocity);

    void flip()
    {
        std::swap(Phi, Phi_new);
    }

    /* Advection */
    void advect(const T_INDEX &index, T timestep, FluidQuantity *velocityField[d]);
    void advect(T timestep, FluidQuantity *velocityField[d]);

    /* */
    // void updateVelocity(const T_INDEX &index, T timestep, T np)
    // {
    //     new_at(index) -= np * timestep;
    // }
};
}


#endif