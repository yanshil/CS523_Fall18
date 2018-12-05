//!#####################################################################
//! \file FluidSolver.h
// https://github.com/OrionQuest/Nova_Examples/blob/master/embedded_deformables/Embedded_Deformables_Example.h
//!#####################################################################
// Class FluidSolver
//######################################################################

#ifndef __FluidSolver__
#define __FluidSolver__

#include "FluidQuantity.h"

namespace Nova
{
template <typename T, int d>
class FluidSolver
{
private:
  using TV = Vector<T, d>;
  using T_INDEX = Vector<int, d>;

  // Density Field
  FluidQuantity<T, d> *_d;
  // Velocity Field
  FluidQuantity<T, d> *_v[d];

  // The inertia grid: m, n, k;
  FluidSimulator_Grid<T, d> *grid;

  // Inertia grid m * n;
  int size;

  T rho;
  T hx;

  T *_rhs;
  T *_p;

  void calculateRHS();
  void setBoundaryCondition();
  void project_GS(int limit, T timestep, bool output = true);
  void applyPressure(T timestep);
  void advection(T timestep);
  void flip();

  int index2offset(const T_INDEX &index) const
  {
    return grid->index2offset(index, grid->counts);
  }

  T_INDEX offset2index(const int os) const
  {
    return grid->offset2index(os, grid->counts);
  }

public:
  FluidSolver()
  {
    std::cout << "??" << std::endl;
  }

  FluidSolver(FluidSimulator_Grid<T, d> &grid, T rho)
      : grid(&grid), rho(rho)
  {

    this->_d = new FluidQuantity<T, d>(grid, -1);

    for (int axis = 0; axis < d; axis++)
      this->_v[axis] = new FluidQuantity<T, d>(grid, axis);

    this->size = grid.counts.Product();
    this->_rhs = new T[size];
    this->_p = new T[size];

    this->hx = grid.hx;
    memset(_p, 0, size * sizeof(double));
  }

  ~FluidSolver()
  {
    delete[] _d;

    for (int axis = 0; axis < d; axis++)
      delete[] _v[axis];
    
    delete[] _rhs;
    delete[] _p;
  }

  T toRGB(T_INDEX &index) const
  {
    int idx = index2offset(index);
    return std::max(std::min(1.0 - _d->Phi()[idx], 1.0), 0.0);
  }
  void update(T timestep);
  void addInflow(const T_INDEX &min_corner, const T_INDEX &max_corner,
                 int axis, T input_value);

  ////===============================================================
  void printPressure()
  {
    for (int i = 0; i < size; i++)
    {

      printf("%.3f", _p[i]);
      if ((i + 1) % grid->counts[0] == 0)
      {
        printf("\n");
      }
      else
      {
        printf(",");
      }
    }
  }
  
  void unitTest_tentProjection(T timestep)
  {
    // RHS
    memset(_rhs, 0, size * sizeof(T));
    _rhs[grid->index2offset(T_INDEX{6,6}, grid->counts)] = 1;
    project_GS(1000, timestep, false);
    applyPressure(timestep);
    printPressure();
  }

  ////===============================================================
};
} // namespace Nova
#endif