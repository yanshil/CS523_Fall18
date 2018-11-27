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
  using T_INDEX = Vector<int, d>;
  using TV = Vector<T, d>;

public:
  int number_of_ghost_cells;
  int size;

  FluidQuantity<T, d> *velocityField[d];
  FluidQuantity<T, d> *density_field;
  Grid<T, d> *grid;

  T_INDEX storing_counts;
  T density;
  T *_rhs;
  T *_p;

  FluidSolver();
  FluidSolver(Grid<T, d> &grid, T density, int number_of_ghost_cells);
  ~FluidSolver();

  T getRGBcolorDensity(T_INDEX &index);

  /* Advection */
  void advection(T timestep);

  /* Calculate the RHS of Poisson Equation */
  void calculateDivergence();

  /* Projection with CG */
  void pressure_solution_Jacobi();
  void projection(int limit);
  void projection(int limit, T timestep = 0.12);
  void Project(int limit);
  /* ================== */
  void SetDivBoundary();
  void SetPressureBoundary();
  void SetVelocityBoundary();
  void SetVelocityBoundary_as0();

  /* Update velocity with pressure */
  void updateVelocity(T timestep);
  /* Make result visulizable */
  void flip();

  //-----------------------------------------------

  void initialize();
  /* UPDATE */
  void addInflow(const T_INDEX &index, const T density, const TV &velocity);

  void update(T timestep);

  //-----------------------------------------------
  T_INDEX Next_Cell(const int axis, const T_INDEX &index);
  T_INDEX Previous_Cell(const int axis, const T_INDEX &index);
  int index2offset(const T_INDEX &index);
  T_INDEX offset2index(const int os);

  void projection_INCSAMPLE(int limit, T timestep);
};

} // namespace Nova

#endif