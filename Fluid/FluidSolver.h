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
  enum Boundary_Condition{Interior, Exteria, Dirichlet};

public:
  int number_of_ghost_cells;
  int size_whole_domain;
  int size_interior_domain;

  // Velocity Field
  FluidQuantity<T, d> *_v[d];
  // Density Field
  FluidQuantity<T, d> *_d;
  // The shared grid
  Grid<T, d> *grid;

  // Domain: Interior Domain + Ghost Domain = Whole Domain
  T_INDEX interior_domain;
  T_INDEX whole_domain;
  T density;
  // flag of Boundary Condition;
  Boundary_Condition *_BCflag;
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
  T_INDEX offset2index(const int os);
  int index2offset(const T_INDEX &index);
};

} // namespace Nova

#endif