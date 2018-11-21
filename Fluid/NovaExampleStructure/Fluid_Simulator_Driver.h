//!#####################################################################
//! \file Fluid_Simulator_Driver.h
// https://github.com/OrionQuest/Nova_Examples/blob/master/embedded_deformables/Embedded_Deformables_Driver.h
//!#####################################################################
// Class Fluid_Simulator_Driver
//######################################################################
#ifndef __Fluid_Simulator_Driver__
#define __Fluid_Simulator_Driver__

// CG_System
#include "Fluid_Solver.h"

namespace Nova
{
template <class T, int d>
class Fluid_Simulator_Driver
{
    using TV = Vector<T, d>;

  public:
    Fluid_Solver<T, d> &solver;

    Fluid_Simulator_Driver(Fluid_Solver<T, d> &solver);
    ~Fluid_Simulator_Driver();
    //######################################################################
    void Initialize();
    void Advance_One_Newton_Iteration(const T target_time, const T dt);
    void Advance_To_Target_Time(const T target_time);
    void Simulate_To_Frame(const int frame);
    void Execute_Main_Program();
    //######################################################################
};
} // namespace Nova

#endif