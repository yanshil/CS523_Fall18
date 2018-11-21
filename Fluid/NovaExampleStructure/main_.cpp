//!#####################################################################
//! \file main.cpp
// https://github.com/OrionQuest/Nova_Examples/blob/master/embedded_deformables/embedded_deformables_2d/main.cpp
//!#####################################################################
#include "Fluid_Simulator_Driver.h"
using namespace Nova;

int main(int argc,char** argv)
{
    enum {d=2};
    typedef float T;typedef Vector<T,d> TV;
    typedef Vector<int,d> T_INDEX;

    Fluid_Solver<T,d> *solver=new Fluid_Solver<T,d>();
    // Parse Arguments for grids
    // solver->Parse(argc,argv);

    Fluid_Simulator_Driver<T,d> driver(*solver);
    driver.Execute_Main_Program();

    delete solver;

    return 0;
}