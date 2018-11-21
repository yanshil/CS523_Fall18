//!#####################################################################
//! \file Fluid_Solver.cpp
//!#####################################################################
#include <nova/Tools/Utilities/File_Utilities.h>
#include <nova/Tools/Parsing/Parse_Args.h>
#include "Fluid_Solver.h"
using namespace Nova;
//######################################################################
// Initialize
//######################################################################
template<class T,int d> void Fluid_Solver<T,d>::
Initialize()
{
    Initialize_Simulation_Mesh();
    // Initialize_Position_Based_State();
    // Initialize_Embedded_Mesh();
    // Initialize_Embedding_Map_And_Weights();
}

//######################################################################
template class Nova::Fluid_Solver<float,2>;
template class Nova::Fluid_Solver<float,3>;
#ifdef COMPILE_WITH_DOUBLE_SUPPORT
template class Nova::Fluid_Solver<double,2>;
template class Nova::Fluid_Solver<double,3>;
#endif