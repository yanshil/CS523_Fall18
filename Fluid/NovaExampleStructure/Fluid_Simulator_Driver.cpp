
//!#####################################################################
//! \file Fluid_Simulator_Driver.cpp
//!#####################################################################
#include <nova/Tools/Krylov_Solvers/Conjugate_Gradient.h>
#include "CG_System.h"
#include "Fluid_Simulator_Driver.h"

using namespace Nova;
//######################################################################
// Constructor
//######################################################################
template<class T,int d> Fluid_Simulator_Driver<T,d>::
Fluid_Simulator_Driver(Fluid_Solver<T, d> &solver_input)
    :Base(solver_input), solver(solver_input)
{}
//######################################################################
// Initialize
//######################################################################
template<class T,int d> void Fluid_Simulator_Driver<T,d>::
Initialize()
{

    Base::Initialize();     // virtual?!
    // solver.Initialize();

    // if(!solver.restart) solver.Initialize();
    // else{solver.Initialize_Position_Based_State();
    //     solver.Read_Output_Files(solver.restart_frame);}

    solver.Initialize_Auxiliary_Structures();
}
//######################################################################
// Advance_One_Newton_Iteration
//######################################################################
template<class T,int d> void Fluid_Simulator_Driver<T,d>::
Advance_One_Newton_Iteration(const T target_time,const T dt)
{
    solver.Set_Kinematic_Positions(time+dt,solver.simulation_mesh->points);
    solver.Interpolate_Embedded_Values(solver.embedded_surface->points,solver.simulation_mesh->points);
    solver.Set_Kinematic_Velocities(time+dt,solver.volume_velocities);
    solver.Interpolate_Embedded_Values(solver.surface_velocities,solver.volume_velocities);

    const size_t number_of_particles=solver.simulation_mesh->points.size();
    Array<TV> rhs(number_of_particles);
    solver.position_based_state->Update_Position_Based_State(solver.simulation_mesh->points);
    solver.Add_Elastic_Forces(solver.simulation_mesh->points,rhs);
    solver.Add_Damping_Forces(solver.volume_velocities,rhs);
    solver.Add_External_Forces(rhs);
    for(size_t i=0;i<number_of_particles;++i) rhs(i)+=(solver.mass(i)*(Vp(i)-solver.volume_velocities(i)))/dt;
    solver.Clear_Values_Of_Kinematic_Particles(rhs);

    Array<TV> delta_X(number_of_particles);
    Array<TV> temp_q(number_of_particles),temp_s(number_of_particles),temp_r(number_of_particles),temp_k(number_of_particles),temp_z(number_of_particles);
    CG_Vector<T,d> cg_x(delta_X),cg_b(rhs),cg_q(temp_q),cg_s(temp_s),cg_r(temp_r),cg_k(temp_k),cg_z(temp_z);

    CG_System<T,d> cg_system(solver,dt);
    Conjugate_Gradient<T> cg;

    T b_norm=cg_system.Convergence_Norm(cg_b);
    Log::cout<<"Norm: "<<b_norm<<std::endl;

    cg.print_residuals=true;
    cg.print_diagnostics=true;
    cg.restart_iterations=solver.cg_restart_iterations;

    // solve
    cg.Solve(cg_system,cg_x,cg_b,cg_q,cg_s,cg_r,cg_k,cg_z,solver.cg_tolerance,0,solver.cg_iterations);

    // update position and velocity
    solver.Clear_Values_Of_Kinematic_Particles(delta_X);
    solver.simulation_mesh->points+=delta_X;solver.volume_velocities+=delta_X/dt;
    solver.Interpolate_Embedded_Values(solver.embedded_surface->points,solver.simulation_mesh->points);
    solver.Interpolate_Embedded_Values(solver.surface_velocities,solver.volume_velocities);
}
//######################################################################
// Advance_To_Target_Time
//######################################################################
template<class T,int d> void Fluid_Simulator_Driver<T,d>::
Advance_To_Target_Time(const T target_time)
{
    bool done=false;
    for(int substep=1;!done;substep++){
        Log::Scope scope("SUBSTEP","substep "+std::to_string(substep));
        T dt=Compute_Dt(time,target_time);
        solver<T,d>::Clamp_Time_Step_With_Target_Time(time,target_time,dt,done);

        Vp=solver.volume_velocities;
        for(size_t i=0;i<solver.simulation_mesh->points.size();++i) solver.simulation_mesh->points(i)+=dt*Vp(i);
        for(int iteration=0;iteration<solver.newton_iterations;++iteration){
            Log::cout<<"Newton Iteration: "<<iteration+1<<"/"<<solver.newton_iterations<<std::endl;
            Advance_One_Newton_Iteration(time,dt);}

        if(!done) solver.Write_Substep("END Substep",substep,0);
        time+=dt;}
}
//######################################################################
// Simulate_To_Frame
//######################################################################
template<class T,int d> void Fluid_Simulator_Driver<T,d>::
Simulate_To_Frame(const int target_frame)
{
    solver.frame_title="Frame "+std::to_string(solver.current_frame);
    if(!solver.restart) Write_Output_Files(solver.current_frame);

    while(solver.current_frame<target_frame){
        Log::Scope scope("FRAME","Frame "+std::to_string(++solver.current_frame));

        Advance_To_Target_Time(solver.Time_At_Frame(solver.current_frame));

        solver.frame_title="Frame "+std::to_string(solver.current_frame);
        Write_Output_Files(++solver.output_number);

        *(solver.output)<<"TIME = "<<time<<std::endl;}
}
//######################################################################
template class Nova::Fluid_Simulator_Driver<float,2>;
template class Nova::Fluid_Simulator_Driver<float,3>;
#ifdef COMPILE_WITH_DOUBLE_SUPPORT
template class Nova::Fluid_Simulator_Driver<double,2>;
template class Nova::Fluid_Simulator_Driver<double,3>;
#endif