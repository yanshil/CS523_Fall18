//!#####################################################################
//! \file CG_Driver.h
//!#####################################################################
// Class CG_Driver
//######################################################################
#ifndef __CG_Driver__
#define __CG_Driver__

#include <nova/Tools/Krylov_Solvers/Conjugate_Gradient.h>
#include <nova/Tools/Utilities/Driver.h>
#include "CG_Storage.h"
#include "CG_System.h"

namespace Nova
{
template <class T, int d>
class CG_Driver : public Driver<T, d>
{
    using TV = Vector<T, d>;
    using Base = Driver<T, d>;
    int size;

  public:
    CG_Storage<T, d> &storage;

    CG_Driver(CG_Storage<T, d> &storage_input)
        :storage(storage_input)
    {
        int size = storage.size;

    }
    ~CG_Driver() {}

    void Execute()
    {
        Array<TV> rhs(size);
        Array<TV> delta_X(size);
        Array<TV> temp_q(size), temp_s(size), temp_r(size), temp_k(size), temp_z(size);
        CG_Vector<T, d> cg_x(delta_X), cg_b(rhs), cg_q(temp_q), cg_s(temp_s), cg_r(temp_r), cg_k(temp_k), cg_z(temp_z);

        CG_System<T, d> cg_system(storage);
        Conjugate_Gradient<T> cg;

        T b_norm = cg_system.Convergence_Norm(cg_b);
        Log::cout << "Norm: " << b_norm << std::endl;

        cg.print_residuals = true;
        cg.print_diagnostics = true;
        cg.restart_iterations = storage.cg_restart_iterations;

        // solve
        cg.Solve(cg_system, cg_x, cg_b, cg_q, cg_s, cg_r, cg_k, cg_z, storage.cg_tolerance, 0, storage.cg_iterations);
    }
};
} // namespace Nova
#endif