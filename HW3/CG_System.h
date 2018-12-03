//!#####################################################################
//! \file CG_System.h
//!#####################################################################
// Class CG_System
//######################################################################
#ifndef __CG_System__
#define __CG_System__

#include <nova/Tools/Krylov_Solvers/Krylov_System_Base.h>
#include "CG_Vector.h"
#include "CG_Storage.h"

namespace Nova
{
template <class T, int d>
class CG_System : public Krylov_System_Base<T>
{
    using TV = Vector<T, d>;
    using Base = Krylov_System_Base<T>;
    using Vector_Base = Krylov_Vector_Base<T>;

    CG_Storage<T, d> &storage;

  public:
    CG_System(CG_Storage<T, d> &storage_input)
        : Base(false, false), storage(storage_input)
    {
    }

    void Multiply(const Vector_Base &v, Vector_Base &result) const
    {
        // Array<TV>& v_array                  = CG_Vector<T,d>::CG_Array(const_cast<Vector_Base&>(v));
        // Array<TV>& result_array             = CG_Vector<T,d>::CG_Array(result);

        // Array<TV> G(v_array.size()),K(v_array.size());
        // result_array.Fill(TV());

        // storage.Add_Damping_Forces(v_array,G);
        // storage.Add_Force_Differentials(v_array,K);

        // for(size_t i=0;i<storage.simulation_mesh->points.size();++i)
        //     result_array(i)=(storage.mass(i)*v_array(i))*one_over_dt_squared-one_over_dt*G(i)-K(i);
    }

    double Inner_Product(const Vector_Base &x, const Vector_Base &y) const
    {
        const Array<TV> &x_array = CG_Vector<T, d>::CG_Array(x);
        const Array<TV> &y_array = CG_Vector<T, d>::CG_Array(y);

        double result = (T)0.;
        for (size_t i = 0; i < storage.size; ++i)
            result += x_array(i).Dot_Product(y_array(i));
        return result;
    }

    T Convergence_Norm(const Vector_Base &x) const
    {
        const Array<TV> &x_array = CG_Vector<T, d>::CG_Array(x);

        T result = (T)0.;
        for (size_t i = 0; i < storage.size; ++i)
            result = std::max(result, x_array(i).Norm());
        return result;
    }

    void Project(Vector_Base &x) const
    {
    }

    void Set_Boundary_Conditions(Vector_Base &v) const {}
    void Project_Nullspace(Vector_Base &x) const {}
};
} // namespace Nova
#endif