
//!#####################################################################
//! \file CG_System.h
//!#####################################################################
// Class CG_System
//######################################################################
#ifndef __CG_System__
#define __CG_System__

#include "CG_Storage.h"
#include "CG_Vector.h"
#include <nova/Tools/Krylov_Solvers/Krylov_System_Base.h>

namespace Nova
{
template <class T, int d>
class CG_System : public Krylov_System_Base<T>
{
  private:
    using TV = Vector<T, d>;
    using Base = Krylov_System_Base<T>;
    using Vector_Base = Krylov_Vector_Base<T>;

    CG_storage<T,d> &storage;

  public:
    CG_System(const CG_storage<T,d> &storage)
        : Base(false, false), storage(storage)
    {
    }

    ~CG_System()
    {
    }

    void Multiply(const Vector_Base &v, Vector_Base &result) const
    {
        // Array<TV> &v_array = CG_Vector<T, d>::CG_Array(const_cast<Vector_Base &>(v));
        // Array<TV> &result_array = CG_Vector<T, d>::CG_Array(result);

        // result_array.Fill(TV());

        // for (int iy = 0; iy < n; iy++)
        // {
        //     for (int ix = 0; ix < m; ix++)
        //     {
        //         int idx = iy * m + ix;
        //         T sum = A.diag(idx) * result_array(idx);

        //         if (ix > 0)
        //         { // if u.previous valid
        //             sum += A.plusi(idx - 1) * result_array(idx - 1);
        //         }

        //         if (iy > 0)
        //         {
        //             sum += A.plusj(idx - m) * result_array(idx - m);
        //         }

        //         if (ix < m - 1)
        //         {
        //             sum += A.plusi(idx) * result_array(idx + 1);
        //         }

        //         if (iy < n - 1)
        //         {
        //             sum += A.plusj(idx) * result_array(idx + m);
        //         }

        //         result_array(idx) = sum;
        //     }
        // }
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
        for (int i = 0; i < storage.size; ++i)
            // result=std::max(result,x_array(i).Norm());
            result = std::max(result, x_array(i).Norm());
        return result;
    }

    void Project(Vector_Base &x) const {}

    void Set_Boundary_Conditions(Vector_Base &v) const {}
    void Project_Nullspace(Vector_Base &x) const {}
    void Apply_Preconditioner(const Krylov_Vector_Base<T>& r,Krylov_Vector_Base<T>& z) const {}
};

} // namespace Nova

#endif
