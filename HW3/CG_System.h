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
#include <nova/Tools/Log/Log.h>

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
        Array<TV> &v_array = CG_Vector<T, d>::CG_Array(const_cast<Vector_Base &>(v));
        Array<TV> &result_array = CG_Vector<T, d>::CG_Array(result);

        // Array<TV> G(v_array.size()),K(v_array.size());
        result_array.Fill(TV());

        for (int iy = 0; iy < storage.n; iy++)
        {
            for (int ix = 0; ix < storage.m; ix++)
            {
                int idx = iy * storage.m + ix;
                result_array(idx) = storage._Adiag(idx) * v_array(idx);

                if (ix > 0)
                { // if u.previous valid
                    result_array(idx) += storage._Aplusi(idx - 1) * v_array(idx - 1);
                }

                if (iy > 0)
                {
                    result_array(idx) += storage._Aplusj(idx - storage.m) * v_array(idx - storage.m);
                }

                if (ix < storage.m - 1)
                {
                    result_array(idx) += storage._Aplusi(idx) * v_array(idx + 1);
                }

                if (iy < storage.n - 1)
                {
                    result_array(idx) += storage._Aplusj(idx) * v_array(idx + storage.m);
                }
            }
        }
    }

    double Inner_Product(const Vector_Base &x, const Vector_Base &y) const
    {
        const Array<TV> &x_array = CG_Vector<T, d>::CG_Array(x);
        const Array<TV> &y_array = CG_Vector<T, d>::CG_Array(y);

        double result = (T)0.;
        for (size_t i = 0; i < storage.size; ++i)
        {
            // TV x0 = x_array(i);
            // TV y0 = x_array(i);
            // result += x0.Max() * y0.Max();
            result += x_array(i).Dot_Product(y_array(i));
        }
        return result;
    }

    T Convergence_Norm(const Vector_Base &x) const
    {
        const Array<TV> &x_array = CG_Vector<T, d>::CG_Array(x);

        T result = (T)0.;
        for (size_t i = 0; i < storage.size; ++i)
            result = std::max(result, x_array(i).Norm());
            // result = std::max(result, x_array(i).Abs().Max());
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