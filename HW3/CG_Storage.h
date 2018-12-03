//!#####################################################################
//! \file CG_Storage.h
//!#####################################################################
// Class CG_Storage
//######################################################################
#ifndef __CG_Storage__
#define __CG_Storage__

#include <nova/Tools/Utilities/Example.h>

namespace Nova
{
template <class T, int d>
class CG_Storage
{
    using TV = Vector<T, d>;
    using T_INDEX = Vector<int, d>;

  public:
    T_INDEX counts;
    int size;
    int newton_iterations, cg_iterations, cg_restart_iterations;
    T cg_tolerance;

    CG_Storage()
    {
        size = counts.Product();
    }

    ~CG_Storage()
    {
    }
};
} // namespace Nova
#endif