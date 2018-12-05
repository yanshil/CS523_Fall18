#include "CG_Driver.h"
#define CELLSIZE 4

using namespace Nova;

enum
{
    d = 2
};
using T = double;
using TV = Vector<T, d>;
using T1 = Vector<T, 1>;
using T_INDEX = Vector<int, d>;

void printArray(Array<T1> array, int size)
{
    for (int i = 0; i < size; i++)
    {
        T1 tmp = array(i);
        std::cout << tmp(0) << ", ";

        if ((i + 1) % CELLSIZE == 0)
        {
            std::cout << "\n";
        }
    }
    std::cout << std::endl;
}

void test1(CG_System<T, d> &cg_system, Array<T1> &result, int size);
void test3(CG_System<T, d> &cg_system, CG_Storage<T, d> &storage, Array<T1> &result);

int main(int argc, char const *argv[])
{
    Range<T, d> range(TV(-0.5), TV(0.5));
    Grid<T, d> *grid = new Grid<T, d>(T_INDEX(CELLSIZE), range);
    CG_Storage<T, d> storage(*grid);

    storage.setting();
    storage.initialize();
    storage.calculateA();
    storage.calculateRHS();

    // CG_Driver<T,d> driver(storage);
    // driver.Execute();
    CG_System<T, d> cg_system(storage);

    Array<T1> result(storage.size);

    // ============Test 1=======================
    // Should Print A
    // test1(cg_system, result, storage.size);

    // ============Test 2=======================

    // ============Test 3=======================
    // Should print RHS
    storage.calculateTrueValue();

    printf("Calculated RHS: \n");
    test3(cg_system, storage, result);

    // print RHS
    printf("True RHS: \n");
    for (int i = 0; i < storage.size; i++)
    {
        std::cout << storage._rhs(i) << ", ";
    }
    std::cout << std::endl;

    // print true value
    printf("True Phi: \n");
    for (int i = 0; i < storage.size; i++)
    {
        printf("%0.4f, ", storage._trueValue(i).Max());
        
        if ((i + 1) %  storage.m == 0) {
            std::cout <<"\n";
        }        
    }
    std::cout << std::endl;

    // ============print A=======================

    // storage.printAFromStorage();

    // ============Aucxiliary Info=======================

    // std::cout << "one_over_dX_square = " << storage.one_over_dX_square << std::endl;
    // std::cout << "grid->one_over_dX = " << grid->one_over_dX << std::endl;

    return 0;
}

// ============ Unit Test Functions =======================

void unitVector(Array<T1> &e, int j, int size)
{
    e.resize(size);
    e.Fill(T1());
    e(j) = 1;
}

void test1(CG_System<T, d> &cg_system, Array<T1> &result, int size)
{
    CG_Vector<T, d> cg_result(result);

    for (int j = 0; j < size; j++)
    {
        for (int i = 0; i < size; i++)
        {
            int idx = j * size + i;

            Array<T1> e1, e2;
            unitVector(e1, i, size);
            unitVector(e2, j, size);

            CG_Vector<T, d> cg_e1(e1), cg_e2(e2);

            cg_system.Multiply(cg_e1, cg_result);
            double re = cg_system.Inner_Product(cg_e2, cg_result);

            printf("%3.0f", re);

            if ((idx + 1) % size == 0)
            {
                std::cout << "\n";
            }
            else
            {
                std::cout << ",";
            }
        }
    }
    std::cout << std::endl;
}

void test2()
{
    // Random x and y

    // x^T A y = y^T A x
    // if A is symmetric, compute x^T A y an y^T A x and check if equal
}

void test3(CG_System<T, d> &cg_system, CG_Storage<T, d> &storage, Array<T1> &result)
{
    // Multiply A and trueValue: A * phi = rhs;
    CG_Vector<T, d> cg_result(result), cg_tv(storage._trueValue);
    // should get RHS = 4;
    cg_system.Multiply(cg_tv, cg_result);

    for (int idx = 0; idx < storage.size; idx++)
    {
        printf("%3.0f", result(idx).Max());
        if ((idx + 1) % storage.size == 0)
        {
            std::cout << "\n";
        }
        else
        {
            std::cout << ",";
        }
    }
    std::cout << std::endl;
}