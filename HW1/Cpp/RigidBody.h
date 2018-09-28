#include <Eigen/Geometry>
#include <iostream>

using namespace std;
using namespace Eigen;

namespace rigidbody {

struct Particle
{
    double mi;
    Vector3d fi;
    Vector3d r0;
    Vector3d ri; // TODO: Low spatial locality when initiate!!
};

class RigidBody
{
public:

    double mass;

    Vector3d v;
    Vector3d omega;

    Vector3d x;    // x(t)
    Quaterniond q; // q(t)
    Vector3d P;    // P(t)
    Vector3d L;    // L(t)

    Matrix3d Ibody;
    Matrix3d Ibodyinv;

    Matrix3d Iinv;
    Matrix3d I;
    Matrix3d R;

    Vector3d force;
    Vector3d torque;

    Particle *vertices;
    int particleNum;

    RigidBody();
    ~RigidBody();

    void initialize();
    void update(double timestep);
    void modelCube();
};

}