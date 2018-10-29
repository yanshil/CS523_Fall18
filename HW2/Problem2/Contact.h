#include "RigidBody.h"

using namespace rigidbody;

struct Contact
{
    RigidBody *a, *b;
    Eigen::Vector3d p, n;
    // bool vf;    /* TRUE if vertex/face contact */
};

// int Ncontacts;
// Contact *Contacts;

Eigen::Vector3d pt_velocity(RigidBody *rb, Eigen::Vector3d p);
bool colliding(Contact *c);
void collision(Contact *c, double epsilon);
void find_all_collisions(Contact contacts[], int ncontacts);
double sdf_box(Vector3d input, Vector3d origin, Vector3d width);

