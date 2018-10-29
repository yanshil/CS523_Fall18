#include "Contact.h"
#define THRESHOLD 0.01

Eigen::Vector3d pt_velocity(RigidBody *rb, Eigen::Vector3d p)
{
    return rb->v + rb->omega.cross(p - rb->x);
}

bool colliding(Contact *c)
{
    Eigen::Vector3d padot = pt_velocity(c->a, c->p);
    Eigen::Vector3d pbdot = pt_velocity(c->b, c->p);
    double vrel = c->n.dot(padot);

    if (vrel > THRESHOLD)
    {
        return false;
    }
    if (vrel > -THRESHOLD)
    {
        return false;
    }
    else if( c->p[2] <= c->b->x[2] )
    {
        std::cout<<"c->p="<<c->p<<std::endl;
        std::cout<<"c->b->x(2)="<<c->b->x[2]<<std::endl;
        std::cout<<"vrel="<<vrel<<std::endl;
        return true;
    }
    else
        return false;
}

void collision(Contact *c, double epsilon)
{
    Eigen::Vector3d padot, pbdot, n, ra, rb;
    padot = pt_velocity(c->a, c->p);
    pbdot = pt_velocity(c->b, c->p);
    n = c->n;
    ra = c->p - c->a->x;
    rb = c->p - c->b->x;

    double vrel = c->n.dot(padot);
    double numerator = -(1 + epsilon) * vrel;

    double term1, term2, term3, term4, j;
    term1 = c->a->one_over_mass;
    term2 = c->b->one_over_mass;
    term3 = n.dot(c->a->Iinv * (ra.cross(n)).cross(ra));
    term4 = n.dot(c->b->Iinv * (rb.cross(n)).cross(rb));

    /* Compute the impulse */
    j = numerator / (term1 + term2 + term3 + term4);
    Eigen::Vector3d impulse_forse = j * n;

    std::cout << "impulse = " << impulse_forse << std::endl;

    /* Apply the impulse to the bodies */
    c->a->P += impulse_forse;
    c->b->P -= impulse_forse;
    c->a->L += ra.cross(impulse_forse);
    c->b->L -= rb.cross(impulse_forse);

    /* Recompute auxiliary variables */
    c->a->v = c->a->P / c->a->mass;
    c->b->v = c->b->P / c->b->mass;

    c->a->omega = c->a->Iinv * c->a->L;
    c->b->omega = c->b->Iinv * c->b->L;
}

void find_all_collisions(Contact contacts[], int ncontacts)
{
    bool had_collision;
    double epsilon = 0.5;

    do
    {
        had_collision = false;

        for (int i = 0; i < ncontacts; i++)
        {
            contacts[i].p = contacts[i].a->vertices[i].ri;
            contacts[i].n = contacts[i].b->n;

            if (colliding(&contacts[i]))
            {
                collision(&contacts[i], epsilon);
                had_collision = true;

                /* Tell the solver we had a collision */
                // TODO: Find the exact collision time and compute update.
            }
        }

    } while (had_collision);
}

double sdf_box(Vector3d input, Vector3d origin, Vector3d width)
{

    Vector3d dTmp;

    for (int i = 0; i < dTmp.size(); i++)
    {
        dTmp[i] = std::abs(input[i] - origin[i]) - width[i];
    }

    // If in the box, this will be the SDF.
    double dTmpMax = dTmp.maxCoeff();

    double dTmpLength = 0;
    for (int i = 0; i < dTmp.size(); i++)
    {
        if (dTmp(i) < 0)
            dTmp(i) = 0;
        dTmpLength += std::pow(dTmp(i), 2.0);
    }

    dTmpLength = std::sqrt(dTmpLength);

    return std::min(dTmpMax, 0.0) + dTmpLength;
}