#include <Eigen/Geometry>
#include <iostream>

#include "RigidBody.h"
#define THRESHOLD 0.001

using namespace rigidbody;

double sdf_box(Vector3d input, Vector3d origin, Vector3d width);



Eigen::Vector3d pt_velocity(RigidBody rb, Eigen::Vector3d point)
{
	return rb.v + rb.omega.cross(point - rb.x);
}

void collision(RigidBody rb, double epsilon, Eigen::Vector3d point)
{
	// Get velocity of collide poin
	Eigen::Vector3d padot, n, ra;
	padot = pt_velocity(rb, point);
	n << 0, 0, 1;
	ra = point - rb.x;

	double vrel, numerator;
	// double vrel = n * (padot - pbdot);
	vrel = n.dot(padot);
	numerator = -(1 + epsilon) * vrel;

	double term1, term2, term3, term4, j;
	term1 = 1 / rb.mass;
	term2 = 0;
	term3 = n.dot(rb.Iinv * ( ra.cross(n) ).cross(ra) );
	term4 = 0;

	/* Compute the impulse */
	j = numerator / (term1 + term2 + term3 + term4);
	Eigen::Vector3d impulse_forse;
	impulse_forse = j * n;
	std::cout<<"impulse = "<<impulse_forse<<std::endl;

	/* Apply the impulse to the bodies */
	rb.P += impulse_forse;
	rb.L += ra.cross(impulse_forse);

	/* Recompute auxiliary variables */
	rb.v = rb.P / rb.mass;
	rb.omega = rb.Iinv * rb.L;
}

bool colliding_with_ground(RigidBody rb, Eigen::Vector3d point, double groundz)
{
	// Find collision point: Min z vector
	Eigen::Vector3d n, padot;
	double vrel;

	/* code */
	n << 0, 0, 1;
	padot = pt_velocity(rb, point);
	vrel = n.dot(padot);

	if (vrel > THRESHOLD)
	{
		return false;
	}
	if (vrel > -THRESHOLD)
	{
		return false;
	}
	else
		return true;
}

void find_all_collisions(RigidBody rb)
{

	bool had_collision;
	double epsilon = 0.5;
	double groundz = -5;

	do
	{
		had_collision = false;

		for (int i = 0; i < 7; i++)
		{
			Eigen::Vector3d point, width, corresponding_point_on_ground;

			point = rb.vertices[i].ri;
			width << 1, 1, 1;
			corresponding_point_on_ground << point(0), point(1), groundz;

			if (sdf_box(corresponding_point_on_ground, rb.x, width) <= 0)
			{
				if (colliding_with_ground(rb, point, groundz))
				{
					collision(rb, epsilon, point);
					had_collision = true;

					/* Tell the solver we had a collision */
					// TODO: Find the exact collision time and compute update.
				}
			}
		}

	} while (had_collision);
}


double sdf_box(Vector3d input, Vector3d origin, Vector3d width)
{
	//   vec3 d = abs(p) - b;
	//   return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));

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