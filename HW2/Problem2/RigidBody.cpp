#define THRESHOLD 0.2
#include "RigidBody.h"

using namespace rigidbody;

double sdf_box(Vector3d input, Vector3d origin, Vector3d width);

RigidBody::RigidBody()
{
	x << 0, 0, 0;
	v << 0, 0, 0;
	omega << 0, 0, 0;

	q = Eigen::Quaterniond{Eigen::AngleAxisd{0, Eigen::Vector3d{1, 1, 1}}};
	R = q.normalized().toRotationMatrix();

	mass = 0;
	particleNum = 0;

	force << 0, 0, 0;
	torque << 0, 0, 0;

	vertices = NULL;

	Ibody = Matrix3d::Zero();
	Ibodyinv = Matrix3d::Zero();

	I = Matrix3d::Zero();
	Iinv = Matrix3d::Zero();

} // End constuctor

RigidBody::~RigidBody()
{
	delete[] vertices;
	cout << "Called Destructor" << endl;
}

void RigidBody::setVelocity(Vector3d velocity)
{
	this->v = velocity;
}

void RigidBody::setOmega(Vector3d omega)
{
	this->omega = omega;
}

void RigidBody::initialize()
{
	for (int i = 0; i < particleNum; i++)
	{
		mass += vertices[i].mi;
		force += vertices[i].fi;
		torque += vertices[i].ri.cross(vertices[i].fi);

		Ibody += vertices[i].mi * (vertices[i].r0.transpose() * vertices[i].r0 * Matrix3d::Identity() - vertices[i].r0 * vertices[i].r0.transpose());
	}

	Ibodyinv = Ibody.inverse();

	P = v * mass;
	I = R * Ibody * R.transpose();
	L = I * omega;
}

void RigidBody::update(double h)
{
	static int frame = 0;
	frame++;

	// // Update Torque
	// // TODO: Bug! Torque should be 0,0,0, in free-fall. But here return value other than 0

	// torque << 0, 0, 0;
	// for (int i = 0; i < particleNum; i++)
	// {
	// 	torque += vertices[i].ri.cross(vertices[i].fi);
	// 	// std::cout<<"vertices[i].ri"<<vertices[i].ri<<std::endl;
	// }
	// std::cout << "torque=" << torque << std::endl;

	P += force * h;
	L += torque * h;
	omega = Iinv * L;

	v = P / mass;
	x += v * h;

	// Update q
	Quaterniond omega_q;
	omega_q.w() = 0;
	omega_q.vec() = omega;
	omega_q.normalize();
	Quaterniond qdot_2 = omega_q * q;
	/*
	    2nd time!!!
	    Wrote 1/2, but well, 1/2 = 0.
	    Should be 1.0 /2 or (double)1/2
	*/
	q.w() = q.w() + h * 0.5 * (qdot_2.w());
	q.x() = q.x() + h * 0.5 * (qdot_2.x());
	q.y() = q.y() + h * 0.5 * (qdot_2.y());
	q.z() = q.z() + h * 0.5 * (qdot_2.z());
	q.normalize();

	R = q.normalized().toRotationMatrix();
	
	I = R * Ibody * R.transpose();
	Iinv = R * Ibodyinv * R.transpose();

	// Update ri!
	for (int j = 0; j < particleNum; j++)
	{
		vertices[j].ri = R * vertices[j].r0 + x;
	}
}

/* Model Example: Cube*/
void RigidBody::modelCube()
{
	particleNum = 8;
	vertices = new Particle[particleNum];

	for (int i = 0; i < particleNum; i++)
	{
		vertices[i].mi = 1.0 / (double)particleNum;
		vertices[i].fi << 0, 0, -9.8 * vertices[i].mi;
	}

	// TODO: Improve locality!
	vertices[0].ri << -1, -1, 1;
	vertices[1].ri << -1, -1, -1;
	vertices[2].ri << -1, 1, -1;
	vertices[3].ri << -1, 1, 1;
	vertices[4].ri << 1, -1, 1;
	vertices[5].ri << 1, -1, -1;
	vertices[6].ri << 1, 1, -1;
	vertices[7].ri << 1, 1, 1;

	// TODO: Copy value to r0 rather than duplicate initialize?
	// Reduce redundancy!
	vertices[0].r0 << -1, -1, 1;
	vertices[1].r0 << -1, -1, -1;
	vertices[2].r0 << -1, 1, -1;
	vertices[3].r0 << -1, 1, 1;
	vertices[4].r0 << 1, -1, 1;
	vertices[5].r0 << 1, -1, -1;
	vertices[6].r0 << 1, 1, -1;
	vertices[7].r0 << 1, 1, 1;
}

Eigen::Vector3d RigidBody::pt_velocity(Eigen::Vector3d point)
{
	return v + omega.cross(point - x);
}

void RigidBody::collision(double epsilon, Eigen::Vector3d point)
{
	// Get velocity of collide poin
	Eigen::Vector3d padot, n, ra;

	padot = pt_velocity(point);

	n << 0, 0, 1;
	ra = point - x;

	double vrel, numerator;
	// double vrel = n * (padot - pbdot);
	vrel = n.dot(padot);
	numerator = -(1 + epsilon) * vrel;

	double term1, term2, term3, term4, j;
	term1 = 1 / mass;
	term2 = 0;
	term3 = n.dot(Iinv * (ra.cross(n)).cross(ra));
	term4 = 0;

	/* Compute the impulse */
	j = numerator / (term1 + term2 + term3 + term4);
	Eigen::Vector3d impulse_forse;
	impulse_forse = j * n;

	std::cout << "v = " << v << std::endl;
	std::cout << "omega = " << omega << std::endl;

	std::cout << "vrel = " << vrel << std::endl;
	std::cout << "impulse = " << impulse_forse << std::endl;

	/* Apply the impulse to the bodies */
	P += impulse_forse;
	L += ra.cross(impulse_forse);

	/* Recompute auxiliary variables */
	v = P / mass;
	omega = Iinv * L;
}

bool RigidBody::colliding_with_ground(Eigen::Vector3d point, double groundz)
{
	// Find collision point: Min z vector
	Eigen::Vector3d n, padot;
	double vrel;

	/* code */
	n << 0, 0, 1;
	padot = pt_velocity(point);
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

void RigidBody::find_all_collisions()
{

	bool had_collision;
	double epsilon = 0.5;
	double groundz = -5;

	do
	{
		had_collision = false;

		for (int i = 0; i < 7; i++)
		{
			Vector3d estimated_collision_point;
			estimated_collision_point << vertices[i].ri(0), vertices[i].ri(1), groundz;

			/* Sensitive detect */
			if (sdf_box(estimated_collision_point, x, Vector3d(1, 1, 1)) <= 0)
			{
				if (colliding_with_ground(vertices[i].ri, groundz))
				{
					std::cout << "----------------" << std::endl;
					std::cout << "collition point:" << vertices[i].ri << std::endl;
					std::cout << "estimated_collision_point:" << estimated_collision_point << std::endl;
					collision(epsilon, estimated_collision_point);
					// collision(epsilon, vertices[i].ri);
					had_collision = true;

					// Tell the solver we had a collison
				}
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