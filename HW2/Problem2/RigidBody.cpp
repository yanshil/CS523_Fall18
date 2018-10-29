#define THRESHOLD 0.2
#include "RigidBody.h"

using namespace rigidbody;

double sdf_box(Vector3d input, Vector3d origin, Vector3d width);

RigidBody::RigidBody()
{
	mass = 0;
	one_over_mass = 0;

	v << 0, 0, 0;
	omega << 0, 0, 0;

	x << 0, 0, 0;
	n << 1, 1, 1;
	q = Eigen::Quaterniond{Eigen::AngleAxisd{0, n}};
	P << 0, 0, 0;
	L << 0, 0, 0;

	Ibody = Matrix3d::Zero();
	Ibodyinv = Matrix3d::Zero();

	I = Matrix3d::Zero();
	Iinv = Matrix3d::Zero();
	R = q.normalized().toRotationMatrix();

	force << 0, 0, 0;
	torque << 0, 0, 0;

	vertices = NULL;
	particleNum = 0;

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

	if (particleNum != 0)
	{
		one_over_mass = 1.0 / mass;
		I = R * Ibody * R.transpose();
		Ibodyinv = Ibody.inverse();
		P = v * mass;
		L = I * omega;
		Iinv = I.inverse();
	}
	else
	{
		// Solid Wall with 1/mass = 0 and Iinv = Zeros()
		one_over_mass = 0;
		Iinv = Matrix3d::Zero();
		n << 0, 0, 1;	// Normal
		x << 0, 0, -5;
	}
}

void RigidBody::setCenterofMass(Vector3d x)
{
	this->x = x;
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

	// std::cout << "omega=" << omega << std::endl;

	v = P / mass;
	x += v * h;

	// Update q
	Quaterniond omega_q;
	omega_q.w() = 0;
	omega_q.vec() = omega;
	omega_q.normalize();
	Quaterniond qdot_2 = omega_q * q;

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

void RigidBody::modelWall()
{
	particleNum = 0;
}