#include <Eigen/Geometry>
#include <iostream>

using namespace std;
using namespace Eigen;

struct Particle
{
    double mi;
    Vector3d fi;
    Vector3d r0;
    Vector3d ri;    // TODO: Low spatial locality when initiate!!
};

class RigidBody
{
    public:

    double mass;

    Vector3d v;
    Vector3d omega;

    Vector3d x;     // x(t)
    Quaterniond q;  // q(t)
    Vector3d P;     // P(t)
    Vector3d L;     // L(t)

    Matrix3d Ibody;
    Matrix3d Ibodyinv;

    Matrix3d Iinv;
    Matrix3d I;
    Matrix3d R;

    Vector3d force;
    Vector3d torque;

    Particle *vertices;
    int particleNum;

    //Particles vertices[8];

    // Constructor
    RigidBody()
    {
        x << 0,0,0;
        v << 0,0,0;
        omega << 0,0,0;

        q = Eigen::Quaterniond{Eigen::AngleAxisd{0, Eigen::Vector3d{1, 1, 1}}};
        R = q.normalized().toRotationMatrix();

        mass = 0;
        particleNum = 0;

        force << 0,0,0;
        torque << 0,0,0;

        vertices = NULL;

        Ibody  = Matrix3d::Zero();
    }

    ~RigidBody(){
        delete[] vertices;
        cout << "Called Destructor" << endl;
    }


    /* TODO: Setter */
    void setV(Vector3d velocity);
    void setOmega(Vector3d angular_velocity);


    /* Model Example: Cube*/
    void modelCube()
    {
        particleNum = 8;
        vertices = new Particle[particleNum];
        
        for(int i = 0; i < particleNum; i++)
        {
            vertices[i].mi = 1.0 / (double)particleNum;
            vertices[i].fi << 0,0,-9.8 * vertices[i].mi;

        }
        // TODO: Improve locality!
        vertices[0].ri << 1,1,1;
        vertices[1].ri << 1,-1,1;
        vertices[2].ri << 1,-1,-1;
        vertices[3].ri << 1,1,-1;
        vertices[4].ri << -1,1,1;
        vertices[5].ri << -1,-1,1;
        vertices[6].ri << -1,-1,-1;
        vertices[7].ri << -1,1,-1;


        // TODO: Copy value to r0 rather than duplicate initialize
        vertices[0].r0 << 1,1,1;
        vertices[1].r0 << 1,-1,1;
        vertices[2].r0 << 1,-1,-1;
        vertices[3].r0 << 1,1,-1;
        vertices[4].r0 << -1,1,1;
        vertices[5].r0 << -1,-1,1;
        vertices[6].r0 << -1,-1,-1;
        vertices[7].r0 << -1,1,-1;
    }


    /*  */
    void initiateRB()
    {


        v << 0,0,10;
        omega << 0.05, 0.02, 0.01;
        
        for(int i = 0; i < particleNum; i++)
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

    void update(double h)
    {
        static int frame = 0;
        frame++;

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
        P += force * h;
        I = R * Ibody * R.transpose();
        // TODO: update torque here?
        for (int i =0; i < particleNum; i++)
        {
            torque += vertices[i].ri.cross(vertices[i].fi);
        }
        L += torque * h;
        omega = I.inverse() * L;

        // Update ri!
        for (int j = 0; j < particleNum; j++)
        {
            vertices[j].ri = R * vertices[j].r0 + x;
        }

    }

    void allocateParticle(int pnum)
    {
        particleNum = pnum;
        vertices = new Particle[particleNum];
        cout << "Allocated " << particleNum << " particlaes" << endl;
    }

    void printStatus()
    {
        std::cout << "x=" << x << std::endl;
        std::cout << "v=" << v << std::endl;
        std::cout << "q.w() = " << q.w() << ", q.vec() = " << std::endl << q.vec() << std::endl;
        std::cout << "R=" << std::endl << R << std::endl;
        
        std::cout << "I=" << std::endl << I << std::endl;
        std::cout << "L=" << L << std::endl;
        std::cout << "Force=" << force << std::endl;
        std::cout << "torque=" << std::endl << torque << std::endl;
        std::cout << "Ibody=" << std::endl << Ibody << std::endl;
    }

    void printQuaternion()
    {
          std::cout << "q.w() = " << q.w() << ", q.vec() = " << std::endl << q.vec() << std::endl;
    }

    void printRotationMatrix()
    {
        std::cout << "R=" << std::endl << R << std::endl;
    }
    
    public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};


int main(int argc, char const *argv[])
{
    RigidBody rb;
    rb.modelCube();
    rb.initiateRB();
    rb.printStatus();
    for(int i =0; i< 10; i++)
    {
        rb.update(0.01);
        rb.printRotationMatrix();
    }
    std::cout << "Run the simulation and display!" << endl;

    return 0;
}
