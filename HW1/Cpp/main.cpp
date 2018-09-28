#include <Eigen/Geometry>
#include <iostream>
#include <GL/glut.h>

// g++ main.cpp -lglut -lGLU -lGL -o rotate
// g++ main.cpp -lglut -lGLU -lGL -I /home/yanshi/Libraries/eigen/ -std=gnu++11


/* Reference:

    OpenGL Cube Example: Cube.c
    https://www.opengl.org/archives/resources/code/samples/glut_examples/examples/examples.html
*/

using namespace std;
using namespace Eigen;

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

    // Constructor
    RigidBody()
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
    }

    ~RigidBody()
    {
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

        for (int i = 0; i < particleNum; i++)
        {
            vertices[i].mi = 1.0 / (double)particleNum;
            vertices[i].fi << 0, 0, -9.8 * vertices[i].mi;
        }

        // TODO: Improve locality!
        vertices[0].ri << -1, -1, 1;
        vertices[1].ri <<  -1, -1, -1;
        vertices[2].ri << -1, 1, -1;
        vertices[3].ri << -1, 1, 1;
        vertices[4].ri << 1, -1, 1;
        vertices[5].ri << 1, -1, -1;
        vertices[6].ri << 1, 1, -1;
        vertices[7].ri <<  1, 1, 1;

        // TODO: Copy value to r0 rather than duplicate initialize?
        // Reduce redundancy!
        vertices[0].r0 << -1, -1, 1;
        vertices[1].r0 <<  -1, -1, -1;
        vertices[2].r0 << -1, 1, -1;
        vertices[3].r0 << -1, 1, 1;
        vertices[4].r0 << 1, -1, 1;
        vertices[5].r0 << 1, -1, -1;
        vertices[6].r0 << 1, 1, -1;
        vertices[7].r0 <<  1, 1, 1;

    }

    /* Initialize Rigid Body Element  */
    void initiateRB()
    {

        v << 0, 0, 10;
        omega << 0.05, 0.02, 0.01;

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

    /* Update RB Element */
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
        
        // Update Torque
        for (int i = 0; i < particleNum; i++)
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


    // Auxiluary printing function
    void printStatus()
    {
        std::cout << "x=" << x << std::endl;
        std::cout << "v=" << v << std::endl;
        std::cout << "q.w() = " << q.w() << ", q.vec() = " << std::endl
                  << q.vec() << std::endl;
        std::cout << "R=" << std::endl
                  << R << std::endl;

        std::cout << "I=" << std::endl
                  << I << std::endl;
        std::cout << "L=" << L << std::endl;
        std::cout << "Force=" << force << std::endl;
        std::cout << "torque=" << std::endl
                  << torque << std::endl;
        std::cout << "Ibody=" << std::endl
                  << Ibody << std::endl;
    }

    void printQuaternion()
    {
        std::cout << "q.w() = " << q.w() << ", q.vec() = " << std::endl
                  << q.vec() << std::endl;
    }

    void printX()
    {
        std::cout << "x=" << std::endl
                  << x << std::endl;
    }

    void printRotationMatrix()
    {
        std::cout << "R=" << std::endl
                  << R << std::endl;
    }

    void printVertices()
    {

        for (int i = 0; i < 8; i++)
        {
            /* code */
            std::cout << "v[" << i << "] =" << std::endl
                      << vertices[i].ri << std::endl;
        }
    }

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

// GLOBAL Rigibody Instance
RigidBody rb;

// Display callback ------------------------------------------------------------
GLfloat n[6][3] = {/* Normals for the 6 faces of a cube. */
                   {-1.0, 0.0, 0.0},
                   {0.0, 1.0, 0.0},
                   {1.0, 0.0, 0.0},
                   {0.0, -1.0, 0.0},
                   {0.0, 0.0, 1.0},
                   {0.0, 0.0, -1.0}};
GLint faces[6][4] = {/* Vertex indices for the 6 faces of a cube. */
                     {0, 1, 2, 3},
                     {3, 2, 6, 7},
                     {7, 6, 5, 4},
                     {4, 5, 1, 0},
                     {5, 6, 2, 1},
                     {7, 4, 0, 3}};
GLfloat v[8][3]; /* Will be filled in with X,Y,Z vertexes. */

void copyToV(RigidBody* rb);

void drawBox(void)
{
    
    int i;

    for (i = 0; i < 6; i++)
    {
        glBegin(GL_QUADS);
        glNormal3fv(&n[i][0]);
        glVertex3fv(&v[faces[i][0]][0]);
        glVertex3fv(&v[faces[i][1]][0]);
        glVertex3fv(&v[faces[i][2]][0]);
        glVertex3fv(&v[faces[i][3]][0]);
        glEnd();
    }
}

void display()
{

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


    // Calculate new position and orientation of vertices
    rb.update(0.01);
    copyToV(&rb);   // Copy the coordinates to global v


    drawBox();
    glFlush();
    glutSwapBuffers();
}

// Keyboard callback function ( called on keyboard event handling )
void keyboard(unsigned char key, int x, int y)
{
    if (key == 'q' || key == 'Q')
        exit(EXIT_SUCCESS);
}


// Light Control
GLfloat light_diffuse[] = {1.0, 0.3, 0.5, 1.0};  /* Red diffuse light. */
GLfloat light_position[] = {10.0, 10.0, 10.0, 10.0}; /* Infinite light location. */

void init(void)
{
    /* Setup cube vertex data. */
    v[0][0] = v[1][0] = v[2][0] = v[3][0] = -1;
    v[4][0] = v[5][0] = v[6][0] = v[7][0] = 1;
    v[0][1] = v[1][1] = v[4][1] = v[5][1] = -1;
    v[2][1] = v[3][1] = v[6][1] = v[7][1] = 1;
    v[0][2] = v[3][2] = v[4][2] = v[7][2] = 1;
    v[1][2] = v[2][2] = v[5][2] = v[6][2] = -1;

    /* Enable a single OpenGL light. */
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHTING);

    /* Use depth buffering for hidden surface elimination. */
    glEnable(GL_DEPTH_TEST);

    /* Setup the view of the cube. */
    glMatrixMode(GL_PROJECTION);
    gluPerspective(/* field of view in degree */ 80.0,
                   /* aspect ratio */ 1.0,
                   /* Z near */ 1.0, /* Z far */ 100.0);
    glMatrixMode(GL_MODELVIEW);    
    gluLookAt(10.0, 10.0, 0.0, /* eye is at (10,10,10) */
              0.0, 0.0, 0.0, /* center is at (0,0,0) */
              0.0, 0.0, 1.0); /* up is in positive Z direction */

}


void timer(int)
{
    /* update animation */
    glutPostRedisplay();
    glutTimerFunc(1000.0/60.0, timer, 0);
}

void copyToV(RigidBody* rb)
{
    for (int i = 0; i < 8; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            v[i][j] = rb->vertices[i].ri[j];
        }
    }
}

int main(int argc, char *argv[])
{

    rb.modelCube();
    rb.initiateRB();

    glutInit(&argc, argv);              // Initialize GLUT
    glutInitWindowSize(1024, 1024);
    glutCreateWindow("Rotating Cude"); // Create a window

    glutDisplayFunc(display);   // Register display callback
    glutKeyboardFunc(keyboard); // Register keyboard callback
    init();
    glutTimerFunc(100, timer, 0);
    glutMainLoop(); // Enter main event loop

    return (EXIT_SUCCESS);
}
