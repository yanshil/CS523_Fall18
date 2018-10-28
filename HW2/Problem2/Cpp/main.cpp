#include <GL/glut.h>
#include "RigidBody.h"

using namespace rigidbody;

// GLOBAL Rigibody Instance
RigidBody rb;

void copyToV(RigidBody* rb);

void init();
void drawBox();
void display();
void keyboard(unsigned char key, int, int);
void timer(int);


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

// Light Control
GLfloat light_diffuse[] = {1.0, 0.3, 0.5, 1.0};  /* Red diffuse light. */
GLfloat light_position[] = {10.0, 10.0, 10.0, 10.0}; /* Infinite light location. */


int main(int argc, char *argv[])
{

    rb.modelCube();
    rb.setVelocity(Vector3d(0, 0, 10));
    rb.setOmega(Vector3d(0.05, 0.02, 0.01));
    rb.initialize();

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

void init(void)
{

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

void drawBox()
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
    copyToV(&rb);   // Copy the coordinates to global v

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    drawBox();
    glFlush();
    glutSwapBuffers();

    // Calculate new position and orientation of vertices
    rb.find_all_collisions();
    rb.update(0.01);
}

// Keyboard callback function ( called on keyboard event handling )
void keyboard(unsigned char key, int x, int y)
{
    if (key == 'q' || key == 'Q')
        exit(EXIT_SUCCESS);
}


void timer(int)
{
    /* update animation */
    glutPostRedisplay();
    glutTimerFunc(1000.0/60.0, timer, 0);
}

void copyToV(RigidBody* rb) {
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 3; j++) {
            v[i][j] = rb->vertices[i].ri[j];
        }
    }
}