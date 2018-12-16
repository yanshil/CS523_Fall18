#include <GL/glut.h>
#include <iostream>
#include "Fluid.h"
#define CELLCOUNTS 128

double timestep = 0.01;
double density = 1;
int iteration = 0;

FluidSolver *solver = new FluidSolver(CELLCOUNTS, CELLCOUNTS, density);
//--------------------OpenGL--------------------

void drawGrid()
{
    int quadCount = CELLCOUNTS;
    float quadSize = 2.0f / static_cast<float>(quadCount);

    // Draw a Red 1x1 Square centered at origin
    glBegin(GL_QUADS);           // Each set of 4 vertices form a quad
    glColor3f(1.0f, 1.0f, 1.0f); // White

    for (int x = 0; x < quadCount; x++)
    {

        float xPos = -1.0 + x * quadSize;

        for (int y = 0; y < quadCount; y++)
        {
            float yPos = -1.0 + y * quadSize;
            GLfloat color = solver->toRGB(x, y);
            // printf("Color(%d, %d) = %f\n", x, y, color);

            glColor3f(color, color, color);

            glVertex2f(xPos, yPos);
            glVertex2f(xPos + quadSize, yPos);
            glVertex2f(xPos + quadSize, yPos + quadSize);
            glVertex2f(xPos, yPos + quadSize);
        }
    }

    glEnd();
}

void display()
{
    //--------------------------------------------------------------

    glClearColor(0.0f, 0.0f, 0.0f, 1.0f); // Set background color to black and opaque
    glClear(GL_COLOR_BUFFER_BIT);         // Clear the color buffer (background)

    drawGrid();

    glFlush(); // Render now

    //-------------------------------------
    // // x0, y0, x1, y1, d, u, v
    // // /* Bottom Mid Single inflow */

    // solver->addInflow(62, 1, 66, 4, -1, 1);
    // solver->addInflow(62, 1, 66, 4, 1, 1);

    
    solver->addInflow(57, 1, 64, 7, -1, 1);
    solver->addInflow(57, 1, 64, 7, 1, 1);
    
    solver->update(timestep);
    
    // /* Top Right and circle test*/
    // solver->addInflow(0.8, 0.8, 0.81, 0.81, 1, 0, 0);
    // solver->update(timestep);

    // Nothing initialized

    iteration++;
    // std::cout << "iteration = " << iteration << std::endl;
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
    glutTimerFunc(1000.0 / 60, timer, 0);
}

// -------------------------------------------------

int main(int argc, char **argv)
{
    // solver->test_initialize_CircleVelocityField();
    // solver->test_initialize_circleDensityField();

    glutInit(&argc, argv);                 // Initialize GLUT
    glutCreateWindow("OpenGL Setup Test"); // Create a window with the given title
    glutInitWindowSize(400, 400);          // Set the window's initial width & height
    // glutInitWindowPosition(160, 160);      // Position the window's initial top-left corner
    glutDisplayFunc(display); // Register display callback handler for window re-paint

    glutKeyboardFunc(keyboard); // Register keyboard callback
    glutTimerFunc(100, timer, 0);
    glutMainLoop(); // Enter the event-processing loop

    return 0;
}