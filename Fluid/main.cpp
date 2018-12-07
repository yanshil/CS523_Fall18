#include <GL/glut.h>
#include <iostream>

#include "FluidSolver.h"
#define CELLCOUNTS 128

using namespace Nova;

/*-------- Global typedef & varaible------------*/
enum
{
    d = 2
};
using T = double;
using T_INDEX = Vector<int, d>;
using TV = Vector<T, d>;

//----------------------------------------------------

T timestep = 0.01;
T density = 1;

int iterations = 0;

// Range<T,d> range(TV(), TV(1));
// FluidSimulator_Grid<T, d> grid(T_INDEX(CELLCOUNTS), range);
FluidSimulator_Grid<T, d> grid(T_INDEX(CELLCOUNTS), Range<T, d>::Unit_Box());
FluidSolver<T, d> *solver = new FluidSolver<T,d>(grid, 1);
// FluidSolver<T,d> *solver = new FluidSolver<T,d>(CELLCOUNTS, CELLCOUNTS, density);

int iteration = 1;
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

            T_INDEX index{x + 1, y + 1};
            GLfloat color = solver->toRGB(index);

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
    // std::cout << "Frame " << iteration << std::endl;
    iteration++;
    //--------------------------------------------------------------

    glClearColor(0.0f, 0.0f, 0.0f, 1.0f); // Set background color to black and opaque
    glClear(GL_COLOR_BUFFER_BIT);         // Clear the color buffer (background)

    drawGrid();

    glFlush(); // Render now

    //-------------------------------------
    solver->addInflow(T_INDEX{58,2}, T_INDEX{64,7}, -1, 1);
    solver->addInflow(T_INDEX{58,2}, T_INDEX{64,7}, 1, 1);
    solver->update(timestep);
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
    solver->output = true;
    solver->initialize();

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