#include <GL/glut.h>
#define SCREEN 100
#define GRID 10
// g++ display.cpp -lGL -lglut

/* Handler for window-repaint event. Call back when the window first appears and
   whenever the window needs to be re-painted. */
void display()
{
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f); // Set background color to black and opaque
    glClear(GL_COLOR_BUFFER_BIT);         // Clear the color buffer (background)

    int quadCount = 16;
    float quadSize = 2.0f / static_cast<float>(quadCount);

    // Draw a Red 1x1 Square centered at origin
    glBegin(GL_QUADS);           // Each set of 4 vertices form a quad
    glColor3f(1.0f, 1.0f, 1.0f); // Red
    bool isWhite = true;

    for (int x = 0; x < quadCount; ++x)
    {
        float xPos = -1.0f + x * quadSize;

        for (int y = 0; y < quadCount; ++y)
        {
            
            if ((x+y) % 2 == 0) {
                glColor3f(1.0f, 1.0f, 1.0f); // White
            }
            else
            {
                glColor3f(0.0f, 0.0f, 0.0f); // Black
            }
            
            float yPos = -1.0f + y * quadSize;

            glVertex2f(xPos, yPos);
            glVertex2f(xPos + quadSize, yPos);
            glVertex2f(xPos + quadSize, yPos + quadSize);
            glVertex2f(xPos, yPos + quadSize);
        }
    }

    glEnd();

    glFlush(); // Render now
}

/* Main function: GLUT runs as a console application starting at main()  */
int main(int argc, char **argv)
{
    glutInit(&argc, argv);                 // Initialize GLUT
    glutCreateWindow("OpenGL Setup Test"); // Create a window with the given title
    glutInitWindowSize(400, 400);        // Set the window's initial width & height
    glutInitWindowPosition(160, 160);      // Position the window's initial top-left corner
    glutDisplayFunc(display);              // Register display callback handler for window re-paint
    glutMainLoop();                        // Enter the event-processing loop
    return 0;
}
