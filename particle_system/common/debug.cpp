
#include "debug.h"

#include "common/globalvar.h"

void drawLine(const glm::vec3 start_pos, const glm::vec3 end_pos,
			  const glm::mat4 &projectionViewModelMatrix,
			  const glm::vec3 color)
{
	GLuint currentShaderProgram = debugLineProgram;
	glUseProgram(currentShaderProgram);
	labhelper::setUniformSlow(currentShaderProgram, "modelViewProjectionMatrix",
							  projectionViewModelMatrix);
	labhelper::setUniformSlow(currentShaderProgram, "line_color", color);

	// draw a line
	static GLuint buffer = 0;
	static GLuint vertexArrayObject = 0;
	static int nofVertices = 2;
	glm::vec3 positions[] = {start_pos, end_pos};
	// do this initialization first time the function is called...
	if (vertexArrayObject == 0)
	{
		glGenVertexArrays(1, &vertexArrayObject);
		glGenBuffers(1, &buffer);
		glBindBuffer(GL_ARRAY_BUFFER, buffer);
		glBufferData(GL_ARRAY_BUFFER, sizeof(positions), positions, GL_STREAM_DRAW);
		CHECK_GL_ERROR();

		// Now attach buffer to vertex array object.
		glBindVertexArray(vertexArrayObject);
		glVertexAttribPointer(0, 3, GL_FLOAT, false, 0, 0);
		glEnableVertexAttribArray(0);
		CHECK_GL_ERROR();
	}
	glBindBuffer(GL_ARRAY_BUFFER, buffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(positions), positions, GL_STREAM_DRAW);

	glBindVertexArray(vertexArrayObject);
	glDrawArrays(GL_LINES, 0, nofVertices);
	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
}