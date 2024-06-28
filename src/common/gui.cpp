#include "gui.h"
#include "common/mmath.h"

void loadShaders(bool is_reload)
{
	GLuint shader;

	shader = labhelper::loadShaderProgram("./simple.vert", "./simple.frag", is_reload);
	if (shader != 0)
	{
		simpleShaderProgram = shader;
	}

	shader = labhelper::loadShaderProgram("./shading.vert", "./shading.frag", is_reload);
	if (shader != 0)
	{
		shaderProgram = shader;
	}

	shader = labhelper::loadShaderProgram("./line.vert", "./line.frag", is_reload);
	if (shader != 0)
	{
		debugLineProgram = shader;
	}
	shader = labhelper::loadShaderProgram("./smoke.vert", "./smoke.frag", is_reload);
	if (shader != 0)
	{
		smokeProgram = shader;
	}
}

///////////////////////////////////////////////////////////////////////////////
/// This function is to hold the general GUI logic
///////////////////////////////////////////////////////////////////////////////
void ControlPanel()
{
	ImGui::Begin("Control panel");
	ImGui::SliderFloat("Light intensity", &g_point_light_intensity, 0.0f, 100000.0f);
	ImGui::SliderFloat("Smoke factor", &g_smoke_factor, 0.0f, 100.0f);
	ImGui::SliderFloat("Env temperature", &g_env_temp, 0.0f, 1000.0f);
	ImGui::SliderFloat("Alpha", &g_alpha, 0, 100);
	ImGui::SliderFloat("Beta", &g_beta, 0, 10);
	ImGui::SliderFloat("Vort Eps", &g_vort_eps, 0, 100);
	ImGui::SliderFloat("Decay Factor", &g_decay_factor, 0.9, 1);
	ImGui::SliderFloat("Dt", &g_dt, 0.001, 0.05);

	if (ImGui::Button("Case 0: Empty"))
	{
		g_case_id = Demo::SMOKE_EMPTY;
	}

	if (ImGui::Button("Case 1: Sphere"))
	{
		g_case_id = Demo::SMOKE_SPHERE;
	}

	if (ImGui::Button("Case 2: Cube"))
	{
		g_case_id = Demo::SMOKE_CUBE;
	}

	if (ImGui::Button("Case 3: Tree"))
	{
		g_case_id = Demo::SMOKE_TREE;
	}

	// set camera position
	if (ImGui::Button("View Front"))
	{
		cameraPosition = glm::vec3(30.f, 5.0f, 5.f);
		cameraDirection = normalize(-glm::vec3(1.f, 0.f, 0.f));
	}
	// set camera position
	if (ImGui::Button("View Top"))
	{
		cameraPosition = glm::vec3(5.f, 30.f, 5.f);
		cameraDirection = normalize(-glm::vec3(0.0f, 1.f, 0.1f));
	}
	// set camera position
	if (ImGui::Button("View Side"))
	{
		cameraPosition = glm::vec3(5.f, 5.0f, 30.f);
		cameraDirection = normalize(-glm::vec3(0.0f, 0.f, 1.f));
	}

	// reset button, to reset simulator
	if (ImGui::Button("Reset"))
	{
		g_simulator_rest = true;
	}

	ImGui::End();
}

void gui()
{
	// ----------------- Set variables --------------------------
	ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate,
				ImGui::GetIO().Framerate);
	{
		std::lock_guard<std::mutex> lock(g_sim_lock);
		ImGui::Text("Simulator Info: %s", simulator_info.c_str());
	}

	// ----------------------------------------------------------
	labhelper::perf::drawEventsWindow();

	ControlPanel();
}

// Mouse input
glm::ivec2 g_prevMouseCoords = {-1, -1};
bool g_isMouseDragging = false;

///////////////////////////////////////////////////////////////////////////////
/// This function is used to update the scene according to user input
///////////////////////////////////////////////////////////////////////////////
bool handleEvents(void)
{
	// check events (keyboard among other)
	SDL_Event event;
	bool quitEvent = false;
	while (SDL_PollEvent(&event))
	{
		labhelper::processEvent(&event);

		if (event.type == SDL_QUIT || (event.type == SDL_KEYUP && event.key.keysym.sym == SDLK_ESCAPE))
		{
			quitEvent = true;
		}
		if (event.type == SDL_KEYUP && event.key.keysym.sym == SDLK_g)
		{
			if (labhelper::isGUIvisible())
			{
				labhelper::hideGUI();
			}
			else
			{
				labhelper::showGUI();
			}
		}
		if (event.type == SDL_MOUSEBUTTONDOWN && event.button.button == SDL_BUTTON_LEFT && (!labhelper::isGUIvisible() || !ImGui::GetIO().WantCaptureMouse))
		{
			g_isMouseDragging = true;
			int x;
			int y;
			SDL_GetMouseState(&x, &y);
			g_prevMouseCoords.x = x;
			g_prevMouseCoords.y = y;
		}

		if (!(SDL_GetMouseState(NULL, NULL) & SDL_BUTTON(SDL_BUTTON_LEFT)))
		{
			g_isMouseDragging = false;
		}

		if (event.type == SDL_MOUSEMOTION && g_isMouseDragging)
		{
			// More info at https://wiki.libsdl.org/SDL_MouseMotionEvent
			int delta_x = event.motion.x - g_prevMouseCoords.x;
			int delta_y = event.motion.y - g_prevMouseCoords.y;
			float rotationSpeed = 0.1f;
			glm::mat4 yaw = rotate(rotationSpeed * deltaTime * -delta_x, worldUp);
			glm::mat4 pitch = rotate(rotationSpeed * deltaTime * -delta_y,
									 normalize(cross(cameraDirection, worldUp)));
			cameraDirection = glm::vec3(pitch * yaw * glm::vec4(cameraDirection, 0.0f));
			g_prevMouseCoords.x = event.motion.x;
			g_prevMouseCoords.y = event.motion.y;
		}
	}

	// check keyboard state (which keys are still pressed)
	const uint8_t *state = SDL_GetKeyboardState(nullptr);
	glm::vec3 cameraRight = cross(cameraDirection, worldUp);

	if (state[SDL_SCANCODE_W])
	{
		cameraPosition += cameraSpeed * deltaTime * cameraDirection;
	}
	if (state[SDL_SCANCODE_S])
	{
		cameraPosition -= cameraSpeed * deltaTime * cameraDirection;
	}
	if (state[SDL_SCANCODE_A])
	{
		cameraPosition -= cameraSpeed * deltaTime * cameraRight;
	}
	if (state[SDL_SCANCODE_D])
	{
		cameraPosition += cameraSpeed * deltaTime * cameraRight;
	}
	if (state[SDL_SCANCODE_Q])
	{
		cameraPosition -= cameraSpeed * deltaTime * worldUp;
	}
	if (state[SDL_SCANCODE_E])
	{
		cameraPosition += cameraSpeed * deltaTime * worldUp;
	}
	return quitEvent;
}

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

void drawCoordinate(const glm::mat4 &projectionViewModelMatrix)
{
	// draw coordinate system
	glDisable(GL_DEPTH_TEST);
	//
	drawLine(glm::vec3(0.0f),
			 glm::vec3(10.0, 0.0, 0.0),
			 projectionViewModelMatrix * glm::mat4(1.0),
			 glm::vec3(1.0f, 0.0f, 0.0f));
	drawLine(glm::vec3(0.0f),
			 glm::vec3(0.0, 10.0, 0.0),
			 projectionViewModelMatrix * glm::mat4(1.0),
			 glm::vec3(0.0f, 1.0f, 0.0f));
	drawLine(glm::vec3(0.0f),
			 glm::vec3(0.0, 0.0, 10.0),
			 projectionViewModelMatrix * glm::mat4(1.0),
			 glm::vec3(0.0f, 0.0f, 1.0f));
	//
	glEnable(GL_DEPTH_TEST);
}

void drawVexel(const std::array<bool, SIZE> &occupied_voxels, const glm::mat4 &viewMatrix,
			   const glm::mat4 &projectionMatrix, const labhelper::Model *model)
{
	glUseProgram(simpleShaderProgram);

	// scale up the cube
	// draw the occupied voxels
	for (int i = 0; i < Nx; i++)
	{
		for (int j = 0; j < Ny; j++)
		{
			for (int k = 0; k < Nz; k++)
			{
				if (occupied_voxels[ACC3D(i, j, k, Ny, Nx)])
				{

					glm::mat4 modelMatrix = glm::mat4(1.0f);
					modelMatrix = glm::scale(glm::vec3(0.2, 0.2, 0.2)) * modelMatrix;
					modelMatrix = glm::translate(
									  glm::vec3(
										  10.f / Nx * i,
										  10.f / Nx * j,
										  10.f / Nx * k)) *
								  modelMatrix;
					labhelper::setUniformSlow(simpleShaderProgram, "modelViewProjectionMatrix",
											  projectionMatrix * viewMatrix * modelMatrix);
					labhelper::render(model);
				}
			}
		}
	}
}
