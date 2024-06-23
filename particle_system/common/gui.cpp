#include "gui.h"

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
