
#ifdef _WIN32
extern "C" _declspec(dllexport) unsigned int NvOptimusEnablement = 0x00000001;
#endif

#include "common/debug.h"
#include "common/globalvar.h"
#include "common/mperf.h"

#include "particle/render.h"
#include "particle/Simulator.h"
#include "particle/CudaRender.h"
#include "particle/Vexelization.h"

#include <thread>
#include <mutex>

double ttime = 0.0;

std::string simulator_info;

std::array<bool, SIZE> occupied_voxels_sphere = generate_vexelized_sphere((int)(10));
std::array<bool, SIZE> occupied_voxels_cube = generate_vexelized_cube((int)(12));
std::array<bool, SIZE> empty_voxels;

float tree_max_length = 0.0f;
std::array<bool, SIZE> occupied_voxels_tree;

int case_id = 0;

// lock for simulator
std::mutex simLock;
bool simulator_rest_trigger = false;
std::array<double, SIZE> density;
double transparency[SIZE];

float g_env_temp = 25.0;

float smoke_factor = 10.f;
int enable_light_tracing = 1;

///////////////////////////////////////////////////////////////////////////////
// Various globals
///////////////////////////////////////////////////////////////////////////////
SDL_Window *g_window = nullptr;
float currentTime = 0.0f;
float previousTime = 0.0f;
float deltaTime = 0.0f;
int windowWidth, windowHeight;
std::vector<FboInfo> FBOList;

// Mouse input
ivec2 g_prevMouseCoords = {-1, -1};
bool g_isMouseDragging = false;

///////////////////////////////////////////////////////////////////////////////
// Camera parameters.
///////////////////////////////////////////////////////////////////////////////
vec3 cameraPosition(45.0f, 45.0f, 45.0f);
vec3 cameraDirection = normalize(vec3(0.0f) - cameraPosition);
float cameraSpeed = 10.f;
vec3 worldUp(0.0f, 1.0f, 0.0f);

///////////////////////////////////////////////////////////////////////////////
// Environment
///////////////////////////////////////////////////////////////////////////////
float environment_multiplier = 1.0f;
GLuint environmentMap;
GLuint irradianceMap;
GLuint reflectionMap;
const std::string envmap_base_name = "001";

///////////////////////////////////////////////////////////////////////////////
// Light source
///////////////////////////////////////////////////////////////////////////////
vec3 lightPosition = vec4(80.0f, 25.0f, 25.0f, 1.0f);
vec3 point_light_color = vec3(1.f, 1.f, 1.f);

float point_light_intensity_multiplier = 15000.0f;

///////////////////////////////////////////////////////////////////////////////
// Models
///////////////////////////////////////////////////////////////////////////////

labhelper::Model *landingpadModel = nullptr;
labhelper::Model *pointLight = nullptr;
labhelper::Model *sphereModel = nullptr;
labhelper::Model *cubeModel = nullptr;
labhelper::Model *generalModel = nullptr;

SmokeRenderer *mmRender = nullptr;

mat4 landingPadModelMatrix;
mat4 testModelMatrix;
mat4 sphereModelMatrix;
mat4 cubeModelMatrix;
mat4 generalModelMatrix;

///////////////////////////////////////////////////////////////////////////////
// Particle system
///////////////////////////////////////////////////////////////////////////////
// Simulator simulator;
std::unique_ptr<Simulator> simulator;
std::unique_ptr<MCUDA::CudaRender> cudaRender;

void loadShaders(bool is_reload)
{
	GLuint shader;

	shader = labhelper::loadShaderProgram("./simple.vert", "./simple.frag", is_reload);
	if (shader != 0)
	{
		simpleShaderProgram = shader;
	}

	shader = labhelper::loadShaderProgram("./background.vert", "./background.frag", is_reload);
	if (shader != 0)
	{
		backgroundProgram = shader;
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
/// This function is called once at the start of the program and never again
///////////////////////////////////////////////////////////////////////////////
void initialize()
{
	ENSURE_INITIALIZE_ONLY_ONCE();

	///////////////////////////////////////////////////////////////////////
	//		Load Shaders
	///////////////////////////////////////////////////////////////////////
	loadShaders(false);

	///////////////////////////////////////////////////////////////////////
	// Load models and set up model matrices
	///////////////////////////////////////////////////////////////////////
	landingpadModel = labhelper::loadModelFromOBJ("../scenes/plane.obj");
	pointLight = labhelper::loadModelFromOBJ("../scenes/point_light.obj");
	sphereModel = labhelper::loadModelFromOBJ("../scenes/particle.obj");
	cubeModel = labhelper::loadModelFromOBJ("../scenes/cube.obj");
	generalModel = labhelper::loadModelFromOBJ("../scenes/Lowpoly_tree_sample.obj");

	float scale_factor = 10.f;
	// testModelMatrix = mat4(1.0f);
	// testModelMatrix = translate(scale_factor * worldUp);
	testModelMatrix = scale(1.f * vec3(scale_factor, scale_factor, scale_factor));

	sphereModelMatrix = translate(scale_factor * vec3(1, 1, 0));
	// sphereModelMatrix *= scale(scale_factor * vec3(8.0 / Nx, 8.0 / Nx, 8.0 / Nx));

	cubeModelMatrix = translate(vec3(-0.5, -0.5, -0.5));
	cubeModelMatrix *= translate((scale_factor)*worldUp);

	cubeModelMatrix *= scale(scale_factor * vec3(12.0 / Nx, 12.0 / Nx, 12.0 / Nx));

	occupied_voxels_tree = generate_vexel(generalModel->m_positions, tree_max_length);

	empty_voxels = std::array<bool, SIZE>();
	std::fill(empty_voxels.begin(), empty_voxels.end(), false);

	// occupied_voxels_sphere = generate_vexel(sphereModel->m_positions, sphere_max_length);

	// occupied_voxels_cube = generate_vexel(cubeModel->m_positions, cube_max_length);

	float offset = (scale_factor * 2.0f / Nx);
	//@todo a bug, need to fix, walk around using offset
	generalModelMatrix = translate((scale_factor / 2.0f + offset) * worldUp);

	generalModelMatrix *= scale((scale_factor / tree_max_length) * vec3(1.f));

	// scale the plane
	landingPadModelMatrix = scale(vec3(50.0f, 50.0f, 50.0f));
	// landingPadModelMatrix = mat4(1.0f);

	mmRender = new SmokeRenderer();
	// mmRender->set_occupied_texture(occupied_voxels_cube);
	// mmRender->set_occupied_texture(occupied_voxels_sphere);
	// mmRender->set_occupied_texture(occupied_voxels_tree);

	///////////////////////////////////////////////////////////////////////
	// Load environment map
	///////////////////////////////////////////////////////////////////////
	environmentMap = labhelper::loadHdrTexture("../scenes/envmaps/" + envmap_base_name + ".hdr");

	irradianceMap = labhelper::loadHdrTexture("../scenes/envmaps/" + envmap_base_name + "_irradiance.hdr");

	// Reflection map
	std::vector<std::string> files;
	const int roughnesses = 8;
	for (int i = 0; i < roughnesses; i++)
	{
		files.push_back("../scenes/envmaps/" + envmap_base_name + "_dl_" + std::to_string(i) + ".hdr");
	}

	reflectionMap = labhelper::loadHdrMipmapTexture(files);

	glEnable(GL_DEPTH_TEST); // enable Z-buffering
	glEnable(GL_CULL_FACE);	 // enables backface culling

	// init FBO
	int w, h;
	SDL_GetWindowSize(g_window, &w, &h);
	const int numFBOs = 1;
	for (int i = 0; i < numFBOs; i++)
	{
		FBOList.push_back(FboInfo());
		FBOList[i].resize(w, h);
	}
}

void debugDrawVexel(const std::array<bool, SIZE> &occupied_voxels, const glm::mat4 &viewMatrix,
					const glm::mat4 &projectionMatrix)
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
					labhelper::render(cubeModel);
				}
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
/// This function is used to draw the main objects on the scene
///////////////////////////////////////////////////////////////////////////////
void drawScene(GLuint currentShaderProgram,
			   const mat4 &viewMatrix,
			   const mat4 &projectionMatrix,
			   const mat4 &lightViewMatrix,
			   const mat4 &lightProjectionMatrix)
{
	glUseProgram(currentShaderProgram);
	// Light source
	vec4 viewSpaceLightPosition = viewMatrix * vec4(lightPosition, 1.0f);
	labhelper::setUniformSlow(currentShaderProgram, "point_light_color", point_light_color);
	labhelper::setUniformSlow(currentShaderProgram, "point_light_intensity_multiplier",
							  point_light_intensity_multiplier);
	labhelper::setUniformSlow(currentShaderProgram, "viewSpaceLightPosition", vec3(viewSpaceLightPosition));
	labhelper::setUniformSlow(currentShaderProgram, "viewSpaceLightDir",
							  normalize(vec3(viewMatrix * vec4(-lightPosition, 0.0f))));

	// Environment
	labhelper::setUniformSlow(currentShaderProgram, "environment_multiplier", environment_multiplier);

	// camera
	labhelper::setUniformSlow(currentShaderProgram, "viewInverse", inverse(viewMatrix));

	// landing pad
	{
		labhelper::perf::Scope s("Landing pad");
		labhelper::setUniformSlow(currentShaderProgram, "modelViewProjectionMatrix",
								  projectionMatrix * viewMatrix * landingPadModelMatrix);
		labhelper::setUniformSlow(currentShaderProgram, "modelViewMatrix", viewMatrix * landingPadModelMatrix);
		labhelper::setUniformSlow(currentShaderProgram, "normalMatrix",
								  inverse(transpose(viewMatrix * landingPadModelMatrix)));
		labhelper::render(landingpadModel);
	}

	// draw an object
	// {
	// 	labhelper::setUniformSlow(currentShaderProgram, "modelViewProjectionMatrix",
	// 							  projectionMatrix * viewMatrix * sphereModelMatrix);
	// 	labhelper::setUniformSlow(currentShaderProgram, "modelViewMatrix", viewMatrix * sphereModelMatrix);
	// 	labhelper::setUniformSlow(currentShaderProgram, "normalMatrix",
	// 							  inverse(transpose(viewMatrix * sphereModelMatrix)));
	// 	labhelper::render(sphereModel);
	// }
	// {
	// 	labhelper::setUniformSlow(currentShaderProgram, "modelViewProjectionMatrix",
	// 							  projectionMatrix * viewMatrix * cubeModelMatrix);
	// 	labhelper::setUniformSlow(currentShaderProgram, "modelViewMatrix", viewMatrix * cubeModelMatrix);
	// 	labhelper::setUniformSlow(currentShaderProgram, "normalMatrix",
	// 							  inverse(transpose(viewMatrix * cubeModelMatrix)));
	// 	labhelper::render(cubeModel);
	// }

	// draw a tree
	// {
	// 	labhelper::setUniformSlow(currentShaderProgram, "modelViewProjectionMatrix",
	// 							  projectionMatrix * viewMatrix * generalModelMatrix);
	// 	labhelper::setUniformSlow(currentShaderProgram, "modelViewMatrix", viewMatrix * generalModelMatrix);
	// 	labhelper::setUniformSlow(currentShaderProgram, "normalMatrix",
	// 							  inverse(transpose(viewMatrix * generalModelMatrix)));
	// 	labhelper::render(generalModel);
	// }

	// debugDrawVexel(occupied_voxels_cube, viewMatrix, projectionMatrix);
	// debugDrawVexel(occupied_voxels_sphere, viewMatrix, projectionMatrix);
	// debugDrawVexel(occupied_voxels_tree, viewMatrix, projectionMatrix);

	switch (case_id)
	{
	case 0:
		break;
	case 1:
		debugDrawVexel(occupied_voxels_sphere, viewMatrix, projectionMatrix);
		break;
	case 2:
		debugDrawVexel(occupied_voxels_cube, viewMatrix, projectionMatrix);
		break;
	case 3:
		debugDrawVexel(occupied_voxels_tree, viewMatrix, projectionMatrix);
		break;

	default:
		break;
	}

	{
		static int current_case = -1;
		if (current_case != case_id)
		{
			current_case = case_id;
			switch (case_id)
			{
			case 0:
				mmRender->set_occupied_texture(empty_voxels);
				break;
			case 1:
				mmRender->set_occupied_texture(occupied_voxels_sphere);
				break;
			case 2:
				mmRender->set_occupied_texture(occupied_voxels_cube);
				break;
			case 3:
				mmRender->set_occupied_texture(occupied_voxels_tree);
				break;
			default:
				break;
			}
		}
	}

	{
		labhelper::perf::Scope s("render smoke");

		glUseProgram(smokeProgram);
		labhelper::setUniformSlow(smokeProgram, "modelViewProjectionMatrix",
								  projectionMatrix * viewMatrix * testModelMatrix);
		labhelper::setUniformSlow(smokeProgram, "modelMatrix", testModelMatrix);
		labhelper::setUniformSlow(smokeProgram, "worldSpaceLightPosition", lightPosition);
		labhelper::setUniformSlow(smokeProgram, "pointLightIntensity", point_light_intensity_multiplier);
		labhelper::setUniformSlow(smokeProgram, "worldSpaceCameraPosition", cameraPosition);
		labhelper::setUniformSlow(smokeProgram, "factor", smoke_factor);
		labhelper::setUniformSlow(smokeProgram, "enable_light_tracing", enable_light_tracing);

		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		// mmRender->render(generateSphereDensity());
		// mmRender->render(generateCubeDensity());
		// try to lock the simulator
		{
			std::lock_guard<std::mutex> lock(simLock);
			mmRender->render(density, transparency);
		}

		glDisable(GL_BLEND);

		// glDisable(GL_DEPTH_TEST);
		// mmRender->render_frame(projectionMatrix * viewMatrix * testModelMatrix);
		// glEnable(GL_DEPTH_TEST);
	}
	// draw a line
	glDisable(GL_DEPTH_TEST);
	//
	drawLine(vec3(0.0f), vec3(10.0, 0.0, 0.0), projectionMatrix * viewMatrix * mat4(1.0), vec3(1.0f, 0.0f, 0.0f));
	drawLine(vec3(0.0f), vec3(0.0, 10.0, 0.0), projectionMatrix * viewMatrix * mat4(1.0), vec3(0.0f, 1.0f, 0.0f));
	drawLine(vec3(0.0f), vec3(0.0, 0.0, 10.0), projectionMatrix * viewMatrix * mat4(1.0), vec3(0.0f, 0.0f, 1.0f));
	//
	glEnable(GL_DEPTH_TEST);
}

///////////////////////////////////////////////////////////////////////////////
/// This function will be called once per frame, so the code to set up
/// the scene for rendering should go here
///////////////////////////////////////////////////////////////////////////////
void display(void)
{
	labhelper::perf::Scope s("Display");

	///////////////////////////////////////////////////////////////////////////
	// Check if window size has changed and resize buffers as needed
	///////////////////////////////////////////////////////////////////////////
	{
		int w, h;
		SDL_GetWindowSize(g_window, &w, &h);
		if (w != windowWidth || h != windowHeight)
		{
			windowWidth = w;
			windowHeight = h;
			for (auto &fbo : FBOList)
			{
				fbo.resize(w, h);
			}
		}
	}

	///////////////////////////////////////////////////////////////////////////
	// setup matrices
	///////////////////////////////////////////////////////////////////////////
	mat4 projMatrix = perspective(radians(40.0f), float(windowWidth) / float(windowHeight), 10.f, 1000.0f);
	mat4 viewMatrix = lookAt(cameraPosition, cameraPosition + cameraDirection, worldUp);

	mat4 lightViewMatrix = lookAt(lightPosition, vec3(0.0f), worldUp);
	static mat4 lightProjMatrix = perspective(radians(45.0f), 1.0f, 25.0f, 100.0f);

	///////////////////////////////////////////////////////////////////////////
	// Bind the environment map(s) to unused texture units
	///////////////////////////////////////////////////////////////////////////
	glActiveTexture(GL_TEXTURE6);
	glBindTexture(GL_TEXTURE_2D, environmentMap);
	glActiveTexture(GL_TEXTURE7);
	glBindTexture(GL_TEXTURE_2D, irradianceMap);
	glActiveTexture(GL_TEXTURE8);
	glBindTexture(GL_TEXTURE_2D, reflectionMap);
	glActiveTexture(GL_TEXTURE0);

	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	glViewport(0, 0, windowWidth, windowHeight);
	glClearColor(0.2f, 0.2f, 0.8f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	{
		labhelper::perf::Scope s("Scene");
		drawScene(shaderProgram, viewMatrix, projMatrix, lightViewMatrix, lightProjMatrix);
	}
}

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
			mat4 yaw = rotate(rotationSpeed * deltaTime * -delta_x, worldUp);
			mat4 pitch = rotate(rotationSpeed * deltaTime * -delta_y,
								normalize(cross(cameraDirection, worldUp)));
			cameraDirection = vec3(pitch * yaw * vec4(cameraDirection, 0.0f));
			g_prevMouseCoords.x = event.motion.x;
			g_prevMouseCoords.y = event.motion.y;
		}
	}

	// check keyboard state (which keys are still pressed)
	const uint8_t *state = SDL_GetKeyboardState(nullptr);
	vec3 cameraRight = cross(cameraDirection, worldUp);

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

///////////////////////////////////////////////////////////////////////////////
/// This function is to hold the general GUI logic
///////////////////////////////////////////////////////////////////////////////
void ControlPanel()
{
	ImGui::Begin("Control panel");
	ImGui::SliderFloat("Light intensity", &point_light_intensity_multiplier, 0.0f, 100000.0f);
	ImGui::SliderFloat("Smoke factor", &smoke_factor, 0.0f, 100.0f);
	ImGui::SliderFloat("Env temperature", &g_env_temp, 0.0f, 1000.0f);

	if (ImGui::Button("Case 0: Empty"))
	{
		case_id = 0;
	}

	if (ImGui::Button("Case 1: Sphere"))
	{
		case_id = 1;
	}

	if (ImGui::Button("Case 2: Cube"))
	{
		case_id = 2;
	}

	if (ImGui::Button("Case 3: Tree"))
	{
		case_id = 3;
	}

	// set camera position
	if (ImGui::Button("View Front"))
	{
		cameraPosition = vec3(30.f, 5.0f, 5.f);
		cameraDirection = normalize(-vec3(1.f, 0.f, 0.f));
	}
	// set camera position
	if (ImGui::Button("View Top"))
	{
		cameraPosition = vec3(5.f, 30.f, 5.f);
		cameraDirection = normalize(-vec3(0.0f, 1.f, 0.1f));
	}
	// set camera position
	if (ImGui::Button("View Side"))
	{
		cameraPosition = vec3(5.f, 5.0f, 30.f);
		cameraDirection = normalize(-vec3(0.0f, 0.f, 1.f));
	}

	// reset button, to reset simulator
	if (ImGui::Button("Reset"))
	{
		simulator_rest_trigger = true;
	}

	ImGui::End();
}

void gui()
{
	// ----------------- Set variables --------------------------
	ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate,
				ImGui::GetIO().Framerate);

	{
		std::lock_guard<std::mutex> lock(simLock);
		ImGui::Text("Simulator Info: %s", simulator_info.c_str());
	}

	// ----------------------------------------------------------
	labhelper::perf::drawEventsWindow();

	ControlPanel();
}

int main(int argc, char *argv[])
{
	// DEBUG_PRINT("MAIN cuda");
	simulator = std::make_unique<Simulator>(ttime);
	cudaRender = std::make_unique<MCUDA::CudaRender>(SIZE, Nx, Ny, Nz);
	cudaRender->init();

	g_window = labhelper::init_window_SDL("OpenGL Project");

	initialize();

	bool stopRendering = false;
	auto startTime = std::chrono::system_clock::now();

	// thred to run simulation
	std::thread simThread([&]()
						  {
		DEBUG_PRINT("Simulation thread started");
		while (!stopRendering)
		{	
			if (simulator_rest_trigger)
			{	
				simulator->setEnvTemperature(g_env_temp);
				simulator->reset();
				simulator_rest_trigger = false;
			}

			static int current_case = -1;
			if(current_case != case_id){
				current_case = case_id;
				switch (case_id)
				{
				case 0:
					simulator->setOccupiedVoxels(empty_voxels);
					break;
				case 1:
					simulator->setOccupiedVoxels(occupied_voxels_sphere);
					break;
				case 2:
					simulator->setOccupiedVoxels(occupied_voxels_cube);
					break;
				case 3:
					simulator->setOccupiedVoxels(occupied_voxels_tree);
					break;
				default:
					break;
				}
				simulator->reset();
			}

			simulator->update();
			// update the simulator
			if (ttime > FINISH_TIME){
				// DEBUG_PRINT("Simulation  finished");
				// sleep for a while
				std::this_thread::sleep_for(std::chrono::milliseconds(100));
				continue;
			}

			{
				std::lock_guard<std::mutex> lock(simLock);
				// copy the density
				std::copy(
					simulator->getDensity(), 
					simulator->getDensity()+SIZE, 
					density.begin()
				);
				// gen transparency
				T_START("genTransparencyMap");
				cudaRender->genTransparencyMap(
					density.data(),
					transparency,
					lightPosition.x, lightPosition.y, lightPosition.z,
					10.0, smoke_factor
				);
				T_END;

				simulator_info = simulator->getPerformanceInfo();
			}
		} });

	while (!stopRendering)
	{
		// update currentTime
		std::chrono::duration<float> timeSinceStart = std::chrono::system_clock::now() - startTime;
		previousTime = currentTime;
		currentTime = timeSinceStart.count();
		deltaTime = currentTime - previousTime;

		// check events (keyboard among other)
		stopRendering = handleEvents();

		// Inform imgui of new frame
		labhelper::newFrame(g_window);

		// Render overlay GUI.
		gui();

		labhelper::perf::Scope scope("Main loop");
		// render to window
		display();

		// Finish the frame and render the GUI
		{
			labhelper::perf::Scope s("finish frame");
			labhelper::finishFrame();
		}
		// Swap front and back buffer. This frame will now been displayed.
		{
			labhelper::perf::Scope s("SDL_GL_SwapWindow");
			SDL_GL_SwapWindow(g_window);
		}
		// limit the frame rate: 30
		std::chrono::duration<float> endTime = std::chrono::system_clock::now() - startTime;
		auto frame_time = endTime.count() - previousTime;
		if (frame_time < 1.0f / 30.0f)
		{
			std::this_thread::sleep_for(
				std::chrono::milliseconds(
					(int)(1000.0f / 30.0f - frame_time * 1000.0f)));
		}
		// DEBUG_PRINT("Frame time:" << frame_time);
	}
	// Free Models
	labhelper::freeModel(landingpadModel);
	labhelper::freeModel(pointLight);
	// Shut down everything. This includes the window and all other subsystems.
	labhelper::shutDown(g_window);
	return 0;
}
