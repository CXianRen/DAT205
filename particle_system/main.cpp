#include <thread>
#include <mutex>

#include "common/debug.h"
#include "common/gui.h"
#include "common/globalvar.h"
#include "common/mperf.h"

#include "example.h"
#include "particle/render.h"
#include "particle/Simulator.h"
#include "particle/CudaRender.h"

using namespace glm;

///////////////////////////////////////////////////////////////////////////////
// Various globals
///////////////////////////////////////////////////////////////////////////////
SDL_Window *g_window = nullptr;
int windowWidth, windowHeight;
std::vector<FboInfo> FBOList;

///////////////////////////////////////////////////////////////////////////////
// Models
///////////////////////////////////////////////////////////////////////////////
labhelper::Model *landingpadModel = nullptr;
labhelper::Model *cubeModel = nullptr;

glm::mat4 landingPadModelMatrix;
glm::mat4 SmokeZoneMatrix;

///////////////////////////////////////////////////////////////////////////////
// Particle system
///////////////////////////////////////////////////////////////////////////////
// Simulator simulator;
double ttime = 0.0;
std::array<double, SIZE> density;
double transparency[SIZE];

std::unique_ptr<Simulator> simulator;
std::unique_ptr<MCUDA::CudaRender> cudaRender;
std::unique_ptr<SmokeRenderer> mmRender;

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
	cubeModel = labhelper::loadModelFromOBJ("../scenes/cube.obj");

	float scale_factor = 10.f;
	SmokeZoneMatrix = scale(1.f * vec3(scale_factor, scale_factor, scale_factor));

	// scale the plane
	landingPadModelMatrix = scale(vec3(50.0f, 50.0f, 50.0f));

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

	// init simulator
	mmRender = std::make_unique<SmokeRenderer>();
	simulator = std::make_unique<Simulator>(ttime);
	cudaRender = std::make_unique<MCUDA::CudaRender>(SIZE, Nx, Ny, Nz);
	cudaRender->init();

	init_examples();
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
	vec4 viewSpaceLightPosition = viewMatrix * vec4(g_light_position, 1.0f);
	labhelper::setUniformSlow(currentShaderProgram, "g_point_light_color", g_point_light_color);
	labhelper::setUniformSlow(currentShaderProgram, "pointLightIntensity",
							  g_point_light_intensity);
	labhelper::setUniformSlow(currentShaderProgram, "viewSpaceLightPosition", vec3(viewSpaceLightPosition));
	labhelper::setUniformSlow(currentShaderProgram, "viewSpaceLightDir",
							  normalize(vec3(viewMatrix * vec4(-g_light_position, 0.0f))));

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

	if (examples.find(g_case_id) != examples.end())
	{
		auto &exp = examples[g_case_id];
		drawVexel(exp->obj_vexels, viewMatrix, projectionMatrix, cubeModel);
	}
	else
	{
		DEBUG_PRINT("Undefined case id");
	}

	{
		labhelper::perf::Scope s("render smoke");

		static Demo current_case = Demo::UNDEFINED;
		if (current_case != g_case_id)
		{
			current_case = g_case_id;
			if (examples.find(g_case_id) != examples.end())
			{
				auto &exp = examples[g_case_id];
				mmRender->set_occupied_texture(exp->obj_vexels);
			}
			else
			{
				DEBUG_PRINT("Undefined case id");
			}
		}

		glUseProgram(smokeProgram);
		labhelper::setUniformSlow(smokeProgram, "modelViewProjectionMatrix",
								  projectionMatrix * viewMatrix * SmokeZoneMatrix);
		labhelper::setUniformSlow(smokeProgram, "modelMatrix", SmokeZoneMatrix);
		labhelper::setUniformSlow(smokeProgram, "worldSpaceLightPosition", g_light_position);
		labhelper::setUniformSlow(smokeProgram, "pointLightIntensity", g_point_light_intensity);
		labhelper::setUniformSlow(smokeProgram, "worldSpaceCameraPosition", cameraPosition);
		labhelper::setUniformSlow(smokeProgram, "factor", g_smoke_factor);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		// try to lock the simulator
		{
			std::lock_guard<std::mutex> lock(g_sim_lock);
			mmRender->render(density, transparency);
		}
		glDisable(GL_BLEND);
	}

	drawCoordinate(projectionMatrix * viewMatrix);
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

	mat4 lightViewMatrix = lookAt(g_light_position, vec3(0.0f), worldUp);
	static mat4 lightProjMatrix = perspective(radians(45.0f), 1.0f, 25.0f, 100.0f);

	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	glViewport(0, 0, windowWidth, windowHeight);
	glClearColor(0.2f, 0.2f, 0.8f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	{
		labhelper::perf::Scope s("Scene");
		drawScene(shaderProgram, viewMatrix, projMatrix, lightViewMatrix, lightProjMatrix);
	}
}

int main(int argc, char *argv[])
{
	// DEBUG_PRINT("MAIN cuda");

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
			simulator->setAlpha(g_alpha);
			simulator->setBeta(g_beta);
			simulator->setVortEps(g_vort_eps);
			simulator->setDecayFactor(g_decay_factor);
			simulator->setDt(g_dt);
			
			if (g_simulator_rest)
			{	
				simulator->setEnvTemperature(g_env_temp);
				simulator->reset();
				g_simulator_rest = false;
			}

			static int current_case = -1;
			if(current_case != g_case_id)
			{
				current_case = g_case_id;
				
				if (examples.find(g_case_id) != examples.end())
				{
					auto &exp = examples[g_case_id];
					simulator->setOccupiedVoxels(exp->obj_vexels);
				}
				else
				{
					DEBUG_PRINT("Undefined case id");
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
				std::lock_guard<std::mutex> lock(g_sim_lock);
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
					g_light_position.x, g_light_position.y, g_light_position.z,
					10.0, g_smoke_factor
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
	}
	// Free Models
	labhelper::freeModel(landingpadModel);
	labhelper::freeModel(cubeModel);

	// Shut down everything. This includes the window and all other subsystems.
	labhelper::shutDown(g_window);
	return 0;
}
