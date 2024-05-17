
#ifdef _WIN32
extern "C" _declspec(dllexport) unsigned int NvOptimusEnablement = 0x00000001;
#endif

#include "common/debug.h"
#include "common/globalvar.h"
#include "particle/render.h"
#include "particle/Simulator.h"
#include "particle/MACGrid.h"


std::shared_ptr<MACGrid> grids = std::make_shared<MACGrid>();
double ttime = 0.0;
std::unique_ptr<Simulator> simulator = std::make_unique<Simulator>(grids, ttime);

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
vec3 lightPosition;
vec3 point_light_color = vec3(1.f, 1.f, 1.f);

float point_light_intensity_multiplier = 100000.0f;
bool step_light = true;
float light_rotation_step = 0.01f;

///////////////////////////////////////////////////////////////////////////////
// Models
///////////////////////////////////////////////////////////////////////////////

labhelper::Model *landingpadModel = nullptr;
labhelper::Model *pointLight = nullptr;
labhelper::Model *particleModel = nullptr;
labhelper::Model *cubeModel = nullptr;
SmokeRenderer *mmRender = nullptr;
mat4 landingPadModelMatrix;
mat4 testModelMatrix;

///////////////////////////////////////////////////////////////////////////////
// Particle system
///////////////////////////////////////////////////////////////////////////////
// Simulator simulator;

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
	particleModel = labhelper::loadModelFromOBJ("../scenes/particle.obj");
	cubeModel = labhelper::loadModelFromOBJ("../scenes/cube.obj");

	testModelMatrix = translate(15.0f * worldUp);
	testModelMatrix *= scale(vec3(10.0f, 10.0f, 10.0f));
	// scale the plane
	landingPadModelMatrix = scale(vec3(50.0f, 50.0f, 50.0f));
	// landingPadModelMatrix = mat4(1.0f);

	mmRender = new SmokeRenderer();

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
	const int numFBOs = 5;
	for (int i = 0; i < numFBOs; i++)
	{
		FBOList.push_back(FboInfo());
		FBOList[i].resize(w, h);
	}
}

void debugDrawLight(const glm::mat4 &viewMatrix,
					const glm::mat4 &projectionMatrix,
					const glm::vec3 &worldSpaceLightPos)
{
	mat4 modelMatrix = glm::translate(worldSpaceLightPos);
	glUseProgram(shaderProgram);
	labhelper::setUniformSlow(shaderProgram, "modelViewProjectionMatrix",
							  projectionMatrix * viewMatrix * modelMatrix);
	// glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	labhelper::render(pointLight);
	// glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

void drawBackground(const mat4 &viewMatrix, const mat4 &projectionMatrix)
{
	glUseProgram(backgroundProgram);
	labhelper::setUniformSlow(backgroundProgram, "environment_multiplier", environment_multiplier);
	labhelper::setUniformSlow(backgroundProgram, "inv_PV", inverse(projectionMatrix * viewMatrix));
	labhelper::setUniformSlow(backgroundProgram, "camera_pos", cameraPosition);
	labhelper::drawFullScreenQuad();
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
	labhelper::setUniformSlow(currentShaderProgram, "modelViewProjectionMatrix",
							  projectionMatrix * viewMatrix * landingPadModelMatrix);
	labhelper::setUniformSlow(currentShaderProgram, "modelViewMatrix", viewMatrix * landingPadModelMatrix);
	labhelper::setUniformSlow(currentShaderProgram, "normalMatrix",
							  inverse(transpose(viewMatrix * landingPadModelMatrix)));

	labhelper::render(landingpadModel);

	// draw a particle
	labhelper::setUniformSlow(currentShaderProgram, "modelViewProjectionMatrix",
							  projectionMatrix * viewMatrix * testModelMatrix);
	labhelper::setUniformSlow(currentShaderProgram, "modelViewMatrix", viewMatrix * testModelMatrix);
	labhelper::setUniformSlow(currentShaderProgram, "normalMatrix",
							  inverse(transpose(viewMatrix * testModelMatrix)));

	// glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	// glDisable(GL_DEPTH_TEST);
	// glDisable(GL_CULL_FACE);
	// labhelper::render(cubeModel);
	// glEnable(GL_CULL_FACE);
	// glEnable(GL_DEPTH_TEST);
	// glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	glUseProgram(smokeProgram);
	labhelper::setUniformSlow(smokeProgram, "modelViewProjectionMatrix",
							  projectionMatrix * viewMatrix * testModelMatrix);
	labhelper::setUniformSlow(smokeProgram, "modelMatrix", testModelMatrix);
	labhelper::setUniformSlow(smokeProgram, "worldSpaceLightPosition", lightPosition);
	labhelper::setUniformSlow(smokeProgram, "pointLightIntensity", point_light_intensity_multiplier);
	labhelper::setUniformSlow(smokeProgram, "worldSpaceCameraPosition", cameraPosition);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	// mmRender->render(generateSphereDensity());
	mmRender->render(grids->density.m_data);
	glDisable(GL_BLEND);

	glDisable(GL_DEPTH_TEST);
	mmRender->render_frame(projectionMatrix * viewMatrix * testModelMatrix);
	glEnable(GL_DEPTH_TEST);

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
	mat4 projMatrix = perspective(radians(45.0f), float(windowWidth) / float(windowHeight), 5.0f, 2000.0f);
	mat4 viewMatrix = lookAt(cameraPosition, cameraPosition + cameraDirection, worldUp);

	static vec3 wheatley_position(-35.0f, 35.0f, -35.0f);

	static vec4 lightStartPosition = vec4(80.0f, 80.0f, 80.0f, 1.0f);
	if (!step_light)
	{

		light_rotation_step += 0.01f;
	}
	lightPosition = vec3(rotate(light_rotation_step, worldUp) * lightStartPosition);

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

	///////////////////////////////////////////////////////////////////////////
	// Draw scene from security camera
	///////////////////////////////////////////////////////////////////////////
	// FboInfo &securityCamFBO = FBOList[0];
	// glBindFramebuffer(GL_FRAMEBUFFER, securityCamFBO.framebufferId);
	// glViewport(0, 0, securityCamFBO.width, securityCamFBO.height);
	// glClearColor(0.2f, 0.2f, 0.8f, 1.0f);
	// glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// drawBackground(viewMatrix, projMatrix);

	// drawScene( shaderProgram,
	// 			securityCamViewMatrix, projMatrix, lightViewMatrix, lightProjMatrix );

	// labhelper::Material& tv_screen = landingpadModel->m_materials[8];
	// tv_screen.m_emission_texture.gl_id = securityCamFBO.colorTextureTargets[0];

	///////////////////////////////////////////////////////////////////////////
	// Draw from camera
	///////////////////////////////////////////////////////////////////////////
	// glBindFramebuffer(GL_FRAMEBUFFER, cameraFBO.framebufferId);
	// glViewport(0, 0, cameraFBO.width, cameraFBO.height);

	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	glViewport(0, 0, windowWidth, windowHeight);
	glClearColor(0.2f, 0.2f, 0.8f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	{
		// labhelper::perf::Scope s( "Background" );
		// drawBackground(viewMatrix, projMatrix);
	}
	{
		labhelper::perf::Scope s("Scene");
		drawScene(shaderProgram, viewMatrix, projMatrix, lightViewMatrix, lightProjMatrix);
	}

	// debugDrawLight(viewMatrix, projMatrix, vec3(lightPosition));
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
	ImGui::Selectable("Step light", &step_light, 0, ImVec2(0, 0));
	if (ImGui::Button("Light move one step") && step_light)
	{
		light_rotation_step += 0.01f;
	}
	// select the light color
	ImGui::ColorEdit3("Light color", &point_light_color[0]);

	ImGui::End();
}

void gui()
{
	// ----------------- Set variables --------------------------
	ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate,
				ImGui::GetIO().Framerate);
	// ----------------------------------------------------------
	labhelper::perf::drawEventsWindow();

	ControlPanel();
}

int main(int argc, char *argv[])
{
	g_window = labhelper::init_window_SDL("OpenGL Project");

	initialize();
	// mmRender.init();


	bool stopRendering = false;
	auto startTime = std::chrono::system_clock::now();

	while (!stopRendering)
	{

		//@tood using multithread to update the particle system
		simulator->update();

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

		// render to window
		display();

		// Finish the frame and render the GUI
		labhelper::finishFrame();

		// Swap front and back buffer. This frame will now been displayed.
		SDL_GL_SwapWindow(g_window);
	}
	// Free Models
	labhelper::freeModel(landingpadModel);
	labhelper::freeModel(pointLight);

	// Shut down everything. This includes the window and all other subsystems.
	labhelper::shutDown(g_window);
	return 0;
}
