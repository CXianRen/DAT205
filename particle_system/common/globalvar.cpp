#include "globalvar.h"

GLuint shaderProgram;       // Shader for rendering the final image
GLuint simpleShaderProgram; // Shader used to draw the shadow map
GLuint debugLineProgram; // Shader used to draw debug lines
GLuint smokeProgram;

///////////////////////////////////////////////////////////////////////////////
// Camera parameters.
///////////////////////////////////////////////////////////////////////////////
glm::vec3 cameraPosition(45.0f, 45.0f, 45.0f);
glm::vec3 cameraDirection = glm::normalize(glm::vec3(0.0f) - cameraPosition);
float cameraSpeed = 10.f;
const glm::vec3 worldUp(0.0f, 1.0f, 0.0f);

// light
float g_point_light_intensity = 15000;
glm::vec3 g_light_position = glm::vec4(80.0f, 25.0f, 25.0f, 1.0f);
glm::vec3 g_point_light_color = glm::vec3(1.f, 1.f, 1.f);

// time setting 
float currentTime = 0.0f;
float previousTime = 0.0f;
float deltaTime = 0.0f;

// case id
int g_case_id = 0;

// simulator parameters
#include "particle/constants.h"
std::string simulator_info;

float g_env_temp = 25.0;
float g_alpha = ALPHA;
float g_beta = BETA;
float g_vort_eps = VORT_EPS;
float g_decay_factor = 0.99;
float g_dt = DT;

std::mutex g_sim_lock;
bool g_simulator_rest = false;

// smoke render parameters
float g_smoke_factor = 10.f;