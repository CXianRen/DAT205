#ifndef __M_GLOBALVAR_H__
#define __M_GLOBALVAR_H__

#include <GL/glew.h>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <chrono>

#include <labhelper.h>
#include <imgui.h>

#include <perf.h>

#include <glm/glm.hpp>
#include <glm/gtx/transform.hpp>

#include <Model.h>
#include "hdr.h"
#include "fbo.h"

///////////////////////////////////////////////////////////////////////////////
// Shader programs
///////////////////////////////////////////////////////////////////////////////
extern GLuint shaderProgram;       // Shader for rendering the final image
extern GLuint simpleShaderProgram; // Shader used to draw the shadow map
extern GLuint debugLineProgram; // Shader used to draw debug lines
extern GLuint smokeProgram;     // Shader used to draw smoke

///////////////////////////////////////////////////////////////////////////////
// Camera parameters.
///////////////////////////////////////////////////////////////////////////////
extern glm::vec3 cameraPosition;
extern glm::vec3 cameraDirection;
extern float cameraSpeed;
extern const glm::vec3 worldUp;

// light
extern float g_point_light_intensity;
extern glm::vec3 g_light_position;
extern glm::vec3 g_point_light_color;

// time setting
extern float currentTime;
extern float previousTime;
extern float deltaTime;

// case id
extern int g_case_id;

// simulator parameters
extern std::string simulator_info;

extern float g_env_temp;
extern float g_alpha;
extern float g_beta;
extern float g_vort_eps;
extern float g_decay_factor;
extern float g_dt;

#include <mutex>
extern std::mutex g_sim_lock;
extern bool g_simulator_rest;

// smoke render parameters
extern float g_smoke_factor;

// demo sence
enum Demo
{
    SMOKE_EMPTY,
    SMOKE_CUBE,
    SMOKE_SPHERE,
    SMOKE_TREE,
    UNDEFINED
};

#endif // __M_GLOBALVAR_H__