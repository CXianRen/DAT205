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
using namespace glm;

#include <Model.h>
#include "hdr.h"
#include "fbo.h"

#include "../particle/simulator.h"

///////////////////////////////////////////////////////////////////////////////
// Shader programs
///////////////////////////////////////////////////////////////////////////////
extern GLuint shaderProgram;	   // Shader for rendering the final image
extern GLuint simpleShaderProgram; // Shader used to draw the shadow map
extern GLuint backgroundProgram;
extern GLuint debugLineProgram;    // Shader used to draw debug lines

#endif // __M_GLOBALVAR_H__