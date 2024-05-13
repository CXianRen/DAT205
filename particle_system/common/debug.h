#ifndef __M_DEBUG_H__
#define __M_DEBUG_H__

#include <iostream>

#include <glm/glm.hpp>

#define DEBUG_PRINT(msg) \
    std::cout << "[DEBUG] :" msg << std::endl;


void drawLine(const glm::vec3 start_pos, const glm::vec3 end_pos,
			  const glm::mat4 &viewMatrix, const glm::mat4 &projectionMatrix, 
			  const glm::vec3 color = glm::vec3(1.0f, 0.0f, 0.0f));

// void draw3DGird(const glm::mat4 &viewMatrix, 
//                 const glm::mat4 &projectionMatrix,
//                 const glm::vec3 center,
//                 int gX, int gY, int gZ, float gCONST_h
//                 ,const glm::vec3 color = glm::vec3(1.0f, 0.0f, 0.0f));

#endif // __M_DEBUG_H__