#ifndef __M_DEBUG_H__
#define __M_DEBUG_H__

#include <iostream>

#include <glm/glm.hpp>

#define DEBUG_PRINT(msg) \
    std::cout << "[DEBUG] :" msg << std::endl;


void drawLine(const glm::vec3 start_pos, const glm::vec3 end_pos,
			  const glm::mat4 &projectionViewModelMatrix, 
			  const glm::vec3 color = glm::vec3(1.0f, 0.0f, 0.0f));

#endif // __M_DEBUG_H__