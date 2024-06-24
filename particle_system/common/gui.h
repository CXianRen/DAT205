#ifndef __MGUI_H__
#define __MGUI_H__
#include "labhelper.h"

#include "common/globalvar.h"
#include "particle/constants.h"

void ControlPanel();

void gui();

void drawLine(const glm::vec3 start_pos, const glm::vec3 end_pos,
              const glm::mat4 &projectionViewModelMatrix,
              const glm::vec3 color = glm::vec3(1.0f, 0.0f, 0.0f));

void drawCoordinate(const glm::mat4 &projectionViewModelMatrix);


void drawVexel(const std::array<bool, SIZE> &occupied_voxels, const glm::mat4 &viewMatrix,
               const glm::mat4 &projectionMatrix, const labhelper::Model *model);

bool handleEvents(void);
#endif // __MGUI_H__