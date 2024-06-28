#ifndef __MOCK_H__
#define __MOCK_H__

#include <array>
#include "constants.h"
#include "common/mmath.h"

namespace MOCK{

std::array<double, SIZE> &
generateSphereDensity();

std::array<double, SIZE> &
generateCubeDensity();

} // namespace MOCK

#endif // __MOCK_H__