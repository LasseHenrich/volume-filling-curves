// volume-filling-curves.h : Include file for standard system include files,
// or project specific include files.

#pragma once

#include <iostream>
#include <vector>
#include <array>
#include <geometrycentral/utilities/vector3.h>

using namespace geometrycentral;

// TODO: Reference additional headers your program requires here.

struct Curve {
    std::vector<Vector3> nodes;
    std::vector<std::array<int, 2>> segments;
    std::vector<double> segmentLengths;
};