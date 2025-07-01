#pragma once

#include <iostream>
#include <vector>
#include <array>
#include <geometrycentral/utilities/vector3.h>
#include <geometrycentral/surface/manifold_surface_mesh.h>

using namespace geometrycentral;

namespace modules {
    struct Curve {
        std::vector<Vector3> nodes;
        std::vector<std::array<int, 2>> segments;
        std::vector<double> segmentLengths;
    };

    struct Surface {
        std::unique_ptr<geometrycentral::surface::ManifoldSurfaceMesh> mesh;
        geometrycentral::surface::VertexData<Vector3> vertexPositions;
    };
}
