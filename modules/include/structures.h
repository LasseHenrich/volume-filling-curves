#pragma once

#include <iostream>
#include <vector>
#include <array>
#include <geometrycentral/utilities/vector3.h>
#include <geometrycentral/surface/manifold_surface_mesh.h>
#include <geometrycentral/surface/vertex_position_geometry.h>

using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace modules {
    struct Curve {
        std::vector<Vector3> nodes;
        std::vector<std::array<int, 2>> segments;
        std::vector<double> segmentLengths;
    };

    struct Surface {
        std::unique_ptr<ManifoldSurfaceMesh> mesh;
        std::unique_ptr<VertexPositionGeometry> geometry;
    };
}
