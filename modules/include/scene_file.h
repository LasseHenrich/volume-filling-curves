#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include <igl/PI.h>
#include <geometrycentral/utilities/vector3.h>

using namespace geometrycentral;

namespace modules {
    namespace scene_file {
        // We may not need all of these props,
        // but they all have a default value assigned,
        // so they don't bother at the moment.
        // -> ToDo: Remove unnecessary props
        struct SceneObject {
            std::string curveFileName = "";
            double radius = 0.1;
            double timestep = 1.;
            double h = igl::PI * radius / 20;
            double p = 2;
            double q = 2;
            double rmax = radius * 5;
            double w_fieldAlignedness = 0;
            double w_curvatureAlignedness = 0;
            double w_bilaplacian = 0;
            bool varyingAlpha = false;
            bool useGeodesicMedialAxis = false;
            bool excecuteOnly = false;
        };
    }

    void splitString(const std::string& str, std::vector<std::string>& cont, char delim = ' ');

    scene_file::SceneObject read_scene(std::string filename);

	std::tuple<
		std::vector<Vector3>, // nodes
		std::vector<std::array<int, 2>> // segments
    > read_nodes(std::string filename);
}
