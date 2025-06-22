#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include <igl/PI.h>
#include <geometrycentral/utilities/vector3.h>
#include <openvdb/openvdb.h>

using namespace geometrycentral;

namespace modules {
    namespace scene_file {
		enum VolumeType {
			PRIMITIVE, // Primitive surface (e.g., sphere, box)
			SDF, // Signed distance field surface (currently not supported)
			MESH // Mesh surface (either to be treated as SDF or as nodes that impact the Medial Axis Energy)
		};

        enum PrimitiveType {
            SPHERE,
            BOX,
            ROUNDBOX,
            TORUS
        };

        struct SceneObject_Volume {
			VolumeType volumeType;

            // Primitive params
            PrimitiveType primitive_type = PrimitiveType::SPHERE;
            std::vector<double> primitive_params = {}; // Parameters for the primitive (e.g., radius for sphere, extents for box)

            // Mesh params
			std::string mesh_filename = "";
            bool convert_to_sdf;
            double mesh_to_sdf_voxelsize = 0.01f;
            float mesh_to_sdf_halfwidth = 3.0f;
			std::vector<Vector3> mesh_points = {}; // initialized through code
            openvdb::FloatGrid::Ptr sdf; // initialized through code
        };

        // We may not need all of these props,
        // but they all have a default value assigned,
        // so they don't bother at the moment.
        // -> ToDo: Remove unnecessary props
        struct SceneObject {
            std::string curveFileName = "";
            double radius = 0.1;
            double rmax = radius * 5;
            double h = igl::PI * radius / 20;
            double p = 2;
            double q = 2;

			SceneObject_Volume volume;
            bool visualizeVolume = false;

            double timestep = 1.;
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
