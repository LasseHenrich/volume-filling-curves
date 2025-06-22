#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include <polyscope/curve_network.h>
#include <polyscope/point_cloud.h>
#include <polyscope/volume_mesh.h>

#include "args/args.hxx"

#include <openvdb/openvdb.h>
#include <openvdb/tools/Interpolation.h>

#include "volume-filling-curves.h"
#include <chrono>
#include <volume_filling_energy.h>
#include <volume_path_evolution.h>
#include <scene_file.h>
#include <mesh_parser.h>

using namespace modules;

scene_file::SceneObject scene;

std::vector<Vector3> nodes = {}, restNodes = {}, initialNodes = {};
std::vector<std::array<int, 2>> segments = {}, restSegments = {}, initialSegments = {};
std::vector<double> segmentLengths = {}, restSegmentLengths = {}, initialSegmentLengths = {};

int iteration = 0;
bool runLoop = false;

std::tuple<
	std::vector<Vector3>, // nodes
	std::vector<std::array<int, 2>>, // segments
	std::vector<double> // norms of the vectors
> formatVectorsForVisualization(
    std::vector<Vector3>& nodes,
    std::vector<Vector3>& directions
) {
	assert(nodes.size() == directions.size(), "nodes and directions must have the same size");

	std::vector<Vector3> dNodes = nodes;
    std::vector<std::array<int, 2>> dSegments = {};
    std::vector<double> dNorm = {};
    
    for (int i = 0; i < directions.size(); i++) {
        dNodes.push_back(nodes[i] + directions[i]);
		dNorm.push_back(norm(directions[i]));
    }

    for (int i = 0; i < directions.size(); i++) {
		int oi = i + directions.size();
		std::array<int, 2> line = { i, oi };
		dSegments.push_back(line);
	}

	return std::make_tuple(dNodes, dSegments, dNorm);
}

void doWork() {
    iteration++;

    restNodes = nodes;
    restSegments = segments;
    restSegmentLengths = segmentLengths;

    std::cout << "===== iteration: " << iteration << "=====" << std::endl;

    std::cout << "numNodes: " << nodes.size() << std::endl;
    std::cout << "numSegments: " << segments.size() << std::endl;
    std::cout << std::endl;

    auto start = std::chrono::high_resolution_clock::now();

    // 1. compute descent direction
	auto [descent, gradient, energy, medialAxis] = modules::volume_filling_energy(
		nodes,
		segments,
		segmentLengths,
		scene
	);

    // 2. evolve curve without self-intersections
	std::tie(nodes, segments, segmentLengths) = modules::volume_path_evolution(
		nodes,
		segments,
		segmentLengths,
        scene.h,
        descent
	);
	

    // visualization
    auto crv = polyscope::registerCurveNetwork("curve", nodes, segments);

    /*
	auto [dNodes, dSegments, dNorm] = formatVectorsForVisualization(nodes, descent);
	auto dn = polyscope::registerCurveNetwork("descent direction", descent, segments);
	dn->addEdgeScalarQuantity("descent norm", dNorm);
	auto gn = polyscope::registerCurveNetwork("gradient direction", gradient, segments);
	gn->addEdgeScalarQuantity("gradient norm", dNorm);
    */
}

void polyscopeCallback() {
    ImGui::Checkbox("run loop", &runLoop);

    if (ImGui::Button("run 1 step") || ImGui::IsKeyDown(' ') || runLoop) {
        doWork();
    }

    if (iteration % 100 == 0) {
        runLoop = false;
    }
}

int main(int argc, char **argv) {
	args::ArgumentParser parser("volume filling curves");
    args::Positional<std::string> inputFilename(parser, "config", "A config file.");
	
	// Parse args
    try {
        parser.ParseCLI(argc, argv);
    }
    catch (args::Help& h) {
        std::cout << parser;
        return 0;
    }
    catch (args::ParseError& e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }

    // Make sure a config file name was given
    if (!inputFilename) {
        std::cerr << "Please specify a config file as argument" << std::endl;
        return EXIT_FAILURE;
    }

    polyscope::init();
	polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;

    polyscope::state::userCallback = polyscopeCallback;
    scene = modules::read_scene(args::get(inputFilename));

    if (scene.rmax == 0) {
        scene.rmax = scene.radius * 10;
    }

    if (scene.volume.volumeType == scene_file::VolumeType::MESH) {
        std::cout << "Initializing OpenVDB" << std::endl;
        openvdb::initialize();

        std::cout << "Reading surface mesh" << std::endl;
        modules::PolyscopeMeshData polyscopeMesh = modules::file_to_polyscope_data(scene.volume.mesh_filename);
        std::cout << "Mesh data created with " << polyscopeMesh.vertices.size() << " vertices and "
            << polyscopeMesh.faces.size() << " faces." << std::endl;

		// visualize volume as mesh, if requested
		// note that this can be very slow for large meshes
        if (scene.visualizeVolume) {
            std::cout << "Registering polyscope mesh" << std::endl;
            polyscope::registerSurfaceMesh("mesh", polyscopeMesh.vertices, polyscopeMesh.faces)
                ->setTransparency(0.25f)
                ->setEnabled(true);
            std::cout << "Registered polyscope mesh" << std::endl;
        }


        if (scene.volume.convert_to_sdf) {
            // Use openvdb to convert mesh to SDF
            std::cout << "Converting Polyscope mesh data to OpenVDB format" << std::endl;
			modules::OpenVDBMeshData openvdbMesh = modules::polyscope_to_openvdb_data(polyscopeMesh);
            std::cout << "Converted mesh data with " << openvdbMesh.vertices.size() << " vertices and "
                << openvdbMesh.faces.size() << " faces." << std::endl;

            std::cout << "Generating SDF" << std::endl;
			openvdb::FloatGrid::Ptr sdf = modules::openvdb_mesh_to_sdf(openvdbMesh, scene.volume.mesh_to_sdf_voxelsize, scene.volume.mesh_to_sdf_halfwidth);
            std::cout << "SDF grid created with " << sdf->activeVoxelCount() << " active voxels." << std::endl;

			if (!modules::sdf_is_watertight(sdf)) {
				std::cerr << "SDF is not watertight, aborting." << std::endl;
				return EXIT_FAILURE;
			}

            if (!sdf) {
                std::cerr << "Error converting mesh to SDF" << std::endl;
                return EXIT_FAILURE;
            }

            // visualize sdf using a point cloud grid
            std::vector<Vector3> sdf_points;
            std::vector<float> sdf_values;

            for (auto iter = sdf->cbeginValueOn(); iter; ++iter) {
                openvdb::Vec3s pos = sdf->indexToWorld(iter.getCoord().asVec3s());
                sdf_points.push_back(Vector3{ pos.x(), pos.y(), pos.z() });
                sdf_values.push_back(iter.getValue() / sdf->background());
            }

            for (auto iter = sdf->cbeginValueOn(); iter; ++iter) {
                if (iter.getValue() < 0) {
                    openvdb::Coord ijk = iter.getCoord();
                    openvdb::Vec3d world_pos = sdf->indexToWorld(ijk);
                    std::cout << "Negative SDF at world: " << world_pos
                        << " value: " << iter.getValue() << std::endl;
                    break; // Just show first one
                }
            }

            polyscope::registerPointCloud("SDF", sdf_points)
                ->addScalarQuantity("SDF value", sdf_values)
                ->setEnabled(true);

            scene.volume.sdf = sdf;
        }
        else {
            scene.volume.mesh_points = polyscopeMesh.vertices;
        }

        // ToDo: Initialize nodes relative to loaded volume?
	}
	else if (scene.volume.volumeType == scene_file::VolumeType::SDF) {
	}

    if (scene.curveFileName != "") {
        std::tie(nodes, segments) = modules::read_nodes(scene.curveFileName);
    }

    // fallback in case nodes have not yet been initialized by any prior logic
    if (nodes.size() == 0) {
        // just create a circle (with noise) around the origin

        float circleRadius = 1.f;
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<float> noiseDist(-0.1f, 0.1f);

        /*
        for (int i = 0; i < 100; i++) {
            float angle = (float)i / 100 * 2 * igl::PI;
            float x = circleRadius * cos(angle) + noiseDist(gen);
            float y = circleRadius * sin(angle) + noiseDist(gen);
			float z = 0.f + noiseDist(gen);
            nodes.push_back(Vector3{ x, y, z });
        }
        */

        for (int i = 0; i < 100; i++) {
			float angle = (float)i / 100 * 2 * igl::PI;
            float x = circleRadius * cos(angle);
            float y = circleRadius * sin(angle);
            float z = circleRadius * sin(2*angle) * 0.3;

			nodes.push_back(Vector3{ x, y, z });
        }

        // if mesh, move nodes to mesh's center
        /*
		if (scene.volume.volumeType == scene_file::VolumeType::MESH && scene.volume.convert_to_sdf) {
			openvdb::Vec3s center = scene.volume.sdf->evalActiveVoxelBoundingBox().getCenter();
			for (auto& node : nodes) {
				node += Vector3{ center.x(), center.y(), center.z() };
			}
		}
        */

        // create segments
        for (int i = 0; i < nodes.size(); i++) {
            int j = (i + 1) % nodes.size();
            segments.push_back(std::array<int, 2>{i, j});
        }
    }

    // calculate segment lengths
	for (int i = 0; i < segments.size(); i++) {
		int j = segments[i][0];
		int k = segments[i][1];
		segmentLengths.push_back(norm(nodes[j] - nodes[k]));
	}

    // set initial values
	initialNodes = nodes;
	initialSegments = segments;
	initialSegmentLengths = segmentLengths;

	polyscope::registerCurveNetwork("initial curve", initialNodes, initialSegments)->setEnabled(false);
    polyscope::registerCurveNetwork("curve", nodes, segments);

    // ToDo: Calculate and render tangents (and normals, if that concept applies)

    polyscope::show();

    return EXIT_SUCCESS;
}
