#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include <polyscope/curve_network.h>
#include <polyscope/point_cloud.h>
#include <polyscope/volume_mesh.h>

#include "args/args.hxx"

#include <openvdb/openvdb.h>
#include <openvdb/tools/Interpolation.h>

#include <chrono>
#include <volume_filling_energy.h>
#include <volume_path_evolution.h>
#include <scene_file.h>
#include <mesh_parser.h>
#include <structures.h>

using namespace modules;

scene_file::SceneObject scene;
Curve currentCurve, restCurve, initialCurve;
Surface currentSurface, restSurface, initialSurface; // 'current'Surface as surface throws an ambiguity with polyscope::SurfaceMesh

int iteration = 0;
bool runLoop = false;

Curve formatVectorsForVisualization(
    const Curve& curve,
    std::vector<Vector3>& directions
) {
	assert(curve.nodes.size() == directions.size(), "nodes and directions must have the same size");

	std::vector<Vector3> dNodes = curve.nodes;
    std::vector<std::array<int, 2>> dSegments = {};
    std::vector<double> dNorm = {};
    
    for (int i = 0; i < directions.size(); i++) {
        dNodes.push_back(curve.nodes[i] + directions[i]);
		dNorm.push_back(norm(directions[i]));
    }

    for (int i = 0; i < directions.size(); i++) {
		int oi = i + directions.size();
		std::array<int, 2> line = { i, oi };
		dSegments.push_back(line);
	}

	return {
		dNodes,
		dSegments,
		dNorm
	};
}

void doWork_curve() {
    restCurve.nodes = currentCurve.nodes;
    restCurve.segments = currentCurve.segments;
    restCurve.segmentLengths = currentCurve.segmentLengths;

    // 1. compute descent direction
    auto [descent, gradient, energy, medialAxis] = modules::volume_filling_energy_curve(
        currentCurve,
        scene
    );

    // 2. evolve curve without self-intersections
    std::tie(currentCurve.nodes, currentCurve.segments, currentCurve.segmentLengths) = modules::volume_path_evolution_curve(
        currentCurve,
        scene.h,
        descent,
        scene
    );


    // visualization
    polyscope::registerCurveNetwork("curve", currentCurve.nodes, currentCurve.segments);
}

void doWork_surface() {
	// ToDo: Clone the current surface to the rest surface
	/*restSurface.vertexPositions = currentSurface.vertexPositions;
	restSurface.mesh = currentSurface.mesh->clone();*/

	// 1. compute descent direction
	auto [descent, gradient, energy, medialAxis] = modules::volume_filling_energy_surface(
		currentSurface,
		scene
	);

	// 2. evolve surface without self-intersections
	modules::volume_path_evolution_surface(
		currentSurface,
		scene.h,
		descent,
		scene
	);

	// visualization
	modules::PolyscopeMeshData currentSurface_polyscope = modules::geometrycentral_to_polyscope_data(&currentSurface);
    polyscope::registerSurfaceMesh("surface", currentSurface_polyscope.vertices, currentSurface_polyscope.faces);
}

void doWork() {
    iteration++;


    auto start = std::chrono::high_resolution_clock::now();

    std::cout << "===== iteration: " << iteration << "=====" << std::endl;

    if (scene.filling_dimension == 1) {
        std::cout << "numNodes: " << currentCurve.nodes.size() << std::endl;
        std::cout << "numSegments: " << currentCurve.segments.size() << std::endl;
        std::cout << std::endl;
        doWork_curve();
    }
    else {
        std::cout << "numVertices: " << currentSurface.mesh->nVertices() << std::endl;
        std::cout << "numFaces: " << currentSurface.mesh->nFaces() << std::endl;
        std::cout << std::endl;
        doWork_surface();
    }
}

void polyscopeCallback() {
    ImGui::Checkbox("run loop", &runLoop);

    if (ImGui::Button("run 1 step") || ImGui::IsKeyDown(' ') || runLoop) {
        doWork();
    }

    if (iteration % 200 == 0) {
        runLoop = false;
    }
}

void initialize_curve(Curve& curve, const scene_file::SceneObject& scene) {
    if (scene.fillingManifoldFileName != "") {
        std::tie(curve.nodes, curve.segments) = modules::read_curve(scene.fillingManifoldFileName);
    }

    // fallback in case nodes have not yet been initialized by any prior logic
    if (curve.nodes.size() == 0) {
        // just create a circle around the origin

        float circleRadius = 1.f;
        for (int i = 0; i < 100; i++) {
			float angle = (float)i / 100 * 2 * igl::PI;
            float x = circleRadius * cos(angle);
            float y = circleRadius * sin(angle);
            float z = circleRadius * sin(2*angle) * 0.3;

			curve.nodes.push_back(Vector3{ x, y, z });
        }

        // create segments
        for (int i = 0; i < curve.nodes.size(); i++) {
            int j = (i + 1) % curve.nodes.size();
            curve.segments.push_back(std::array<int, 2>{i, j});
        }
    }

    // calculate segment lengths
	for (int i = 0; i < curve.segments.size(); i++) {
		int j = curve.segments[i][0];
		int k = curve.segments[i][1];
		curve.segmentLengths.push_back(norm(curve.nodes[j] - curve.nodes[k]));
	}

    // set initial values
    initialCurve.nodes = curve.nodes;
    initialCurve.segments = curve.segments;
    initialCurve.segmentLengths = curve.segmentLengths;

    polyscope::registerCurveNetwork("initial curve", initialCurve.nodes, initialCurve.segments)->setEnabled(false);
    polyscope::registerCurveNetwork("curve", curve.nodes, curve.segments);
}


void initialize_surface(Surface& surface, const scene_file::SceneObject& scene) {
    if (scene.fillingManifoldFileName == "") {
		std::cerr << "No initial surface mesh file specified in the scene configuration." << std::endl;
        std::abort();
    }

	surface = modules::file_to_geometrycentral_data(scene.fillingManifoldFileName);

    // we just read it two times, as the copying is difficult and I failed to find a solutino
	initialSurface = modules::file_to_geometrycentral_data(scene.fillingManifoldFileName);

	modules::PolyscopeMeshData currentSurface_polyscope = modules::geometrycentral_to_polyscope_data(&surface);
	modules::PolyscopeMeshData initialSurface_polyscope = modules::geometrycentral_to_polyscope_data(&initialSurface);
	polyscope::registerSurfaceMesh("initial surface", initialSurface_polyscope.vertices, initialSurface_polyscope.faces)->setEnabled(false);
	polyscope::registerSurfaceMesh("surface", currentSurface_polyscope.vertices, currentSurface_polyscope.faces);
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
		modules::Surface geometrycentralMeshData = modules::file_to_geometrycentral_data(scene.volume.mesh_filename);
		std::cout << "Surface mesh read with " << geometrycentralMeshData.mesh->nVertices() << " vertices and "
			<< geometrycentralMeshData.mesh->nFaces() << " faces." << std::endl;

		std::cout << "Converting to Polyscope mesh data" << std::endl;
		modules::PolyscopeMeshData polyscopeMesh = modules::geometrycentral_to_polyscope_data(&geometrycentralMeshData);
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
			modules::OpenVDBMeshData openvdbMesh = modules::polyscope_to_openvdb_data(&polyscopeMesh);
            std::cout << "Converted mesh data with " << openvdbMesh.vertices.size() << " vertices and "
                << openvdbMesh.faces.size() << " faces." << std::endl;

			if (!modules::mesh_is_manifold(&geometrycentralMeshData.mesh)) {
				std::cerr << "Mesh is not watertight, aborting." << std::endl;
				return EXIT_FAILURE;
			}

            std::cout << "Generating SDF" << std::endl;
			openvdb::FloatGrid::Ptr sdf = modules::openvdb_mesh_to_sdf(&openvdbMesh, scene.volume.mesh_to_sdf_voxelsize, scene.volume.mesh_to_sdf_halfwidth);
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

            auto sdf_vis = polyscope::registerPointCloud("SDF", sdf_points);
            sdf_vis->addScalarQuantity("SDF value", sdf_values);
            sdf_vis->setEnabled(false);

            scene.volume.sdf = sdf;
        }
        else {
            scene.volume.mesh_points = polyscopeMesh.vertices;
        }

        // ToDo: Initialize nodes relative to loaded volume?
	}
	else if (scene.volume.volumeType == scene_file::VolumeType::SDF) {
	}

    if (scene.filling_dimension == 1) {
        initialize_curve(currentCurve, scene);
    }
    else {
        initialize_surface(currentSurface, scene);
    }

    polyscope::show();

    return EXIT_SUCCESS;
}