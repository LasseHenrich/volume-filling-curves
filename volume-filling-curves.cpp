#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include <polyscope/curve_network.h>
#include <polyscope/point_cloud.h>
#include <polyscope/volume_mesh.h>

#include "args/args.hxx"

#include "volume-filling-curves.h"
#include <chrono>
#include <volume_filling_energy.h>
#include <volume_path_evolution.h>
#include <scene_file.h>

using namespace modules;

float radius = 30;
float timestep = 1;
float h = igl::PI * radius / 25;
float rmax = 0.5;

modules::SceneObject scene;

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

    // 0. initialize options object
    modules::VolumeFillingEnergy::Options options;  
    options.p = scene.p;
    options.q = scene.q;
    options.w_bilaplacian = scene.w_bilaplacian;
	// ToDo: Assign rest of options from scene file

    // 1. compute descent direction
	auto [descent, gradient, energy, medialAxis] = modules::volume_filling_energy(
		nodes,
		segments,
		segmentLengths,
		radius,
		rmax,
		options
	);

    // 2. evolve curve without self-intersections
	std::tie(nodes, segments, segmentLengths) = modules::volume_path_evolution(
		nodes,
		segments,
		segmentLengths,
        h,
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
    polyscope::state::userCallback = polyscopeCallback;
    scene = modules::read_scene(args::get(inputFilename));

    radius = scene.radius;
    h = scene.h;
    timestep = scene.timestep;
    rmax = scene.rmax == 0 ? radius * 10 : scene.rmax;

    // ToDo: Use specified volume

    if (scene.curveFileName != "") {
        std::tie(nodes, segments) = modules::read_nodes(scene.curveFileName);

        /*
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<float> noiseDist(-0.1f, 0.1f);

        // add noise to nodes
        for (int i = 0; i < nodes.size(); i++) {
            nodes[i].x += noiseDist(gen);
            nodes[i].y += noiseDist(gen);
            nodes[i].z += noiseDist(gen);
        }
        */
    }

    if (nodes.size() == 0) {
        float circleRadius = 1.f;
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<float> noiseDist(-0.1f, 0.1f);

        // just create a circle (with noise) around the origin
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

	polyscope::registerCurveNetwork("initial curve", initialNodes, initialSegments);
    polyscope::registerCurveNetwork("curve", nodes, segments);

    // ToDo: Calculate and render tangents (and normals, if that concept applies)

    polyscope::show();

    return EXIT_SUCCESS;
}
