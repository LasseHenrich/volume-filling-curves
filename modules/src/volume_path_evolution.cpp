#include "volume_path_evolution.h"
#include <tuple>
#include <vector>
#include <geometrycentral/utilities/vector3.h>
#include <chrono>
#include <remesh_curve.h>

using namespace geometrycentral;

namespace modules {
	const int max_iters = 100;
	const double shrink = .8;

	std::tuple<
		std::vector<Vector3>, // nodes
		std::vector<std::array<int, 2>>, // segments
		std::vector<double> // segment lengths
	> volume_path_evolution(
		const std::vector<Vector3>& nodes,
		const std::vector<std::array<int, 2>>& segments,
		const std::vector<double>& segmentLengths,
		const double h,
		const std::vector<Vector3>& descentDirections
	) {
		auto newNodes = nodes;
		auto newSegments = segments;
		auto newSegmentLengths = segmentLengths;

		auto start = std::chrono::high_resolution_clock::now();

		double step = 1.0f;
		
		for (int i = 0; i < max_iters; i++) {
			// 1. Move nodes according to the descent direction
			for (size_t j = 0; j < nodes.size(); j++) {
				newNodes[j] = nodes[j] + step * descentDirections[j];
			}

			// 2. Recalculate segments lengths
			for (size_t j = 0; j < segments.size(); j++) {
				const auto& segment = segments[j];
				int startIdx = segment[0];
				int endIdx = segment[1];
				newSegmentLengths[j] = norm(newNodes[endIdx] - newNodes[startIdx]);
			}

			// 3. Remesh
			std::tie(newNodes, newSegments, newSegmentLengths) = modules::remesh_curve(
				newNodes,
				newSegments,
				newSegmentLengths,
				h
			);

			// 4. Check for self-intersections (ToDo: Make more sophisticated)
			bool isValid = check_path_validity(newNodes, newSegments);

			if (isValid) {
				return std::make_tuple(newNodes, newSegments, newSegmentLengths);
			}

			step *= shrink;
		}

		// if we reach here, it means we couldn't find a valid path,
		// so we return the original path
		return std::make_tuple(nodes, segments, segmentLengths);
	}

	bool check_path_validity(
		const std::vector<Vector3>& nodes,
		const std::vector<std::array<int, 2>>& segments)
	{
		// ToDo: Implement a function to check for self-intersections
		return true;
	}
}