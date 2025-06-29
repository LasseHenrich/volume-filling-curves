#include "volume_path_evolution.h"
#include <tuple>
#include <vector>
#include <geometrycentral/utilities/vector3.h>
#include <chrono>
#include <remesh_curve.h>
#include <openvdb/tools/Interpolation.h>

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
		const std::vector<Vector3>& descentDirections,
		const scene_file::SceneObject& options
	) {
		auto newNodes = nodes;
		auto newSegments = segments;
		auto newSegmentLengths = segmentLengths;

		auto start = std::chrono::high_resolution_clock::now();

		double step = 1.0f;
		
		for (int i = 0; i < max_iters; i++) {
			// 1. Move nodes according to the descent direction and clamp to boundary
			clamp_to_boundary(newNodes, nodes, step, descentDirections, options);

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

	void clamp_to_boundary(
		std::vector<Vector3>& newNodes,
		const std::vector<Vector3>& old_nodes,
		double step,
		const std::vector<Vector3>& descentDirections,
		const scene_file::SceneObject& options) {

		// If the volume is not represented by an SDF, skip clamping.
		if (!options.volume.convert_to_sdf) {
			for (size_t j = 0; j < old_nodes.size(); j++) {
				newNodes[j] = old_nodes[j] + step * descentDirections[j];
			}
			return;
		}

		openvdb::tools::GridSampler<openvdb::FloatGrid, openvdb::tools::BoxSampler> sampler(*options.volume.sdf);

		for (size_t j = 0; j < old_nodes.size(); j++) {
			Vector3 p = old_nodes[j];
			Vector3 v = step * descentDirections[j];
			Vector3 p_target = p + v;

			double d_target = sampler.wsSample(openvdb::Vec3d(p_target.x, p_target.y, p_target.z));

			// if the target point is outside the volume (SDF > 0), find the intersection point.
			if (options.use_backprojection && d_target > 0) {
				// bisection search to find the parameter t in [0, 1] such that p + t*v is on the boundary.
				double t_low = 0.0;
				double t_high = 1.0;

				// 10 iterations should be sufficient for good precision... could be adjusted based on the step size though
				for (int k = 0; k < 10; ++k) {
					double t_mid = (t_low + t_high) / 2.0;
					Vector3 p_mid = p + t_mid * v;
					double d_mid = sampler.wsSample(openvdb::Vec3d(p_mid.x, p_mid.y, p_mid.z));
					if (d_mid < 0) { // midpoint is inside, so the intersection is in the upper half of the interval.
						t_low = t_mid;
					} else { // midpoint is outside or on the boundary, so the intersection is in the lower half.
						t_high = t_mid;
					}
				}
				newNodes[j] = p + t_low * v;
			} else {
				newNodes[j] = p_target;
			}
		}
	}

	bool check_path_validity(
		const std::vector<Vector3>& nodes,
		const std::vector<std::array<int, 2>>& segments)
	{
		// ToDo: Implement a function to check for self-intersections
		return true;
	}
}
