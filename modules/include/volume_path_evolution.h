#pragma once
#include <geometrycentral/utilities/vector3.h>
#include "scene_file.h"
#include "../../volume-filling-curves.h"
#include <vector>
#include <array>
#include <tuple>

using namespace geometrycentral;

namespace modules {
	std::tuple<
		std::vector<Vector3>, // nodes
		std::vector<std::array<int, 2>>, // segments
		std::vector<double> // segment lengths
	> volume_path_evolution_curve(
		const Curve& curve,
		const double h,
		const std::vector<Vector3>& descentDirections,
		const scene_file::SceneObject& options
	);

	void clamp_to_boundary(
		std::vector<Vector3>& nodes,
		const std::vector<Vector3>& old_nodes,
		double step,
		const std::vector<Vector3>& descentDirections,
		const scene_file::SceneObject& options
	);

	bool check_path_validity(
		const std::vector<Vector3>& nodes,
		const std::vector<std::array<int, 2>>& segments
	);
}
