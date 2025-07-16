#include <vector>
#include <geometrycentral/utilities/vector3.h>
#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/GridOperators.h>

using namespace geometrycentral;

namespace modules {
	/**
	 * Largely corresponds to medial_axis_euclidian in surface-filling-curves
	 */
	std::vector<std::vector<Vector3>> medial_axis_curve(
		const std::vector<Vector3>& nodes,
		const std::vector<std::array<int, 2>>& segments,
		const std::vector<Vector3>& nodeTangents,
		const std::vector<Vector3>& nodeNormals,
		const std::vector<Vector3>& nodeBitangents,
		const double maxRadius,
		const std::vector<Vector3>& meshPoints
	);

	std::vector<std::vector<Vector3>> medial_axis_surface(
		const std::vector<Vector3>& nodes,
		const std::vector<Vector3>& nodeNormals,
		const double maxRadius
	);

	void backproject_centers_to_hull(
		std::vector<std::vector<Vector3>>& nodeMedialAxis,
		const std::vector<Vector3>& nodes,
		openvdb::tools::GridSampler<openvdb::FloatGrid, openvdb::tools::BoxSampler>& sampler,
		double radius,
		double maxRadius
	);
}