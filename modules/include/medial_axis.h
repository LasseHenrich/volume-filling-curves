#include <vector>
#include <geometrycentral/utilities/vector3.h>

using namespace geometrycentral;

namespace modules {
	/**
	 * Largely corresponds to medial_axis_euclidian in surface-filling-curves
	 */
	std::vector<std::vector<Vector3>> medial_axis(
		const std::vector<Vector3>& nodes,
		const std::vector<std::array<int, 2>>& segments,
		const std::vector<Vector3>& nodeTangents,
		const std::vector<Vector3>& nodeNormals,
		const std::vector<Vector3>& nodeBitangents,
		const double maxRadius
	);
}