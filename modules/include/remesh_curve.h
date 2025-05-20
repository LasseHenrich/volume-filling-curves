#include <tuple>
#include <geometrycentral/utilities/vector3.h>

using namespace geometrycentral;

namespace modules {
	std::tuple<
		std::vector<Vector3>, // nodes
		std::vector<std::array<int, 2>>, // segments
		std::vector<double> // segment lengths
	> remesh_curve(
		const std::vector<Vector3> nodes,
		const std::vector<std::array<int, 2>>& segments,
		const std::vector<double>& segmentLengths,
		const double h
	);
}