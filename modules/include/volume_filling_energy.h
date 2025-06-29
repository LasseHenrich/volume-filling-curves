#include <tuple>
#include <vector>
#include <glm/glm.hpp>
#include <geometrycentral/utilities/vector3.h>
#include <scene_file.h>
#include <TinyAD/ScalarFunction.hh>
#include "../../volume-filling-curves.h"

using namespace geometrycentral;

namespace modules {
	std::tuple<
		std::vector<Vector3>, // descent direction
		std::vector<Vector3>, // gradient direction
		double, // energy
		std::vector<std::vector<Vector3>> // medial axis
	> volume_filling_energy_curve(
		const Curve& curve,
		const scene_file::SceneObject& options = {}
	);

	std::tuple<
		std::vector<Vector3>, // descent direction
		std::vector<Vector3>, // gradient direction
		double, // energy
		std::vector<std::vector<Vector3>> // medial axis
	> dirichlet_energy(
		const std::vector<Vector3>& nodes,
		const std::vector<std::array<int, 2>>& segments,
		const std::vector<double>& segmentLengths,
		const double radius,
		const double maxRadius,
		const scene_file::SceneObject& options = {}
	);

	Vector3 arbitrary_normal(const Vector3& t);

	template <typename ReturnType, typename ElementType>
	ReturnType volume_penalty(const ElementType& element, scene_file::SceneObject_Volume &surface) {

	}
}