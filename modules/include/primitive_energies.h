#include <TinyAD/ScalarFunction.hh>
#include <scene_file.h>

using namespace geometrycentral;

namespace modules {
	template <typename ReturnType, typename ElementType>
	ReturnType box_penalty(const ElementType& element, const Eigen::Vector3d& boxExtents, double totalCurveLength);

	template <typename ReturnType, typename ElementType>
	ReturnType primitive_penalty(const ElementType& element, scene_file::SceneObject_Volume& primitive);
}