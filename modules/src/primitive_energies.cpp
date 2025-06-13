#include <primitive_energies.h>

namespace modules {
	template <typename T, typename ElementType>
	T box_penalty(const ElementType& element, const Eigen::Vector3d& boxExtents, double totalCurveLength)
	{
		
	}

	template <typename ReturnType, typename ElementType>
	ReturnType primitive_penalty(const ElementType& element, scene_file::SceneObject_Volume& primitive)
	{
		if (primitive.primitive_type == PrimitiveType::SPHERE) {

		}
		else if (primitive.primitive_type == PrimitiveType::BOX) {
			if (primitive.primitive_params.size() < 3) {
				std::cerr << "Error: Box primitive requires 3 parameters (half-extents)" << std::endl;
				std::abort();
			}
			Eigen::Vector3d boxExtents(
				primitive.primitive_params[0],
				primitive.primitive_params[1],
				primitive.primitive_params[2]
			);
			double totalCurveLength = element.variables.size(); // Assuming this is the total length of the curve
			return box_penalty(element, boxExtents, totalCurveLength);
		}

		std::abort(); // Unsupported primitive type
	}
}