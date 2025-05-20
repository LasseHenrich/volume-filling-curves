#include <tuple>
#include <vector>
#include <glm/glm.hpp>
#include <geometrycentral/utilities/vector3.h>

using namespace geometrycentral;

namespace modules {
	namespace VolumeFillingEnergy {

		struct Options {
			double p = 2; // exponent for length deviation penalty
			double q = 2; // exponent for direction alignment
			//double stepSizeFactor = 0.01; // scaling factor for step size
			//double w_medialAxis = 1.0; // weight for medial axis energy
			double w_bilaplacian = 0; // weight for bilaplacian energy
		};
	}

	std::tuple<
		std::vector<Vector3>, // descent direction
		std::vector<Vector3>, // gradient direction
		double, // energy
		std::vector<std::vector<Vector3>> // medial axis
	> volume_filling_energy(
		const std::vector<Vector3>& nodes,
		const std::vector<std::array<int, 2>>& segments,
		const std::vector<double>& segmentLengths,
		const double radius,
		const double maxRadius,
		const VolumeFillingEnergy::Options& options = {}
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
		const VolumeFillingEnergy::Options& options = {}
	);
}