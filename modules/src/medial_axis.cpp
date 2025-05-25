#include "medial_axis.h"
#include <knn-cpp/knncpp.h>

namespace modules {
	int closestPointIndex(
		const Vector3& p,
		const knncpp::KDTreeMinkowskiX<double, knncpp::EuclideanDistance<double>>& kdtree
	) {
		Eigen::MatrixXd P(3, 1);
		Eigen::MatrixXd distances;
		knncpp::Matrixi indices;
		P << p.x, p.y, p.z;
		kdtree.query(P, 1, indices, distances);

		assert(indices.rows() == 1, "KDTree query failed");

		return indices(0, 0);
	}

	double maximumBallRadius(
		const Vector3& x,
		const Vector3& b,
		const int i,
		const double maxRadius,
		const std::vector<Vector3>& nodes,
		knncpp::KDTreeMinkowskiX<double, knncpp::EuclideanDistance<double>>& kdtree
	) {
		auto r = maxRadius;
		auto c = x + r * b;

		int nn = closestPointIndex(c, kdtree);
		bool finished = nn == i;

		double bsMax = 1.0, bsMin = 0.0;
		int itrc = 0;

		while (!finished) {
			itrc++;

			r = maxRadius * (bsMax + bsMin) / 2.0;

			auto c = x + r * b;
			int nn = closestPointIndex(c, kdtree);

			if (nn == i) {
				bsMin = (bsMax + bsMin) / 2;
			}
			else {
				// check if there is a point inside the ball
				Vector3 xy = nodes[nn] - nodes[i];
				r = pow(norm(xy), 2) / (2 * dot(xy, b));

				auto c = x + r * b;
				int _nn = closestPointIndex(c, kdtree);

				if (_nn == nn || _nn == i) {
					// if no one inside, break
					finished = true;
				}
				else {
					// else, continue the binary search
					bsMax = (bsMax + bsMin) / 2;
					assert(bsMax > bsMin);
				}
			}

			if (itrc > 100) {
				break;
			}
		}

		return r;
	}

	std::vector<std::vector<Vector3>> medial_axis(
		const std::vector<Vector3>& nodes,
		const std::vector<std::array<int, 2>>& segments,
		const std::vector<Vector3>& nodeTangents,
		const std::vector<Vector3>& nodeNormals,
		const std::vector<Vector3>& nodeBitangents,
		const double maxRadius
	) {
		// create a closes point query
		Eigen::MatrixXd _V(3, nodes.size());
		for (int i = 0; i < nodes.size(); i++) {
			_V.col(i) << nodes[i].x, nodes[i].y, nodes[i].z;
		}
		knncpp::KDTreeMinkowskiX<double, knncpp::EuclideanDistance<double>> kdtree(_V);
		kdtree.build();

		std::vector<std::vector<Vector3>> nodeMedialAxis(nodes.size());

		for (int i = 0; i < nodes.size(); i++) {
			// create a plane perpendicular to the tangent at this node
			Vector3 t = nodeTangents[i];
			Vector3 n = nodeNormals[i];
			Vector3 b = nodeBitangents[i];
			Vector3 x = nodes[i];

			double r_min_minus = std::numeric_limits<double>::infinity();
			double r_min_plus = std::numeric_limits<double>::infinity();

			r_min_plus = std::min(maximumBallRadius(x, b, i, maxRadius, nodes, kdtree), maxRadius);
			r_min_minus = std::min(maximumBallRadius(x, -b, i, maxRadius, nodes, kdtree), maxRadius);

			nodeMedialAxis[i].emplace_back(nodes[i] - r_min_minus * b);
			nodeMedialAxis[i].emplace_back(nodes[i] + r_min_plus * b);
		}

		return nodeMedialAxis;
	}
}