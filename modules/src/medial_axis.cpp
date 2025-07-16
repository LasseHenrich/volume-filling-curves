#include "medial_axis.h"
#include <knn-cpp/knncpp.h>
#include <chrono>
#include <cassert>

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

	std::vector<std::vector<Vector3>> medial_axis_curve(
		const std::vector<Vector3>& nodes,
		const std::vector<std::array<int, 2>>& segments,
		const std::vector<Vector3>& nodeTangents,
		const std::vector<Vector3>& nodeNormals,
		const std::vector<Vector3>& nodeBitangents,
		const double maxRadius,
		const std::vector<Vector3>& meshPoints
	) {
		auto kdTreeStart = std::chrono::high_resolution_clock::now();

		std::vector<Vector3> allpoints = nodes;
		std::cout << "Number of nodes: " << nodes.size() << std::endl;
		allpoints.insert(allpoints.end(), meshPoints.begin(), meshPoints.end());
		std::cout << "Number of all points: " << allpoints.size() << std::endl;

		// create a closest point query
		Eigen::MatrixXd _V(3, allpoints.size());
		for (int i = 0; i < allpoints.size(); i++) {
			_V.col(i) << allpoints[i].x, allpoints[i].y, allpoints[i].z;
		}
		knncpp::KDTreeMinkowskiX<double, knncpp::EuclideanDistance<double>> kdtree(_V);
		kdtree.build();

		auto kdTreeEnd = std::chrono::high_resolution_clock::now();

		std::vector<std::vector<Vector3>> nodeMedialAxis(nodes.size());

		#pragma omp parallel for
		for (int i = 0; i < nodes.size(); i++) {
			// create a plane perpendicular to the tangent at this node
			Vector3 t = nodeTangents[i];
			Vector3 n = nodeNormals[i];
			Vector3 b = nodeBitangents[i];
			Vector3 x = nodes[i];

			double r_min_minus_b = std::numeric_limits<double>::infinity();
			double r_min_plus_b = std::numeric_limits<double>::infinity();
			double r_min_minus_n = std::numeric_limits<double>::infinity();
			double r_min_plus_n = std::numeric_limits<double>::infinity();

			r_min_minus_b = std::min(maximumBallRadius(x, -b, i, maxRadius, allpoints, kdtree), maxRadius);
			r_min_plus_b = std::min(maximumBallRadius(x, b, i, maxRadius, allpoints, kdtree), maxRadius);
			r_min_minus_n = std::min(maximumBallRadius(x, -n, i, maxRadius, allpoints, kdtree), maxRadius);
			r_min_plus_n = std::min(maximumBallRadius(x, n, i, maxRadius, allpoints, kdtree), maxRadius);

			nodeMedialAxis[i].emplace_back(nodes[i] + r_min_minus_b * -b);
			nodeMedialAxis[i].emplace_back(nodes[i] + r_min_plus_b * b);
			nodeMedialAxis[i].emplace_back(nodes[i] + r_min_minus_n * -n);
			nodeMedialAxis[i].emplace_back(nodes[i] + r_min_plus_n * n);
		}

		return nodeMedialAxis;
	}

	std::vector<std::vector<Vector3>> medial_axis_surface(
		const std::vector<Vector3>& nodes,
		const std::vector<Vector3>& nodeNormals,
		const double maxRadius
	) {
		auto kdTreeStart = std::chrono::high_resolution_clock::now();

		// create a closest point query
		Eigen::MatrixXd _V(3, nodes.size());
		for (int i = 0; i < nodes.size(); i++) {
			_V.col(i) << nodes[i].x, nodes[i].y, nodes[i].z;
		}
		knncpp::KDTreeMinkowskiX<double, knncpp::EuclideanDistance<double>> kdtree(_V);
		kdtree.build();

		auto kdTreeEnd = std::chrono::high_resolution_clock::now();

		std::vector<std::vector<Vector3>> nodeMedialAxis(nodes.size());

		#pragma omp parallel for
		for (int i = 0; i < nodes.size(); i++) {
			Vector3 n = nodeNormals[i];
			Vector3 x = nodes[i];

			double r_min_minus_n = std::numeric_limits<double>::infinity();
			double r_min_plus_n = std::numeric_limits<double>::infinity();
			
			r_min_plus_n = std::min(maximumBallRadius(x, n, i, maxRadius, nodes, kdtree), maxRadius);
			r_min_minus_n = std::min(maximumBallRadius(x, -n, i, maxRadius, nodes, kdtree), maxRadius);
			
			nodeMedialAxis[i].emplace_back(nodes[i] - r_min_minus_n * n);
			nodeMedialAxis[i].emplace_back(nodes[i] + r_min_plus_n * n);
		}

		return nodeMedialAxis;
	}

	void backproject_centers_to_hull(
		std::vector<std::vector<Vector3>>& nodeMedialAxis,
		const std::vector<Vector3>& nodes,
		openvdb::tools::GridSampler<openvdb::FloatGrid, openvdb::tools::BoxSampler>& sampler,
		double radius,
		double maxRadius
	) {
		const size_t numNodes = nodeMedialAxis.size();

		#pragma omp parallel for
		for (int nodeId = 0; nodeId < numNodes; ++nodeId) {
			Vector3 pos = nodes[nodeId];

			std::vector<Vector3>& centers = nodeMedialAxis[nodeId];
			const size_t numPoints = centers.size();

			std::vector<Vector3> u_vectors;
			std::vector<double> lengths;

			u_vectors.reserve(numPoints);
			lengths.reserve(numPoints);

			for (size_t i = 0; i < numPoints; ++i) {
				Vector3 u = pos - centers[i];
				u_vectors.push_back(u);
				lengths.push_back(norm(u));
			}

			openvdb::Vec3d pos_ovdb(pos.x, pos.y, pos.z);

			// Sample SDF value at current position
			double sdf_value = sampler.wsSample(pos_ovdb);

			// Compute gradient using finite differences
			const double h = 1e-6;
			double sdf_px = sampler.wsSample(pos_ovdb + openvdb::Vec3d(h, 0, 0));
			double sdf_mx = sampler.wsSample(pos_ovdb - openvdb::Vec3d(h, 0, 0));
			double sdf_py = sampler.wsSample(pos_ovdb + openvdb::Vec3d(0, h, 0));
			double sdf_my = sampler.wsSample(pos_ovdb - openvdb::Vec3d(0, h, 0));
			double sdf_pz = sampler.wsSample(pos_ovdb + openvdb::Vec3d(0, 0, h));
			double sdf_mz = sampler.wsSample(pos_ovdb - openvdb::Vec3d(0, 0, h));

			double grad_x = (sdf_px - sdf_mx) / (2.0 * h);
			double grad_y = (sdf_py - sdf_my) / (2.0 * h);
			double grad_z = (sdf_pz - sdf_mz) / (2.0 * h);

			double grad_length = std::sqrt(grad_x * grad_x + grad_y * grad_y + grad_z * grad_z);
			if (grad_length < 1e-6) {
				continue;
			}

			// gradient should already have length 1,
			// but we normalize it just in case
			grad_x /= grad_length;
			grad_y /= grad_length;
			grad_z /= grad_length;

			std::vector<double> sdf_approx_at_centers;
			sdf_approx_at_centers.reserve(numPoints);

			for (size_t i = 0; i < numPoints; ++i) {
				const auto& c = centers[i];
				double sdf_approx = sdf_value +
					grad_x * (c.x - pos.x) +
					grad_y * (c.y - pos.y) +
					grad_z * (c.z - pos.z);
				sdf_approx_at_centers.push_back(sdf_approx);
			}

			// note that the gradient already has length 1
			Vector3 closest_surface_point = pos - // sdf is negative inside, so '-'
				sdf_approx_at_centers[0] * Vector3{ grad_x, grad_y, grad_z };

			Vector3 v = pos - closest_surface_point;
			if (sdf_value > 0) {
				v = -v; // flip direction if outside
			}

			auto backproject = [&](double sdf_at_c, Vector3 u, double l) -> double {
				double b = sdf_at_c - radius;
				if (b <= 0) {
					return l;
				}

				double b_tilde = (b * u.norm() * v.norm()) / dot(u, v);
				if (dot(u, v) <= 0) { // if non-positive, the angle is obtuse, so u is not pointing towards the surface
					return l;
				}

				return std::max(
					std::min(l - b_tilde, maxRadius),
					0.0
				);
			};

			for (size_t i = 0; i < numPoints; ++i) {
				double prev_length = lengths[i];
				double new_length = backproject(sdf_approx_at_centers[i], u_vectors[i], lengths[i]);

				// move the center point to the new position
				if (std::abs(new_length - prev_length) > 1e-6) {
					Vector3 u = u_vectors[i];
					nodeMedialAxis[nodeId][i] = pos - new_length * (u / u.norm());
				}
			}
		}
	}
}