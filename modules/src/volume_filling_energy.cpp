#include "volume_filling_energy.h"
#include <chrono>
#include <iostream>
#include <TinyAD/ScalarFunction.hh>
#include <glm/glm.hpp>
#include <Eigen/Dense>
#include <array>
#include <medial_axis.h>
#include <TinyAD/Utils/LinearSolver.hh>
#include <TinyAD/Utils/NewtonDirection.hh>

namespace modules {
	const double branchRatio = std::sqrt(std::sqrt(2));

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
		const VolumeFillingEnergy::Options& options
	) {
		auto start = std::chrono::high_resolution_clock::now();

		double p = options.p;
		double q = options.q;
		double branchRadius = radius * branchRatio;
		double alpha = 4 / (pow(branchRadius, 2));

		std::cout << "alpha: " << alpha << ", p: " << p << ", q: " << q << std::endl;

		std::vector<Vector3> descent(nodes.size(), Vector3{ .0, .0, .0 }), gradient(nodes.size(), Vector3{ .0, .0, .0 });

		double totalCurveLength = .0;
		for (int i = 0; i < segments.size(); i++) {
			totalCurveLength += segmentLengths[i];
		}

		auto func = TinyAD::scalar_function<3>(TinyAD::range(nodes.size()));

		std::vector<Eigen::Matrix3d> rotationMatrix(segments.size());
		std::vector<Vector3> segmentTangents(segments.size());
		std::vector<Vector3> segmentNormal(segments.size());
		std::vector<Vector3> segmentBitangent(segments.size());
		
		// 1. Compute tangents (verified CORRECT with HARDCODED normal)
		for (size_t i = 0; i < segments.size(); i++) {
			int j = segments[i][0];
			int k = segments[i][1];

			Vector3 t = normalize(nodes[k] - nodes[j]);

			/*
			Vector3 n;
			if (std::abs(t.z) > 0.707) {
				// if tangent is more along z-axis, use x-axis for the cross product
				n = normalize(cross(Vector3{ 1, 0, 0 }, t));
			}
			else {
				// otherwise, use z-axis for the cross product
				n = normalize(cross(Vector3{ 0, 0, 1 }, t));
			}

			// ensure normal is perpendicular to tangent
			n = normalize(n - dot(n, t) * t);

			// note that with this, we sometimes have a sign flip
			// in comparison to surface-filling-curves...
			*/
			Vector3 n = Vector3{ 0, 1, 0 }; // DEBUG for 2D case. ToDo

			Vector3 b = cross(n, t);

			segmentTangents[i] = t;
			segmentNormal[i] = n;
			segmentBitangent[i] = b;

			Eigen::Matrix3d R;
			R << t.x, b.x, n.x,
				t.y, b.y, n.y,
				t.z, b.z, n.z;

			rotationMatrix[i] = R;
		}

		//std::cout << "field alignedness: " << options.w_fieldAlignedness << ", curvature alignedness: " << options.w_curvatureAlignedness << ", bilaplacian weight: " << options.w_bilaplacian << std::endl;
		// no support for these at the moment, all zero

		// 2. Dirichlet term (verified CORRECT)
		func.add_elements<2>(TinyAD::range(nodes.size()), [&](auto& element)->TINYAD_SCALAR_TYPE(element) {
			using T = TINYAD_SCALAR_TYPE(element);
			int segmentId = element.handle;

			int _v0 = segments[segmentId][0];
			int _v1 = segments[segmentId][1];

			Eigen::Vector3<T> x0 = element.variables(_v0);
			Eigen::Vector3<T> x1 = element.variables(_v1);

			Eigen::Vector3d v0(
				nodes[_v0].x,
				nodes[_v0].y,
				nodes[_v0].z
			);

			Eigen::Vector3d v1(
				nodes[_v1].x,
				nodes[_v1].y,
				nodes[_v1].z
			);

			double edgeLen = segmentLengths[segmentId];

			Eigen::Matrix3d R = rotationMatrix[segmentId];

			Eigen::Vector3<T> p0 = R.transpose() * (x0 - v0);
			Eigen::Vector3<T> p1 = R.transpose() * (x1 - v1) + Eigen::Vector3d(edgeLen, 0, 0);

			T dx = abs(p0(0) - p1(0));
			T dy = abs(p0(1) - p1(1));
			T dz = abs(p0(2) - p1(2));

			// crs0 & crs1 as well as dot0 & dot1 are zero
			// without field-aligned and curvature energies

			auto result = (
				pow(pow(dx, 2) + pow(dy, 2) + pow(dz, 2), p / 2)
				) / (edgeLen * totalCurveLength);

			return result;
		});

		auto dirichletEnd = std::chrono::high_resolution_clock::now();

		// 3. Medial axis term
		auto tetTopologyEnd = std::chrono::high_resolution_clock::now();

		// 3.1 Node-based tangents (verified CORRECT)
		std::vector<Vector3> nodeTangents(nodes.size());
		std::vector<Vector3> nodeNormals(nodes.size());
		std::vector<Vector3> nodeBitangents(nodes.size());
		std::vector<double> nodeWeight(nodes.size(), .0);
		for (int i = 0; i < segments.size(); i++) {
			for (int j = 0; j < 2; j++) {
				auto t = segmentTangents[i];
				auto n = segmentNormal[i];
				auto b = segmentBitangent[i];

				int v = segments[i][j];

				nodeTangents[v] += normalize(t);
				nodeNormals[v] += normalize(n);
				nodeBitangents[v] += normalize(b);

				if (norm(nodeTangents[v]) < .00001) {
					nodeTangents[v] -= .01 * normalize(n);
				}
				if (norm(nodeNormals[v]) < .00001) {
					nodeNormals[v] -= .01 * normalize(t);
				}
				if (norm(nodeBitangents[v]) < .00001) {
					nodeBitangents[v] -= .01 * normalize(t);
				}

				nodeWeight[v] += segmentLengths[i] * 0.5;
			}
		}

		for (int i = 0; i < nodes.size(); i++) {
			nodeTangents[i] = normalize(nodeTangents[i]);
			nodeNormals[i] = normalize(nodeNormals[i]);
			nodeBitangents[i] = normalize(nodeBitangents[i]);
		}

		auto tangentEvalEnd = std::chrono::high_resolution_clock::now();

		// 3.2 Medial axis (verified CORRECT)
		auto nodeMedialAxis = modules::medial_axis(
			nodes,
			segments,
			nodeTangents,
			nodeNormals,
			nodeBitangents,
			maxRadius
		);

		auto medialAxisEnd = std::chrono::high_resolution_clock::now();

		// 3.3 Add medial axis energy term
		// note that with our current setup we just have one uniform alpha
		func.add_elements<1>(TinyAD::range(nodes.size()), [&](auto& element)->TINYAD_SCALAR_TYPE(element) {
			using T = TINYAD_SCALAR_TYPE(element);
			int nodeId = element.handle;

			auto _c_0 = nodeMedialAxis[nodeId][0];
			auto _c_1 = nodeMedialAxis[nodeId][1];

			Eigen::Vector3d c_0(_c_0.x, _c_0.y, _c_0.z);
			Eigen::Vector3d c_1(_c_1.x, _c_1.y, _c_1.z);

			T l_0 = (element.variables(nodeId) - c_0).norm();
			T l_1 = (element.variables(nodeId) - c_1).norm();

			auto result = alpha * nodeWeight[nodeId] * (
				pow(pow(l_0, 2) + pow(l_1, 2), q / 2)
				) / totalCurveLength;

			return result;
		});

		auto repulsiveEnd = std::chrono::high_resolution_clock::now();

		auto x = func.x_from_data([&](int v_idx) {
			auto v = nodes[v_idx];
			return Eigen::Vector3d(v.x, v.y, v.z);
			});

		TinyAD::LinearSolver solver;
		// calls all lambda functions that were added via func.add_elements
		auto [f, g, H_proj] = func.eval_with_hessian_proj(x);

		auto hessianEnd = std::chrono::high_resolution_clock::now();

		auto d = TinyAD::newton_direction(g, H_proj, solver);

		auto newtonEnd = std::chrono::high_resolution_clock::now();

		func.x_to_data(d, [&](int idx, const Eigen::Vector3d& p) {
			descent[idx] = Vector3 { p(0), p(1), p(2) };
			});

		func.x_to_data(g, [&](int idx, const Eigen::Vector3d& p) {
			gradient[idx] = Vector3{ p(0), p(1), p(2) };
			});

		auto end = std::chrono::high_resolution_clock::now();

		// verified CORRECT until here

		/*
		Log time recordings...
		*/

		return std::make_tuple(
			descent,
			gradient,
			f,
			nodeMedialAxis
		);
	}

	/*
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
		const VolumeFillingEnergy::Options& options
	) {
		auto start = std::chrono::high_resolution_clock::now();

		double p = options.p;
		double q = options.q;
		double branchRadius = radius * branchRatio;
		double alpha = 4.0 / (branchRadius * branchRadius);

		std::cout << "alpha: " << alpha << ", p: " << p << ", q: " << q << std::endl;

		std::vector<Vector3> descent(nodes.size(), Vector3());
		std::vector<Vector3> gradient(nodes.size(), Vector3());

		double energy = 0.0;
		double totalCurveLength = 0.0;

		for (int i = 0; i < segments.size(); i++) {
			totalCurveLength += segmentLengths[i];
		}

		std::vector<Vector3> nodeTangents(nodes.size());
		std::vector<Vector3> nodeNormals(nodes.size());
		std::vector<Vector3> nodeBitangents(nodes.size());
		std::vector<double> nodeWeights(nodes.size());

		// 1. Compute tangents
		for (size_t i = 0; i < segments.size(); i++) {
			int j = segments[i][0];
			int k = segments[i][1];

			Vector3 segment = nodes[k] - nodes[j];
			double length = norm(segment);

			if (length < 1e-10) {
				continue;
			}

			Vector3 tangent = segment / length;
			nodeTangents[j] += tangent;
			nodeTangents[k] += tangent;
			nodeWeights[j] += length * 0.5;
			nodeWeights[k] += length * 0.5;
		}

		// 2. Normalize tangents
		for (size_t i = 0; i < nodes.size(); i++) {
			if (nodeWeights[i] > 1e-10) {
				nodeTangents[i] = normalize(nodeTangents[i]);
			}
		}

		// 3. Compute consistent normals and bitangents
		for (size_t i = 0; i < nodes.size(); i++) {
			Vector3 tangent = nodeTangents[i];

			if (norm(tangent) < 1e-10) {
				continue;
			}

			if (std::abs(tangent.z) > 0.707) {
				// if tangent is more along z-axis, use x-axis for the cross product
				nodeNormals[i] = normalize(cross(Vector3{ 1, 0, 0 }, tangent));
			}
			else {
				// otherwise, use z-axis for the cross product
				nodeNormals[i] = normalize(cross(Vector3{ 0, 0, 1 }, tangent));
			}

			nodeBitangents[i] = normalize(cross(tangent, nodeNormals[i]));
			nodeNormals[i] = normalize(nodeNormals[i] - dot(nodeNormals[i], tangent) * tangent);
		}

		// 4. Compute medial axis
		std::vector<std::vector<Vector3>> nodeMedialAxis =
			medial_axis(nodes, segments, nodeTangents, nodeNormals, nodeBitangents, maxRadius);

		// 5. Caculate energy terms and gradients for each segment (dirichlet, field, curvature, bilaplacian)
		for (size_t i = 0; i < segments.size(); i++) {
			int j = segments[i][0];
			int k = segments[i][1];

			Vector3 segment = nodes[k] - nodes[j];
			double currentLength = norm(segment);

			if (currentLength < 1e-6) {
				continue;
			}

			Vector3 tangent = segment / currentLength;

			Vector3 normal;
			if (std::abs(tangent.z) > 0.707) {
				// if tangent is more along z-axis, use x-axis for the cross product
				normal = normalize(cross(Vector3{ 1, 0, 0 }, tangent));
			}
			else {
				// otherwise, use z-axis for the cross product
				normal = normalize(cross(Vector3{ 0, 0, 1 }, tangent));
			}

			// ensure normal is perpendicular to tangent
			normal = normalize(normal - dot(normal, tangent) * tangent);

			Vector3 bitangent = cross(tangent, normal);

			double lengthRatio = currentLength / segmentLengths[i];

			// energy terms contribution
			double dirichletTerm = 0;// pow(lengthRatio - 1.0, p); // modified to penalize deviation from rest length
			double fieldAlignedness = 0;// pow(dot(tangent, normal), q);
			double curvatureAlignedness = 0;//pow(dot(tangent, bitangent), q);
			double bilaplacianTerm = pow(lengthRatio, p) * pow(norm(bitangent), q);

			// accumulate energy
			energy += dirichletTerm + fieldAlignedness + curvatureAlignedness + bilaplacianTerm;

			// compute gradient (note: gradient points in the direction of energy increase)
			Vector3 dirichletGrad = Vector3(); alpha* p* pow(lengthRatio - 1.0, p - 1)* (1.0 / segmentLengths[i])* tangent;
			Vector3 fieldAlignednessGrad = Vector3();// alpha* q* pow(dot(tangent, normal), q - 1)* normal;
			Vector3 curvatureAlignednessGrad = Vector3();// alpha* q* pow(dot(tangent, bitangent), q - 1)* bitangent;
			Vector3 bilaplacianGrad = alpha * (
				p * pow(lengthRatio, p - 1) * (1.0 / segmentLengths[i]) * pow(bitangent.norm(), q) * tangent +
				q * pow(lengthRatio, p) * pow(bitangent.norm(), q - 1) * normalize(bitangent)
			);

			// apply gradients to nodes
			gradient[j] += dirichletGrad + fieldAlignednessGrad + curvatureAlignednessGrad + bilaplacianGrad;
			gradient[k] -= dirichletGrad + fieldAlignednessGrad + curvatureAlignednessGrad + bilaplacianGrad;
		}

		// 6. Add medial axis energy term
		for (size_t i = 0; i < nodes.size(); i++) {
			if (nodeMedialAxis[i].size() != 2) {
				continue;
			}

			Vector3 c_0 = nodeMedialAxis[i][0];
			Vector3 c_1 = nodeMedialAxis[i][1];

			double l_0 = norm(nodes[i] - c_0);
			double l_1 = norm(nodes[i] - c_1);

			double medialAxisEnergy = pow(pow(l_0, 2) + pow(l_1, 2), q * 0.5);
			energy += alpha * nodeWeights[i] * medialAxisEnergy / totalCurveLength;

			Vector3 gradL0 = (nodes[i] - c_0) / (l_0 > 1e-10 ? l_0 : 1e-10);
			Vector3 gradL1 = (nodes[i] - c_1) / (l_1 > 1e-10 ? l_1 : 1e-10);

			Vector3 medialAxisGrad = alpha * nodeWeights[i] * q *
				pow(pow(l_0, 2) + pow(l_1, 2), q * 0.5 - 1) *
				(l_0 * gradL0 + l_1 * gradL1) /
				totalCurveLength;

			gradient[i] += medialAxisGrad;
		}

		// 7. Compute descent direction (negative gradient, normalized and scaled)
		double maxGradNorm = 0.0;
		for (size_t i = 0; i < nodes.size(); i++) {
			double gradNorm = gradient[i].norm();
			if (gradNorm > maxGradNorm) {
				maxGradNorm = gradNorm;
			}
		}

		// 8. Normalize descent direction to have a reasonable step size
		double scaleFactor = (maxGradNorm > 1e-10) ? options.stepSizeFactor / maxGradNorm : 0.0;
		for (size_t i = 0; i < nodes.size(); i++) {
			descent[i] = scaleFactor * gradient[i];
		}

		auto end = std::chrono::high_resolution_clock::now();
		std::cout << "Energy calculation time: "
			<< std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
			<< "ms, Energy: " << energy << std::endl;

		return std::make_tuple(
			descent,
			gradient,
			energy,
			nodeMedialAxis
		);
	}
	*/

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
		const VolumeFillingEnergy::Options& options
	) {

		auto start = std::chrono::high_resolution_clock::now();

		double p = options.p;
		double q = options.q;
		double branchRadius = radius * branchRatio;
		double alpha = 4.0 / (branchRadius * branchRadius);

		std::cout << "alpha: " << alpha << ", p: " << p << ", q: " << q << std::endl;

		std::vector<Vector3> descent(nodes.size(), Vector3());
		std::vector<Vector3> gradient(nodes.size(), Vector3());

		double energy = 0.0;
		double totalCurveLength = 0.0;

		for (int i = 0; i < segments.size(); i++) {
			totalCurveLength += segmentLengths[i];
		}

		// 1. Calculate energy terms and gradients for each segment
		for (size_t i = 0; i < segments.size(); i++) {
			int j = segments[i][0];
			int k = segments[i][1];

			Vector3 segment = nodes[k] - nodes[j];
			double currentLength = norm(segment);

			if (currentLength < 1e-6) {
				continue;
			}

			Vector3 tangent = segment / currentLength;

			Vector3 normal;
			if (std::abs(tangent.z) > 0.707) {
				// if tangent is more along z-axis, use x-axis for the cross product
				normal = normalize(cross(Vector3{1, 0, 0}, tangent));
			}
			else {
				// otherwise, use z-axis for the cross product
				normal = normalize(cross(Vector3{ 0, 0, 1 }, tangent));
			}

			// ensure normal is perpendicular to tangent
			normal = normalize(normal - dot(normal, tangent) * tangent);

			Vector3 bitangent = cross(tangent, normal);

			double lengthRatio = currentLength / segmentLengths[i];

			// energy terms contribution
			double dirichletTerm = pow(lengthRatio - 1.0, p); // modified to penalize deviation from rest length
			double fieldAlignedness = pow(dot(tangent, normal), q);
			double curvatureAlignedness = pow(dot(tangent, bitangent), q);
			double bilaplacianTerm = pow(lengthRatio, p) * pow(norm(bitangent), q);

			// accumulate energy
			energy += dirichletTerm + fieldAlignedness + curvatureAlignedness + bilaplacianTerm;

			// compute gradient (note: gradient points in the direction of energy increase)
			Vector3 dirichletGrad = alpha * p * pow(lengthRatio - 1.0, p - 1) * (1.0 / segmentLengths[i]) * tangent;
			Vector3 fieldAlignednessGrad = alpha * q * pow(dot(tangent, normal), q - 1) * normal;
			Vector3 curvatureAlignednessGrad = alpha * q * pow(dot(tangent, bitangent), q - 1) * bitangent;
			Vector3 bilaplacianGrad = alpha * (
				p * pow(lengthRatio, p - 1) * (1.0 / segmentLengths[i]) * pow(bitangent.norm(), q) * tangent +
				q * pow(lengthRatio, p) * pow(bitangent.norm(), q - 1) * normalize(bitangent)
			);

			// apply gradients to nodes
			gradient[j] += dirichletGrad + fieldAlignednessGrad + curvatureAlignednessGrad + bilaplacianGrad;
			gradient[k] -= dirichletGrad + fieldAlignednessGrad + curvatureAlignednessGrad + bilaplacianGrad;
		}


		// 2. Compute descent direction (negative gradient, normalized and scaled)
		double maxGradNorm = 0.0;
		for (size_t i = 0; i < nodes.size(); i++) {
			double gradNorm = gradient[i].norm();
			if (gradNorm > maxGradNorm) {
				maxGradNorm = gradNorm;
			}
		}
			
		const double stepSizeFactor = 0.01; // scaling factor for step size
		double scaleFactor = (maxGradNorm > 1e-10) ? stepSizeFactor / maxGradNorm : 0.0;

		for (size_t i = 0; i < nodes.size(); i++) {
			descent[i] = scaleFactor * gradient[i];
		}

		auto end = std::chrono::high_resolution_clock::now();
		std::cout << "Energy calculation time: "
			<< std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
			<< "ms, Energy: " << energy << std::endl;

		return std::make_tuple(
			descent,
			gradient,
			energy,
			std::vector<std::vector<Vector3>>{}
		);
	}
}