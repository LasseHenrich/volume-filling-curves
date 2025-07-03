#include "volume_filling_energy.h"
#include <chrono>
#include <iostream>
#include <glm/glm.hpp>
#include <Eigen/Dense>
#include <array>
#include <medial_axis.h>
#include <primitive_energies.h>
#include <TinyAD/Utils/LinearSolver.hh>
#include <TinyAD/Utils/NewtonDirection.hh>
#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/GridOperators.h>

namespace modules {
	template<typename T>
	T min(T a, T b) {
		return (a < b) ? a : b;
	}

	const double branchRatio = std::sqrt(std::sqrt(2));

	template<typename Element>
	auto calculate_medial_axis_sdf_energy(
		Element& element,
		const std::vector<std::vector<Vector3>>& nodeMedialAxis,
		const scene_file::SceneObject& options,
		openvdb::tools::GridSampler<openvdb::FloatGrid, openvdb::tools::BoxSampler>& sampler,
		double radius,
		double maxRadius,
		double alpha,
		const std::vector<double>& nodeWeight,
		double q,
		double totalCurveLength) -> TINYAD_SCALAR_TYPE(element)
	{
		using T = TINYAD_SCALAR_TYPE(element);
		int nodeId = element.handle;
		Eigen::Vector3<T> pos = element.variables(nodeId);

		auto _c_0 = nodeMedialAxis[nodeId][0];
		auto _c_1 = nodeMedialAxis[nodeId][1];

		Eigen::Vector3d c_0(_c_0.x, _c_0.y, _c_0.z);
		Eigen::Vector3d c_1(_c_1.x, _c_1.y, _c_1.z);

		T l_0 = (pos - c_0).norm();
		T l_1 = (pos - c_1).norm();


		if (options.volume.volumeType == scene_file::VolumeType::MESH && options.volume.convert_to_sdf && options.use_volumetric_energy) {
			// Get current passive position for sampling
			openvdb::Vec3d world_pos(TinyAD::to_passive(pos(0)),
				TinyAD::to_passive(pos(1)),
				TinyAD::to_passive(pos(2)));

			// Sample SDF value at current position
			double sdf_value = sampler.wsSample(world_pos);

			// Compute gradient using finite differences
			const double h = 1e-1;
			double sdf_px = sampler.wsSample(world_pos + openvdb::Vec3d(h, 0, 0));
			double sdf_mx = sampler.wsSample(world_pos - openvdb::Vec3d(h, 0, 0));
			double sdf_py = sampler.wsSample(world_pos + openvdb::Vec3d(0, h, 0));
			double sdf_my = sampler.wsSample(world_pos - openvdb::Vec3d(0, h, 0));
			double sdf_pz = sampler.wsSample(world_pos + openvdb::Vec3d(0, 0, h));
			double sdf_mz = sampler.wsSample(world_pos - openvdb::Vec3d(0, 0, h));

			double grad_x = (sdf_px - sdf_mx) / (2.0 * h);
			double grad_y = (sdf_py - sdf_my) / (2.0 * h);
			double grad_z = (sdf_pz - sdf_mz) / (2.0 * h);

			// Create differentiable SDF using linear approximation
			/*T sdf_approx_at_node = T(sdf_value) +
				T(grad_x) * (pos(0) - T(world_pos[0])) +
				T(grad_y) * (pos(1) - T(world_pos[1])) +
				T(grad_z) * (pos(2) - T(world_pos[2]));*/

			// we don't just wsSample at c_0 and c_1, as this would be insensitive
			// to "jumping" across very tight surface regions.
			// Instead, we use the local gradient to approximate the SDF at c_0 and c_1,
			// assuming a linear gradient field around the current position.
			// (Also, TinyAD needs a notion of a gradient, which is why would need to
			// compute it anyway.)
			T sdf_approx_at_c_0 = T(sdf_value) +
				T(grad_x) * (c_0(0) - T(world_pos[0])) +
				T(grad_y) * (c_0(1) - T(world_pos[1])) +
				T(grad_z) * (c_0(2) - T(world_pos[2]));

			T sdf_approx_at_c_1 = T(sdf_value) +
				T(grad_x) * (c_1(0) - T(world_pos[0])) +
				T(grad_y) * (c_1(1) - T(world_pos[1])) +
				T(grad_z) * (c_1(2) - T(world_pos[2]));

			// backproject... note that I'm not sure whether this is the
			// correct way to backproject
			T overshoot_0 = sdf_approx_at_c_0 - radius;
			if (overshoot_0 > 0) {
				l_0 = min(l_0 - overshoot_0, (T)maxRadius);
			}
			T overshoot_1 = sdf_approx_at_c_1 - radius;
			if (overshoot_1 > 0) {
				l_1 = min(l_1 - overshoot_1, (T)maxRadius);
			}
		}

		auto result = alpha * nodeWeight[nodeId] * (
			pow(
				pow(l_0, 2) +
				pow(l_1, 2),
			q / 2)
			) / totalCurveLength;

		return result;
	}

	std::tuple<
		std::vector<Vector3>, // descent direction
		std::vector<Vector3>, // gradient direction
		double, // energy
		std::vector<std::vector<Vector3>> // medial axis
	> volume_filling_energy_curve(
		const Curve& curve,
		const scene_file::SceneObject& options
	) {
		std::vector<Vector3> nodes = curve.nodes;
		std::vector<std::array<int, 2>> segments = curve.segments;
		std::vector<double> segmentLengths = curve.segmentLengths;

		// making sure that we have at least three segments
		if (segments.size() < 3) {
			std::cerr << "Error: At least three segments are required." << std::endl;
			std::abort();
		}

		auto start = std::chrono::high_resolution_clock::now();

		double radius = options.radius;
		double maxRadius = options.rmax;
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
		
		// 1. Compute tangents
		for (size_t i = 0; i < segments.size(); i++) {
			int j = segments[i][0];
			int k = segments[i][1];

			Vector3 t = normalize(nodes[k] - nodes[j]);
			segmentTangents[i] = t;
		}

		for (size_t i = 0; i < segments.size(); i++) {
			Vector3 t = segmentTangents[i];
			Vector3 b;
			Vector3 n;

			if (true) { // 3D
				Vector3 t_prev, t_next;
				if (i == 0) {
					t_prev = segmentTangents[segments.size() - 1];
					t_next = segmentTangents[i + 1];
				}
				else if (i == segments.size() - 1) {
					t_prev = segmentTangents[i - 1];
					t_next = segmentTangents[0];
				}
				else {
					t_prev = segmentTangents[i - 1];
					t_next = segmentTangents[i + 1];
				}
				n = normalize(cross(t_prev, t_next));
				if (norm(t_prev - t_next) < 1e-5) {  // curvature=NaN iff t_prev \approx t_next
					n = arbitrary_normal(t);
				}
				b = normalize(cross(t, n));
			}
			else { // 2D
				n = Vector3{ 0, 1, 0 }; // DEBUG for 2D case. ToDo
				b = cross(n, t);
			}

			segmentNormal[i] = n;
			segmentBitangent[i] = b;

			Eigen::Matrix3d R;
			R << t.x, b.x, n.x,
				t.y, b.y, n.y,
				t.z, b.z, n.z;

			rotationMatrix[i] = R;
		}

		// 2. Dirichlet term (Penalize length)
		if (options.use_length_energy) {
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

				// R.transpose() multiplication: Differenz wird auch Richtung von Kante projiziert
				// andere zwei Richtungen haben keinen groﬂen Einfluss
				Eigen::Vector3<T> p0 = R.transpose() * (x0 - v0);
				Eigen::Vector3<T> p1 = R.transpose() * (x1 - v1);

				T dx = abs(p0(0) - (p1(0) + edgeLen));
				T dy = abs(p0(1) - p1(1));
				T dz = abs(p0(2) - p1(2));

				auto result = (
					pow(pow(dx, 2) + pow(dy, 2) + pow(dz, 2), p / 2)
					) / (edgeLen * totalCurveLength);

				return result;
			});
		}

		auto dirichletEnd = std::chrono::high_resolution_clock::now();

		// 3. Medial axis term
		auto tetTopologyEnd = std::chrono::high_resolution_clock::now();

		// 3.1 Node-based tangents
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

		// 3.2 Medial axis
		auto nodeMedialAxis = modules::medial_axis_curve(
			nodes,
			segments,
			nodeTangents,
			nodeNormals,
			nodeBitangents,
			maxRadius,
			options.volume.mesh_points
		);

		auto medialAxisEnd = std::chrono::high_resolution_clock::now();

		// 3.3 Add medial axis energy term AS WELL AS volumetric energy term when using SDF
		// note that with our current setup we just have one uniform alpha
		double totalMedialAxisEnergy = 0;
		openvdb::tools::GridSampler<openvdb::FloatGrid, openvdb::tools::BoxSampler> sampler(*options.volume.sdf);
		func.add_elements<1>(TinyAD::range(nodes.size()), [&](auto& element) -> TINYAD_SCALAR_TYPE(element) {
			return calculate_medial_axis_sdf_energy(
				element,
				nodeMedialAxis,
				options,
				sampler,
				radius,
				maxRadius,
				alpha,
				nodeWeight,
				q,
				totalCurveLength
			);
		});

		std::cout << "Total Medial Axis Energy: " << totalMedialAxisEnergy << std::endl;

		auto repulsiveEnd = std::chrono::high_resolution_clock::now();

		// 4. Volume-constraint terms
		// commented-out as they seem redundant with the backprojection

		auto volume = options.volume;
		if (volume.volumeType == scene_file::VolumeType::PRIMITIVE) {
			if (volume.primitive_type == scene_file::PrimitiveType::SPHERE) {
				if (volume.primitive_params.size() < 1) {
					std::cerr << "Error: Sphere primitive requires 1 parameter (radius)" << std::endl;
					std::abort();
				}
				double sphereRadius = volume.primitive_params[0];
				func.add_elements<1>(TinyAD::range(nodes.size()), [&](auto& element)->TINYAD_SCALAR_TYPE(element) {
					using T = TINYAD_SCALAR_TYPE(element);
					int nodeId = element.handle;
					Eigen::Vector3<T> pos = element.variables(nodeId);

					// Signed distance to sphere
					T sdf = pos.norm() - sphereRadius;

					return (sdf > 0.0 ? 1000.0 * pow(sdf, 2) : 0) / totalCurveLength;
				});
			}
			else if (volume.primitive_type == scene_file::PrimitiveType::BOX) {
				if (volume.primitive_params.size() < 3) {
					std::cerr << "Error: Box primitive requires 3 parameters (half-extents)" << std::endl;
					std::abort();
				}
				Eigen::Vector3d boxExtents(
					volume.primitive_params[0],
					volume.primitive_params[1],
					volume.primitive_params[2]
				);

				func.add_elements<1>(TinyAD::range(nodes.size()), [&](auto& element)->TINYAD_SCALAR_TYPE(element) {
					using T = TINYAD_SCALAR_TYPE(element);
					int nodeId = element.handle;
					Eigen::Vector3<T> pos = element.variables(nodeId);

					// Signed distance to box
					Eigen::Vector3<T> q = pos.cwiseAbs() - boxExtents;
					T exterior_dist = (q.cwiseMax(0.0)).norm();
					T interior_dist = fmin(fmax(q(0), fmax(q(1), q(2))), T(0.0));
					T sdf = exterior_dist + interior_dist;

					return (sdf > 0.0 ? 1000.0 * pow(sdf, 2) : 0) / totalCurveLength;
				});
			}
			else if (volume.primitive_type == scene_file::PrimitiveType::ROUNDBOX) {
				if (volume.primitive_params.size() < 4) {
					std::cerr << "Error: Roundbox primitive requires 4 parameters (half-extents and radius)" << std::endl;
					std::abort();
				}
				Eigen::Vector3d boxExtents(
					volume.primitive_params[0],
					volume.primitive_params[1],
					volume.primitive_params[2]
				);
				double radius = volume.primitive_params[3];
				func.add_elements<1>(TinyAD::range(nodes.size()), [&](auto& element)->TINYAD_SCALAR_TYPE(element) {
					using T = TINYAD_SCALAR_TYPE(element);
					int nodeId = element.handle;
					Eigen::Vector3<T> pos = element.variables(nodeId);

					// Signed distance to roundbox
					Eigen::Vector3<T> q = pos.cwiseAbs() - boxExtents;
					T exterior_dist = (q.cwiseMax(0.0)).norm();
					T interior_dist = fmin(fmax(q(0), fmax(q(1), q(2))), T(0.0));
					T sdf = exterior_dist + interior_dist;
					T roundbox_sdf = sdf - radius;

					return (roundbox_sdf > 0.0 ? 1000.0 * pow(roundbox_sdf, 2) : 0) / totalCurveLength;
				});
			}
			else if (volume.primitive_type == scene_file::PrimitiveType::TORUS) {
				if (volume.primitive_params.size() < 2) {
					std::cerr << "Error: Torus primitive requires 2 parameters (major and minor radii)" << std::endl;
					std::abort();
				}
				double R = volume.primitive_params[0]; // Major radius
				double r = volume.primitive_params[1]; // Minor radius
				func.add_elements<1>(TinyAD::range(nodes.size()), [&](auto& element)->TINYAD_SCALAR_TYPE(element) {
					using T = TINYAD_SCALAR_TYPE(element);
					int nodeId = element.handle;
					Eigen::Vector3<T> pos = element.variables(nodeId);

					// Signed distance to torus
					Eigen::Vector2<T> q = Eigen::Vector2<T>(Eigen::Vector2<T>(pos(0), pos(1)).norm() - R, pos(2));
					T torus_sdf = q.norm() - r;

					return (torus_sdf > 0.0 ? 1000.0 * pow(torus_sdf, 2) : 0) / totalCurveLength;
				});
			}
			else {
				std::abort(); // Unsupported primitive type
			}
		}
		else if (volume.volumeType == scene_file::VolumeType::MESH && volume.convert_to_sdf) {
			// Case handled in medial axis energy term
		}
		else if (volume.volumeType == scene_file::VolumeType::MESH) {
			// To nothing, handled via medial_axis calculation
		}
		else if (volume.volumeType == scene_file::VolumeType::SDF) {
			// Handle sdf volume type if needed
			std::cerr << "SDF volume type not implemented yet." << std::endl;
		}
		else {
			std::cerr << "Unknown volume type or no volume" << std::endl;
		}


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

		return std::make_tuple(
			descent,
			gradient,
			f,
			nodeMedialAxis
		);
	}

	std::tuple<
		std::vector<Vector3>, // descent direction
		std::vector<Vector3>, // gradient direction
		double,              // energy
		std::vector<std::vector<Vector3>> // medial axis (if applicable)
	> volume_filling_energy_surface(
		const Surface& surface,
		const scene_file::SceneObject& options
	) {
		geometrycentral::surface::ManifoldSurfaceMesh* mesh = surface.mesh.get();
		geometrycentral::surface::VertexPositionGeometry* geometry = surface.geometry.get();

		size_t numNodes = mesh->nVertices();
		size_t numFaces = mesh->nFaces();

		std::vector<Vector3> nodes(numNodes);
		for (size_t i = 0; i < numNodes; i++) {
			nodes[i] = geometry->vertexPositions[i];
		}

		// making sure that we have at least four faces
		if (numFaces < 4) {
			std::cerr << "Error: At least four faces are required." << std::endl;
			std::abort();
		}

		auto start = std::chrono::high_resolution_clock::now();

		double radius = options.radius;
		double maxRadius = options.rmax;
		double p = options.p;
		double q = options.q;
		double branchRadius = radius * branchRatio;
		double alpha = 4 / (pow(branchRadius, 2));

		std::cout << "alpha: " << alpha << ", p: " << p << ", q: " << q << std::endl;

		std::vector<Vector3> descent(nodes.size(), Vector3{ .0, .0, .0 }), gradient(nodes.size(), Vector3{ .0, .0, .0 });

		double totalSurfaceArea = 0.0;
		for (const auto& face : mesh->faces()) {
			totalSurfaceArea += geometry->faceArea(face);
		}

		auto func = TinyAD::scalar_function<3>(TinyAD::range(nodes.size()));

		if (options.use_length_energy) {
			// ToDo
		}

		auto dirichletEnd = std::chrono::high_resolution_clock::now();

		// 3. Medial axis term
		auto tetTopologyEnd = std::chrono::high_resolution_clock::now();

		// 3.1 Node-based normals
		std::vector<Vector3> nodeNormals(nodes.size());
		geometry->requireVertexNormals(); // computes the vertexNormals array
		VertexData<Vector3>& vertexNormals = geometry->vertexNormals;

		 for (size_t i = 0; i < numNodes; i++) {
			 Vertex v = mesh->vertex(i);
			 nodeNormals[i] = vertexNormals[v];
		 }

		auto tangentEvalEnd = std::chrono::high_resolution_clock::now();

		// 3.2 Medial axis
		auto nodeMedialAxis = modules::medial_axis_surface(
			nodes,
			nodeNormals,
			maxRadius
		);

		auto medialAxisEnd = std::chrono::high_resolution_clock::now();

		// 3.3 Add medial axis energy term AS WELL AS volumetric energy term when using SDF
		// note that with our current setup we just have one uniform alpha
		openvdb::tools::GridSampler<openvdb::FloatGrid, openvdb::tools::BoxSampler> sampler(*options.volume.sdf);
		func.add_elements<1>(TinyAD::range(nodes.size()), [&](auto& element)->TINYAD_SCALAR_TYPE(element) {
			return calculate_medial_axis_sdf_energy(
				element,
				nodeMedialAxis,
				options,
				sampler,
				radius,
				maxRadius,
				alpha,
				std::vector<double>(nodes.size(), 1.0), // uniform weight for all nodes. ToDo
				q,
				totalSurfaceArea
			);
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
			descent[idx] = Vector3{ p(0), p(1), p(2) };
			});

		func.x_to_data(g, [&](int idx, const Eigen::Vector3d& p) {
			gradient[idx] = Vector3{ p(0), p(1), p(2) };
			});

		auto end = std::chrono::high_resolution_clock::now();

		return std::make_tuple(
			descent,
			gradient,
			f,
			std::vector<std::vector<Vector3>>{}
		);
	}


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
		const scene_file::SceneObject& options
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

	Vector3 arbitrary_normal(const Vector3& t) {
		if (std::abs(t.x) <= std::abs(t.y) && std::abs(t.x) <= std::abs(t.z)) return normalize(cross(t, Vector3{ 1,0,0 }));
		if (std::abs(t.y) <= std::abs(t.z)) return normalize(cross(t, Vector3{ 0,1,0 }));
		return normalize(cross(t, Vector3{ 0,0,1 }));
	}
}