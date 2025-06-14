#include <openvdb/openvdb.h>
#include <openvdb/tools/MeshToVolume.h>

#include <geometrycentral/surface/surface_mesh.h>
#include <geometrycentral/surface/meshio.h>

using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace modules {
	std::vector<Vector3> mesh_to_nodes(std::string filename) {
		std::cout << "Initializing OpenVDB" << filename << std::endl;
		openvdb::initialize();

		std::unique_ptr<SurfaceMesh> mesh;
		std::unique_ptr<VertexPositionGeometry> geometry;
		std::tie(mesh, geometry) = readSurfaceMesh(filename);

		if (!mesh || !geometry) {
			std::cerr << "Error reading mesh from file: " << filename << std::endl;
			std::abort();
		}

		std::vector<Vector3> nodes;
		std::cout << "Loading vertices" << std::endl;
		for (Vertex v : mesh->vertices()) {
			auto p = geometry->inputVertexPositions[v];
			nodes.emplace_back(Vector3{ p.x, p.y, p.z });
		}

		std::cout << "Loaded " << nodes.size() << " nodes from mesh." << std::endl;
		return nodes;
	}

	openvdb::FloatGrid::Ptr mesh_to_sdf(std::string filename, double voxelsize) {
		std::cout << "Initializing OpenVDB" << std::endl;
		openvdb::initialize();

		std::unique_ptr<SurfaceMesh> mesh;
		std::unique_ptr<VertexPositionGeometry> geometry;
		std::tie(mesh, geometry) = readSurfaceMesh(filename);

		if (!mesh || !geometry) {
			std::cerr << "Error reading mesh from file: " << filename << std::endl;
			std::abort();
		}


		std::vector<openvdb::Vec3s> points;
        std::vector<openvdb::Vec3I> triangles;

        // Convert vertices
		std::cout << "Loading vertices" << std::endl;
        for (Vertex v : mesh->vertices()) {
            auto p = geometry->inputVertexPositions[v];
            points.emplace_back(p.x, p.y, p.z);

		}

        // Convert faces (only triangles)
		std::cout << "Loading triangles" << std::endl;
		for (Face f : mesh->faces()) {
			if (f.degree() != 3) {
				std::cerr << "Mesh contains non-triangular faces, which are not supported." << std::endl;
				std::abort();
			}

			// Collect vertex indices by iterating through adjacent vertices
			std::vector<size_t> vertexIndices;
			for (Vertex v : f.adjacentVertices()) {
				vertexIndices.push_back(v.getIndex());
			}

			triangles.emplace_back(
				openvdb::Vec3I(vertexIndices[0], vertexIndices[1], vertexIndices[2]));
		}

        // Generate SDF
		std::cout << "Generating SDF" << std::endl;
        openvdb::math::Transform::Ptr transform =
            openvdb::math::Transform::createLinearTransform(voxelsize);

        openvdb::FloatGrid::Ptr sdfGrid =
            openvdb::tools::meshToLevelSet<openvdb::FloatGrid>(*transform, points, triangles);

		std::cout << "SDF grid created with " << sdfGrid->activeVoxelCount() << " active voxels." << std::endl;
		return sdfGrid;
	}
}