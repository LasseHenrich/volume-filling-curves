#include "mesh_parser.h"
#include <openvdb/openvdb.h>
#include <openvdb/tools/MeshToVolume.h>

namespace modules {
	PolyscopeMeshData file_to_polyscope_data(std::string filename) {
		std::unique_ptr<SurfaceMesh> mesh;
		std::unique_ptr<VertexPositionGeometry> geometry;
		std::tie(mesh, geometry) = readSurfaceMesh(filename);

		if (!mesh || !geometry) {
			std::cerr << "Error reading mesh from file: " << filename << std::endl;
			std::abort();
		}

		std::vector<Vector3> vertices;
		std::cout << "Loading vertices" << std::endl;
		for (Vertex v : mesh->vertices()) {
			auto p = geometry->inputVertexPositions[v];
			vertices.emplace_back(Vector3{ p.x, p.y, p.z });
		}
		std::cout << "Loaded " << vertices.size() << " vertices from mesh." << std::endl;

		std::vector<std::array<size_t, 3>> faces;
		std::cout << "Loading triangles" << std::endl;
		for (Face f : mesh->faces()) {
			if (f.degree() != 3) {
				std::cerr << "Mesh contains non-triangular faces, which are not supported." << std::endl;
				std::abort();
			}
			std::vector<size_t> vertexIndices;
			for (Vertex v : f.adjacentVertices()) {
				vertexIndices.push_back(v.getIndex());
			}
			faces.emplace_back(
				std::array<size_t, 3>{vertexIndices[0], vertexIndices[1], vertexIndices[2]});
		}
		std::cout << "Loaded " << faces.size() << " triangles from mesh." << std::endl;

		PolyscopeMeshData meshData;
		meshData.vertices = vertices;
		meshData.faces = faces;
		return meshData;
	}

	OpenVDBMeshData polyscope_to_openvdb_data(PolyscopeMeshData polyscopeMesh) {
		OpenVDBMeshData openvdbMesh;
		openvdbMesh.vertices.reserve(polyscopeMesh.vertices.size());
		openvdbMesh.faces.reserve(polyscopeMesh.faces.size());
		for (const auto& vertex : polyscopeMesh.vertices) {
			openvdbMesh.vertices.emplace_back(openvdb::Vec3s(vertex.x, vertex.y, vertex.z));
		}
		for (const auto& face : polyscopeMesh.faces) {
			openvdbMesh.faces.emplace_back(
				openvdb::Vec3I(face[0], face[1], face[2]));
		}
		return openvdbMesh;
	}

	openvdb::FloatGrid::Ptr openvdb_mesh_to_sdf(OpenVDBMeshData openvdbMesh, double voxelsize) {
		auto points = openvdbMesh.vertices;
		auto triangles = openvdbMesh.faces;

        openvdb::math::Transform::Ptr transform =
            openvdb::math::Transform::createLinearTransform(voxelsize);

        openvdb::FloatGrid::Ptr sdfGrid =
            openvdb::tools::meshToLevelSet<openvdb::FloatGrid>(*transform, points, triangles);

		return sdfGrid;
	}

	bool sdf_is_watertight(openvdb::FloatGrid::Ptr sdfGrid) {
		// count negative values in SDF
		int negativeCount = 0;
		for (auto iter = sdfGrid->cbeginValueOn(); iter; ++iter) {
			if (iter.getValue() < 0) {
				negativeCount++;
			}
		}
		// check if there are significantly less negative values than positive values
		int positiveCount = sdfGrid->activeVoxelCount() - negativeCount;
		if (negativeCount < positiveCount / 10) {
			std::cerr << "SDF grid has significantly fewer negative values than positive values and is probably not watertight" << std::endl;
			return false;
		}
		
		std::cout << "SDF grid has " << negativeCount << " negative values and "
			<< positiveCount << " positive values." << std::endl;
		return true;
	}

	// doesn't work correctly right now
	bool mesh_is_watertight(std::unique_ptr<SurfaceMesh>* mesh) {
		for (Edge e : (*mesh)->edges()) {
			if (e.isBoundary()) {
				std::cout << "Mesh is not watertight: found boundary edge." << std::endl;
				return false;
			}
		}

		for (Edge e : (*mesh)->edges()) {
			int faceCount = 0;
			for (Face f : e.adjacentFaces()) {
				faceCount++;
			}
			if (faceCount != 2) {
				std::cout << "Mesh is not watertight: edge " << e.getIndex() << " has " << faceCount << " adjacent faces." << std::endl;
				return false;
			}
		}

		for (Vertex v : (*mesh)->vertices()) {
			if (v.degree() == 0) {
				std::cout << "Mesh is not watertight: vertex " << v.getIndex() << " has degree 0." << std::endl;
				return false;
			}
		}

		std::cout << "Mesh is watertight." << std::endl;
		return true;
	}
}