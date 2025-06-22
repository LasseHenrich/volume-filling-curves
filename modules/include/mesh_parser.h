#include <openvdb/openvdb.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/surface_mesh.h>

using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace modules {
	struct PolyscopeMeshData {
		std::vector<Vector3> vertices;
		std::vector<std::array<size_t, 3>> faces;
	};
	struct OpenVDBMeshData {
		std::vector<openvdb::Vec3s> vertices;
		std::vector<openvdb::Vec3I> faces;
	};
	struct OpenVDBMeshAndGeometry {
		std::unique_ptr<SurfaceMesh> mesh;
		std::unique_ptr<VertexPositionGeometry> geometry;
	};

	PolyscopeMeshData file_to_polyscope_data(std::string filename);
	OpenVDBMeshData polyscope_to_openvdb_data(PolyscopeMeshData polyscopeMesh);
	openvdb::FloatGrid::Ptr openvdb_mesh_to_sdf(OpenVDBMeshData openvdbMesh, double voxelsize, float halfwidth);
	bool sdf_is_watertight(openvdb::FloatGrid::Ptr sdfGrid);
	bool mesh_is_watertight(std::unique_ptr<SurfaceMesh>* mesh);
}