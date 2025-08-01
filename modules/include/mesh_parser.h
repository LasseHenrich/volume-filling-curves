#include <openvdb/openvdb.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/surface_mesh.h>
#include <structures.h>

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

	Surface file_to_geometrycentral_data(std::string filename);
	PolyscopeMeshData geometrycentral_to_polyscope_data(Surface const* geometrycentralMeshData);
	OpenVDBMeshData polyscope_to_openvdb_data(PolyscopeMeshData const* polyscopeMeshData);
	openvdb::FloatGrid::Ptr openvdb_mesh_to_sdf(OpenVDBMeshData const* openvdbMeshData, double voxelsize, float halfwidth);
	bool sdf_is_watertight(openvdb::FloatGrid::Ptr sdfGrid);
	bool mesh_is_manifold(std::unique_ptr<ManifoldSurfaceMesh>* mesh);
}