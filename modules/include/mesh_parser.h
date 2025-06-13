#include <openvdb/openvdb.h>

namespace modules {
	std::vector<Vector3> mesh_to_nodes(std::string filename);
	openvdb::FloatGrid::Ptr mesh_to_sdf(std::string filename, double voxelsize);
}