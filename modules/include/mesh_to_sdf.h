#include <openvdb/openvdb.h>

namespace modules {
	openvdb::FloatGrid::Ptr mesh_to_sdf(std::string filename, double voxelsize);
}