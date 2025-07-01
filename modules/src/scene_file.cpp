#include "scene_file.h"
#include "mesh_parser.h"
#include <regex>

namespace modules {
	using namespace std;

    void splitString(const std::string& str, std::vector<string>& cont, char delim) {
        std::stringstream ss(str);
        std::string token;
        while (std::getline(ss, token, delim)) {
            cont.push_back(token);
        }
    }

    std::string getDirectoryFromPath(std::string str) {
        using namespace std;
        vector<string> parts;
        splitString(str, parts, '/');

        int nParts = parts.size();
        if (nParts == 1) return "./";

        string path = "";

        for (int i = 0; i < nParts - 1; i++) {
            path = path + parts[i] + "/";
        }

        return path;
    }

    void processLine(scene_file::SceneObject& scene, std::string directory, std::vector<std::string>& parts) {
        string key = parts[0];

        if (key == "#" || key == "//") {
            return;
        }

        if (key == "filling_manifold") {
            scene.fillingManifoldFileName = directory + parts[1];
        }
        else if (key == "radius") {
            scene.radius = stod(parts[1]);
            scene.h = scene.radius * igl::PI / 20; // Define h after radius !!!
            scene.rmax = scene.radius * 5; // Define rmax after radius !!!
        }
        else if (key == "h") {
            scene.h = stod(parts[1]);
        }
        else if (key == "p") {
            scene.p = stod(parts[1]);
        }
        else if (key == "q") {
            scene.q = stod(parts[1]);
        }
        else if (key == "rmax") {
            scene.rmax = stod(parts[1]);
        }
        else if (key == "use_backprojection") {
            scene.use_backprojection = (parts.size() < 2 || parts[1] == "true");
        }
		else if (key == "filling_dimension") {
			if (parts.size() < 2) {
				cerr << "Error: filling dimension not specified" << endl;
				exit(1);
			}
			scene.filling_dimension = stoi(parts[1]);
		}

        else if (key == "volume") {
            if (parts.size() < 2) {
                cerr << "Error: surface type not specified" << endl;
                exit(1);
            }
            scene.volume = scene_file::SceneObject_Volume();
            if (parts[1] == "primitive") {
				scene.volume.volumeType = scene_file::VolumeType::PRIMITIVE;
				if (parts.size() < 3) {
					cerr << "Error: primitive type not specified" << endl;
					exit(1);
				}
				if (parts[2] == "sphere") {
					scene.volume.primitive_type = scene_file::PrimitiveType::SPHERE;
					if (parts.size() != 4) {
						cerr << "Error: sphere radius not specified" << endl;
						exit(1);
					}
					scene.volume.primitive_params.push_back(stod(parts[3]));
				}
				else if (parts[2] == "box") {
					scene.volume.primitive_type = scene_file::PrimitiveType::BOX;
					if (parts.size() != 6) {
						cerr << "Error: box extents not specified" << endl;
						exit(1);
					}
					for (size_t i = 3; i < 6; i++) {
						scene.volume.primitive_params.push_back(stod(parts[i]));
					}
				}
				else if (parts[2] == "roundbox") {
					scene.volume.primitive_type = scene_file::PrimitiveType::ROUNDBOX;
					if (parts.size() != 7) {
						cerr << "Error: roundbox extents or radius not specified" << endl;
						exit(1);
					}
					for (size_t i = 3; i < 7; i++) {
						scene.volume.primitive_params.push_back(stod(parts[i]));
					}
				}
				else if (parts[2] == "torus") {
					scene.volume.primitive_type = scene_file::PrimitiveType::TORUS;
					if (parts.size() != 5) {
						cerr << "Error: torus radii not specified" << endl;
						exit(1);
					}
					for (size_t i = 3; i < 5; i++) {
						scene.volume.primitive_params.push_back(stod(parts[i]));
					}
				}
				else {
					cerr << "Error: unknown primitive type '" << parts[2] << "'" << endl;
					exit(1);
				}
			}
			else if (parts[1] == "sdf") {
				scene.volume.volumeType = scene_file::VolumeType::SDF;
			}
			else if (parts[1] == "mesh") {
				scene.volume.volumeType = scene_file::VolumeType::MESH;
				if (parts.size() < 3) {
					cerr << "Error: mesh filename not specified" << endl;
					exit(1);
				}
				scene.volume.mesh_filename = directory + parts[2];

                if (parts.size() < 4) {
					scene.volume.convert_to_sdf = false; // Default to true
				}
                else {
                    scene.volume.convert_to_sdf = (parts[3] == "true");
                }

                if (parts.size() > 4 && scene.volume.convert_to_sdf) {
                    scene.volume.mesh_to_sdf_voxelsize = stod(parts[4]);
                }

				if (parts.size() > 5 && scene.volume.convert_to_sdf) {
					scene.volume.mesh_to_sdf_halfwidth = stod(parts[5]);
				}
			}
			else {
				cerr << "Error: unknown surface type '" << parts[1] << "'" << endl;
				exit(1);
            }
        }
        else if (key == "visualize_volume") {
			scene.visualizeVolume = parts.size() < 2 || (parts[1] == "true");
		}

        else if (key == "timestep") {
            scene.timestep = stod(parts[1]);
        }
        else if (key == "field_aligned") {
            scene.w_fieldAlignedness = stod(parts[1]);
        }
        else if (key == "curxvature_aligned") {
            scene.w_curvatureAlignedness = stod(parts[1]);
        }
        else if (key == "bilaplacian") {
            scene.w_bilaplacian = stod(parts[1]);
        }
        else if (key == "varying_alpha") {
            scene.varyingAlpha = true;
        }
        else if (key == "geodesic_medial_axis") {
            scene.useGeodesicMedialAxis = true;
        }
        else if (key == "excecute_only") {
            scene.excecuteOnly = true;
        }
    }

    scene_file::SceneObject read_scene(std::string filename) {
        string directory = getDirectoryFromPath(filename);

        ifstream inFile;
        inFile.open(filename);

        if (!inFile) {
            cerr << "Could not open file " << filename << endl;
            exit(1);
        }

        scene_file::SceneObject scene;

        std::vector<std::string> parts;
        for (std::string line; std::getline(inFile, line); ) {
            if (line == "" || line == "\n") continue;
            parts.clear();
            splitString(line, parts, ' ');
            processLine(scene, directory, parts);
        }

        inFile.close();
        return scene;
    }

    std::tuple<
        std::vector<Vector3>, // nodes
        std::vector<std::array<int, 2>> // segments
    > read_curve(std::string filename) {
        std::vector<Vector3> nodes;
        std::vector<std::array<int, 2>> segments;
        std::ifstream file(filename);

        if (!file) {
            cerr << "Could not open file " << filename << endl;
            exit(1);
        }

        std::string line;
        // Regular expression to match the format <x, y, z>
        std::regex coordPattern("<([^,]+),\\s*([^,]+),\\s*([^>]+)>");
        // Regular expression to match the format "segment: a, b"
        std::regex segmentPattern("segment:\\s*(\\d+),\\s*(\\d+)");
        std::smatch matches;

        while (std::getline(file, line)) {
            // Check if line contains node coordinates
            if (std::regex_search(line, matches, coordPattern) && matches.size() == 4) {
                try {
                    float x = std::stof(matches[1].str());
                    float y = std::stof(matches[2].str());
                    float z = std::stof(matches[3].str());

                    nodes.push_back(Vector3{ x, y, z });
                }
                catch (const std::exception& e) {
                    std::cerr << "Error parsing coordinates: " << line << std::endl;
                    std::cerr << "Exception: " << e.what() << std::endl;
                }
            }
            // Check if line contains segment information
            else if (std::regex_search(line, matches, segmentPattern) && matches.size() == 3) {
                try {
                    int start = std::stoi(matches[1].str());
                    int end = std::stoi(matches[2].str());

                    segments.push_back({ start, end });
                }
                catch (const std::exception& e) {
                    std::cerr << "Error parsing segment: " << line << std::endl;
                    std::cerr << "Exception: " << e.what() << std::endl;
                }
            }
        }

        file.close();
        return { nodes, segments };
    }

    Surface read_surface(std::string filename) {
        GeometryCentralMeshData mesh_data = file_to_geometrycentral_data(filename);

        auto faces = mesh_data.mesh->getFaceVertexList();
        std::unique_ptr<ManifoldSurfaceMesh> manifold_mesh = std::make_unique<ManifoldSurfaceMesh>(faces);

        VertexData<Vector3> vertex_positions(*manifold_mesh);
        for (Vertex v : manifold_mesh->vertices()) {
            // Use the integer index of the vertex to look up the position data
            // from the original mesh's geometry container. This is the correct way
            // to transfer data between two different mesh objects.
            // Otherwise some internal geometrycentral assertion fails...
            vertex_positions[v] = mesh_data.geometry->inputVertexPositions[v.getIndex()];
        }

        Surface surface;
        surface.mesh = std::move(manifold_mesh);
        surface.vertexPositions = std::move(vertex_positions);

        return surface;
    }
}