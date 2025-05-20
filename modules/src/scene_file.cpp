#include "scene_file.h"
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

    void processLine(SceneObject& scene, std::string directory, std::vector<std::string>& parts) {
        string key = parts[0];

        if (key == "#" || key == "//") {
            return;
        }

        if (key == "curve") {
            scene.curveFileName = directory + parts[1];
        }
        else if (key == "radius") {
            scene.radius = stod(parts[1]);
            scene.h = scene.radius * igl::PI / 20;
            scene.rmax = scene.radius * 5;
        }
        else if (key == "timestep") {
            scene.timestep = stod(parts[1]);
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

    SceneObject read_scene(std::string filename) {
        string directory = getDirectoryFromPath(filename);

        ifstream inFile;
        inFile.open(filename);

        if (!inFile) {
            cerr << "Could not open file " << filename << endl;
            exit(1);
        }

        SceneObject scene;

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
    > read_nodes(std::string filename) {
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
}