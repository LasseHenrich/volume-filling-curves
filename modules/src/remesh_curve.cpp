#include <remesh_curve.h>
#include <vector>
#include <map>

namespace modules {
	std::tuple<
		std::vector<Vector3>, // nodes
		std::vector<std::array<int, 2>>, // segments
		std::vector<double> // segment lengths
	> removeShortEdges(
		const std::vector<Vector3>& nodes,
		const std::vector<std::array<int, 2>>& segments,
		const std::vector<double>& segmentLengths,
		const double h
	) {
		std::map<int, bool> deletingNodes;

		std::vector<std::vector<int>> node2Segments(nodes.size(), std::vector<int>{});
		for (int i = 0; i < segments.size(); i++) {
			node2Segments[segments[i][0]].emplace_back(i);
			node2Segments[segments[i][1]].emplace_back(i);
		}

		// 1. determine which nodes to delete
		for (int i = 0; i < segments.size(); i++) {
			double edgeLen = segmentLengths[i];

			if (edgeLen < h) {
				double shortestLen = INFINITY;
				int shortestV = -1;

				for (int j = 0; j < 2; j++) {
					int v = segments[i][j];

					int otherS = node2Segments[v][0] == i ? node2Segments[v][1] : node2Segments[v][0];

					if (segmentLengths[otherS] < shortestLen && !deletingNodes[v]) {
						shortestLen = segmentLengths[otherS];
						shortestV = v;
					}
				}

				if (shortestV != -1) {
					deletingNodes[shortestV] = true;
				}
			}
		}

		std::map<int, int> node2NewNode;
		std::vector<Vector3> newNodes = {};

		for (int i = 0; i < nodes.size(); i++) {
			if (deletingNodes[i]) {
				continue;
			}
			node2NewNode[i] = newNodes.size();
			newNodes.push_back(nodes[i]);
		}

		// 2. update segments
		std::vector<std::array<int, 2>> newSegments = {};
		std::vector<double> newSegmentLengths = {};

		for (int i = 0; i < segments.size(); i++) {
			int v0 = segments[i][0], v1 = segments[i][1];

			if (!deletingNodes[v0] && !deletingNodes[v1]) {
				newSegments.push_back({ node2NewNode[v0], node2NewNode[v1] });
				newSegmentLengths.push_back(segmentLengths[i]);
			}
			else if (!deletingNodes[v0] && deletingNodes[v1]) {
				// traverse the curve until we find a node that is not deleted
				int v = v1;
				int currentSegment = i;
				while (deletingNodes[v]) {
					currentSegment = node2Segments[v][0] == currentSegment ? node2Segments[v][1] : node2Segments[v][0];
					v = segments[currentSegment][0] == v ? segments[currentSegment][1] : segments[currentSegment][0];
				}

				newSegments.push_back({ node2NewNode[v0], node2NewNode[v] });
				newSegmentLengths.push_back(segmentLengths[currentSegment]);
			}
		}

		return {
			newNodes,
			newSegments,
			newSegmentLengths
		};
	}

	std::tuple<
		std::vector<Vector3>, // nodes
		std::vector<std::array<int, 2>>, // segments
		std::vector<double> // segment lengths
	> subdivideSegments(
		const std::vector<Vector3>& nodes,
		const std::vector<std::array<int, 2>>& segments,
		const std::vector<double>& segmentLengths,
		const double h
	) {
		auto newNodes = nodes;
		auto newSegments = segments;
		auto newSegmentLengths = segmentLengths;

		for (int i = 0; i < segments.size(); i++) {
			double edgeLen = segmentLengths[i];

			if (edgeLen > 2 * h) {
				int divisionCount = std::ceil(edgeLen / h);
				double lenPerDivision = edgeLen / divisionCount;

				// add the new nodes
				std::vector<int> newNodeIndices = {};
				for (int j = 1; j < divisionCount; j++) {
					Vector3 newNode = nodes[segments[i][0]] + (nodes[segments[i][1]] - nodes[segments[i][0]]) * (j / (double)divisionCount);
					newNodes.push_back(newNode);
					newNodeIndices.push_back(newNodes.size() - 1);
				}
				// not adding the first and last node, as they are already in the newNodes vector

				// adapt the segments
				// first segment
				newSegments[i][1] = newNodeIndices[0];
				newSegmentLengths[i] = lenPerDivision;
				// intermediary segments
				for (int j = 0; j < newNodeIndices.size() - 1; j++) {
					newSegments.push_back({ newNodeIndices[j], newNodeIndices[j + 1] });
					newSegmentLengths.push_back(lenPerDivision);
				}
				// last segment
				newSegments.push_back({ newNodeIndices[newNodeIndices.size() - 1], segments[i][1] });
				newSegmentLengths.push_back(lenPerDivision);
			}
		}

		return {
			newNodes,
			newSegments,
			newSegmentLengths
		};
	}

	std::tuple<
		std::vector<Vector3>, // nodes
		std::vector<std::array<int, 2>>, // segments
		std::vector<double> // segment lengths
	> remesh_curve(
		const std::vector<Vector3> nodes,
		const std::vector<std::array<int, 2>>& segments,
		const std::vector<double>& segmentLengths,
		const double h
	) {
		// 1. remove nodes that are too close
		auto [newNodes, newSegments, newSegmentLengths] = removeShortEdges(
			nodes,
			segments,
			segmentLengths,
			h
		);

		// 2. subdivide segments
		// ToDo: This is bugged at the moment; it adds too many nodes
		std::tie(newNodes, newSegments, newSegmentLengths) = subdivideSegments(
			newNodes,
			newSegments,
			newSegmentLengths,
			h
		);

		// in case we have too little nodes or segments, we return the original data
		if (newNodes.size() < 3 || newSegments.size() < 3) {
			return {
			  nodes,
			  segments,
			  segmentLengths
			};
		}

		return {
		  newNodes,
		  newSegments,
		  newSegmentLengths
		};
	}
}