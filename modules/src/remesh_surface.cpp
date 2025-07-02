#include <remesh_surface.h>

namespace modules {

	void remesh_surface(
		const Surface& surface,
		const double h
	) {
		// According to https://geometry-central.net/tutorials/basic_mutation/

        // Work directly with the input mesh and geometry (in-place modification)
        ManifoldSurfaceMesh* mesh = surface.mesh.get();
        VertexPositionGeometry* geometry = surface.geometry.get();

        // Track original vertices and edges
        VertexData<bool> isOrigVert(*mesh, true);
        EdgeData<bool> isOrigEdge(*mesh, true);
        std::vector<Edge> toFlip;

        // Phase 1: Split all original edges
        for (Edge e : mesh->edges()) {
            if (!isOrigEdge[e]) continue; // Skip newly created edges

            // Get edge endpoints and their positions
            Vertex oldA = e.halfedge().tipVertex();
            Vertex oldB = e.halfedge().tailVertex();
            Vector3 oldAPos = geometry->vertexPositions[oldA];
            Vector3 oldBPos = geometry->vertexPositions[oldB];

            // Split the edge
            Vertex newV = mesh->splitEdgeTriangular(e).vertex();
            isOrigVert[newV] = false;

            // Position new vertex at edge midpoint
            Vector3 newPos = 0.5 * (oldAPos + oldBPos);
            geometry->inputVertexPositions[newV] = newPos;

            // Mark new edges and collect edges to flip
            for (Edge newE : newV.adjacentEdges()) {
                isOrigEdge[newE] = false;
                Vertex otherV = newE.otherVertex(newV);

                // Mark edges between old and new vertices for flipping
                if (isOrigVert[otherV] && otherV != oldA && otherV != oldB) {
                    toFlip.push_back(newE);
                }
            }
        }

        // Phase 2: Flip edges connecting old and new vertices
        for (Edge e : toFlip) {
            mesh->flip(e);
        }

        // ToDo: combine too large faces

        // Refresh geometry
        geometry->refreshQuantities();
	}
}