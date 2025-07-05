#include <remesh_surface.h>
#include "geometrycentral/surface/remeshing.h"

namespace modules {

	void remesh_surface(
		const Surface& surface,
		const double h
	) {
		RemeshOptions remeshOptions;
		remeshOptions.targetEdgeLength = h;
		adjustEdgeLengths(
			*surface.mesh,
			*surface.geometry,
			remeshOptions
		);
		surface.geometry->refreshQuantities();
	}
}