#include "vector.h"

struct plane {
	Point d; //point on plane
	Vector pNorm; //surface normal of the plane
} Plane;

plane(Point p, Vector v);