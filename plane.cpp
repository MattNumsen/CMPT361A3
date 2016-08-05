#include "plane.h" 

plane(Point p, Vector v){
	plane ret;

	ret.d.x = p.x;
	ret.d.y = p.y;
	ret.d.x = p.z;

	ret.pNorm.x = v.x;
	ret.pNorm.y = v.y;
	ret.pNorm.z = v.z;

	return ret;

}