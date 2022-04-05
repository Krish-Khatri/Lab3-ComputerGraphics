#include "Common/Common.h"
#include "Rendering/ShadeAndShapes.h"

#include "Common/Transform.h"

// Returns the rotation matrix for rotation around the "axis" by "d" radians
Mat4 axisRotationMatrix(double d, int axis) {
	int id0, id1;
	switch (axis) {
	case OP_XAXIS: id0 = 1; id1 = 2; break; // X -> Y,Z
	case OP_YAXIS: id0 = 2; id1 = 0; break; // Y -> Z,X
	case OP_ZAXIS: id0 = 0; id1 = 1; break; // Z -> X,Y
	default: assert(false);
	}

	Mat4 m; // identity
	double c = cos(d);
	double s = sin(d);
	m[id0][id0] = c;
	m[id0][id1] = -s;
	m[id1][id0] = s;
	m[id1][id1] = c;
	return move(m);
}

// Returns the rotation matrix for rotation around the "axis" by "d" radians
Mat4 axisRotationMatrix(double cT, double sT, const Vec3& axis) {
	Mat4 m; // identity
	double u1 = axis[0], u2 = axis[1], u3 = axis[2];
	m[0][0] = cT + (1 - cT) * u1 * u1;
	m[0][1] = (1 - cT) * u1 * u2 + sT * u3;
	m[0][2] = (1 - cT) * u1 * u3 - sT * u2;
	m[1][0] = (1 - cT) * u1 * u2 - sT * u3;
	m[1][1] = cT + (1 - cT) * u2 * u2;
	m[1][2] = (1 - cT) * u2 * u3 + sT * u1;
	m[2][0] = (1 - cT) * u1 * u3 + sT * u2;
	m[2][1] = (1 - cT) * u2 * u3 - sT * u1;
	m[2][2] = cT + (1 - cT) * u3 * u3;
	return move(m);
}

//========================================================================
// translate() and rotate()
//========================================================================

void Sphere::translate(const Vec3& trans) {
	_center += trans;
	// No need to update the matrices for spheres
}
// Rotation doesn't affect a sphere
void Sphere::rotate(double r, int axis) {}

void Box::translate(const Vec3& trans) {
	_corner += trans;
	updateTransform(); // Must call this to update the inner matrices
}

void Box::rotate(double d, int axis) {
	// NOTE: The rotation code for the other objects should look very similar.
	Mat4 m = axisRotationMatrix(d, axis);
	Vec3 cv = _corner - _center;

	// rotate the three linearly independent vectors that determine the shape
	_lengthv = _lengthv * m;
	_widthv = _widthv * m;
	_heightv = _heightv * m;

	// rotation is about the center, so points that are not the center will also be rotated;
	_corner = _center + (cv * m);

	updateTransform();
}

// TODO: fill in these functions
void Cylinder::translate(const Vec3& trans) {
	_q += trans;
	updateTransform();
}
void Cylinder::rotate(double d, int axis) {
	Mat4 m = axisRotationMatrix(d, axis);
	Pt3 center = getCenter();
	Vec3 ax = _q - center;
	ax = ax * m;
	_q = center + ax;

	_u = _u * m; _u.normalize();
	_v = _v * m; _v.normalize();
	_a = _a * m; _a.normalize();

	updateTransform();
}
void Cone::translate(const Vec3& trans) {
	_q += trans;
	updateTransform();
}
void Cone::rotate(double d, int axis) {
	Mat4 m = axisRotationMatrix(d, axis);
	Pt3 center = getCenter();
	Vec3 ax = _q - center;
	ax = ax * m;
	_q = center + ax;

	_u = _u * m; _u.normalize();
	_v = _v * m; _v.normalize();
	_a = _a * m; _a.normalize();

	updateTransform();
}
void Ellipsoid::translate(const Vec3& trans) {
	_center += trans;
	updateTransform();
}
void Ellipsoid::rotate(double d, int axis) {
	// NOTE: The rotation code for the other objects should look very similar.
	Mat4 m = axisRotationMatrix(d, axis);

	// rotate the three linearly independent vectors that determine the shape
	_u = _u * m; _u.normalize();
	_v = _v * m; _v.normalize();
	_w = _w * m; _w.normalize();

	updateTransform();
}


//========================================================================
// updateTransform() and Intersector::visit()
//========================================================================


// The operator is the widget that allows you to translate and rotate a geometric object
// It is colored as red/green/blue.  When one of the axis is highlighted, it becomes yellow.
void Intersector::visit(Operator* op, void* ret) {
	IsectAxisData* iret = (IsectAxisData*)ret;
	Pt3 center = op->getPrimaryOp()->getCenter();
	iret->hit = false;
	const int axes[3] = { OP_XAXIS, OP_YAXIS, OP_ZAXIS };

	if (op->getState() == OP_TRANSLATE) {
		const Ray rays[3] = {
			Ray(center, op->getDirX()), // X
			Ray(center, op->getDirY()), // Y
			Ray(center, op->getDirZ()) // Z
		};

		double bestTime = 0.1;
		for (int i = 0; i < 3; i++) {
			double hitTime = GeometryUtils::pointRayDist(rays[i].p + 0.5 * rays[i].dir, _ray);
			if (bestTime > hitTime) {
				bestTime = hitTime;
				iret->hit = true;
				iret->axis = axes[i];
			}
		}
	}
	else if (op->getState() == OP_ROTATE) {
		const Plane planes[3] = {
			Plane(center, Vec3(1, 0, 0, 0)), // X: YZ-plane
			Plane(center, Vec3(0, 1, 0, 0)), // Y: XZ-plane
			Plane(center, Vec3(0, 0, 1, 0)) // Z: XY-plane
		};

		double bestTime = 0.1;
		for (int i = 0; i < 3; i++) {
			double timeTemp = GeometryUtils::planeRay(planes[i], _ray);
			double hitTime = abs(mag(_ray.at(timeTemp) - center) - OP_STEP);
			if (bestTime > hitTime) {
				bestTime = hitTime;
				iret->hit = true;
				iret->axis = axes[i];
			}
		}
	}
}


// The definition of a sphere can be pretty sparse,
// so you don't need to define the transform associated with a sphere.
void Sphere::updateTransform() {}

// TODO: rewrite this function as described in the textbook
void Intersector::visit(Sphere* sphere, void* ret) {
	IsectData* iret = (IsectData*)ret;

	iret->hit = false;
	iret->t = 0;
	iret->inside = false;

	Pt3 center = sphere->getCenter();
	double radius = sphere->getRadius();

	double closest = GeometryUtils::pointRayClosest(center, _ray);
	Vec3 C2P = _ray.at(closest) - center;
	double r2 = radius * radius;
	double d2 = C2P * C2P;

	if (d2 > r2 + EPS) {
		return;
	}

	Vec3 v = _ray.dir / mag(_ray.dir);

	double h = sqrt(r2 - d2);
	Pt3 a = center + C2P - h * v;
	double t = _ray.param(a);

	if (t < 0) {
		Pt3 a = center + C2P + h * v;
		t = _ray.param(a);
		iret->inside = true;
	}

	if (t < 0) {
		return;
	}

	iret->t = t;
	iret->hit = true;
	iret->normal = (_ray.at(iret->t) - center);
	iret->normal.normalize();
}


// The updateTransform functions should properly set up the transformation of this geometry.
// The transformation takes a unit shape into the shape described by the parameters.
// This function also computes the inverse of the transformation.
void Box::updateTransform() {
	Vec3 ncenter = _corner;
	ncenter += _length / 2 * _lengthv;
	ncenter += _width / 2 * _widthv;
	ncenter += _height / 2 * _heightv;
	_center = ncenter;

	_mat[0][0] = _lengthv[0] * _length;
	_mat[0][1] = _lengthv[1] * _length;
	_mat[0][2] = _lengthv[2] * _length;
	_mat[0][3] = 0;
	_mat[1][0] = _widthv[0] * _width;
	_mat[1][1] = _widthv[1] * _width;
	_mat[1][2] = _widthv[2] * _width;
	_mat[1][3] = 0;
	_mat[2][0] = _heightv[0] * _height;
	_mat[2][1] = _heightv[1] * _height;
	_mat[2][2] = _heightv[2] * _height;
	_mat[2][3] = 0;
	_mat[3][0] = _corner[0];
	_mat[3][1] = _corner[1];
	_mat[3][2] = _corner[2];
	_mat[3][3] = 1;

	// NOTE: These two lines are required at the end
	_imat = !_mat;
	Geometry::updateTransform();
}

// A box has six faces, which are basically six planes with rectangular boundaries.
// TODO: Compute the face normal
void Intersector::visit(Box* op, void* ret) {

	// Initialization
	IsectData* iret = (IsectData*)ret;
	iret->hit = false; // no collision
	iret->t = DINF;
	iret->inside = false;

	// Convert (original ray, current box) to (converted ray, canonical box)
	Mat4 invMat = op->getInverseMat();
	Pt3 newPoint = _ray.p * invMat; newPoint /= newPoint[3];
	// NOTE: Do not normalize this vector, or hit time will be wrong
	Vec3 newDir = _ray.dir * invMat; newDir[3] = 0;
	Ray ray(newPoint, newDir);

	// Canonical box: unit cube, axis-aligned, all coordinates within [0,1]
	const Plane planes[6] = {
		Plane(Pt3(.5, .5, 0), Vec3(0, 0, -1, 0)), // Bottom (-z)
		Plane(Pt3(.5, .5, 1), Vec3(0, 0, 1, 0)), // Top (+z)
		Plane(Pt3(.5, 0, .5), Vec3(0, -1, 0, 0)), // Left (-y)
		Plane(Pt3(.5, 1, .5), Vec3(0, 1, 0, 0)), // Right (+y)
		Plane(Pt3(0, .5, .5), Vec3(-1, 0, 0, 0)), // Back (-x)
		Plane(Pt3(1, .5, .5), Vec3(1, 0, 0, 0)) // Front (+x)
	};
	// Non-constant axes: {x,y} for z, {x,z} for y, and {y,z} for x
	const int Axis0[6] = { 0, 0, 0, 0, 1, 1 };
	const int Axis1[6] = { 1, 1, 2, 2, 2, 2 };

	// Ray-intersection with the canonical box
	for (int i = 0; i < 6; i++) {
		const Plane& pl = planes[i];
		double hitTime = GeometryUtils::planeRay(pl, ray);
		Pt3 hitPoint = ray.at(hitTime);
		// positive hit time, closer to the "eye point", and hit point within the unit square
		if (hitTime > 0 && iret->t > hitTime &&
			hitPoint[Axis0[i]] >= 0 && hitPoint[Axis0[i]] <= 1 &&
			hitPoint[Axis1[i]] >= 0 && hitPoint[Axis1[i]] <= 1)
		{
			iret->hit = true;
			iret->t = hitTime;
			iret->normal = pl.n; // Canonical-box normal
		}
	}

	if (iret->hit) {
		// TODO: Compute the face normal for the hit plane (canonical box -> current box)
		iret->normal = iret->normal * op->getForwardMat();
		iret->normal.normalize();
		iret->inside = 0 < newPoint[0] && newPoint[0] < 1 &&
			0 < newPoint[1] && newPoint[1] < 1 &&
			0 < newPoint[2] && newPoint[2] < 1;
	}
}

// TODO: need to fill in the following functions
void Ellipsoid::updateTransform() {
	_mat = Transform::scale(ORIGIN, _u, _su) *
		Transform::scale(ORIGIN, _v, _sv) *
		Transform::scale(ORIGIN, _w, _sw) *
		Transform::translate(_center);
	// NOTE: These two lines are required at the end
	_imat = !_mat;
	Geometry::updateTransform();
}

void Intersector::visit(Ellipsoid* ellipsoid, void* ret) {
	// Initialization
	IsectData* iret = (IsectData*)ret;
	iret->hit = false; // no collision
	iret->t = DINF;
	Mat4 invMat = ellipsoid->getInverseMat();
	Pt3 newPoint = _ray.p * invMat;
	// NOTE: Do not normalize this vector, or hit time will be wrong
	Vec3 newDir = _ray.dir * invMat;
	Ray ray(newPoint, newDir);

	Intersector intersector(ray);

	Sphere sphere(ORIGIN, 1.0);

	intersector.visit(&sphere, iret);

	if (iret->hit) {
		Mat4 m_u = Transform::rotOnly(ellipsoid->getForwardMat());
		iret->normal = iret->normal * transpose(!m_u);
		iret->normal.normalize();
	}
}

void Cylinder::updateTransform() {
	Vec3 axis = cross(Vec3(0, 0, 1, 0), _a);
	double sT = mag(axis);
	double cT = Vec3(0, 0, 1, 0) * _a;
	axis.normalize();
	_mat = axisRotationMatrix(cT, sT, axis)
		* Transform::scale(ORIGIN, _u, _su)
		* Transform::scale(ORIGIN, _v, _sv)
		* Transform::scale(ORIGIN, _a, _h)
		* Transform::translate(_q);
	_imat = !_mat;
	Geometry::updateTransform();
}

void intersectUnitCylinder(const Ray& _ray, IsectData* iret) {
	iret->hit = false; // no collision
	iret->t = DINF;

	iret->inside = 0 < _ray.p[2] && _ray.p[2] < 1 && // Within XY planes
		_ray.p[0] * _ray.p[0] + _ray.p[1] * _ray.p[1] < 1;	// Within circle on XY plane

	Pt3 P0 = ORIGIN;
	Pt3 P1 = P0 + Vec3(0, 0, 1, 0);
	Vec3 A = P1 - P0;
	double h = mag(A); A.normalize();
	double radius = 1.0;
	Ray axisRay(P0, A);
	if (abs(_ray.dir * A) <= EPS) {
		Plane centerPlane((P0 + P1) / 2, A);
		double d = (centerPlane.p - _ray.p) * centerPlane.n;
		if (abs(d) > h / 2 - EPS)return;
		Plane pl(_ray.p, A);
		double t = GeometryUtils::planeRay(pl, axisRay);
		Pt3 c = axisRay.at(t);
		Vec3 dirP_unit = _ray.dir / mag(_ray.dir);
		Vec3 w = (_ray.p - c) - ((_ray.p - c) * dirP_unit) * dirP_unit;
		double x = radius * radius - w * w;
		if (x < EPS)return;
		x = sqrt(x);
		Pt3 Q = c + w;
		Pt3 D1 = Q - x * dirP_unit;
		Pt3 D2 = Q + x * dirP_unit;
		t = _ray.param(D1);
		if (t < iret->t && t > 0) {
			iret->hit = true;
			iret->t = t;
			iret->normal = D1 - c;
			iret->normal.normalize();
		}
		t = _ray.param(D2);
		if (t < iret->t && t > 0) {
			iret->hit = true;
			iret->t = t;
			iret->normal = D2 - c;
			iret->normal.normalize();
		}
		return;
	}

	Plane bottom(P0, -A);
	Plane top(P1, A);

	double bottomHT = GeometryUtils::planeRay(bottom, _ray);
	Pt3 E1 = _ray.at(bottomHT);
	Vec3 distB = E1 - P0;
	if (distB * distB <= radius * radius + EPS && bottomHT > 0) {
		iret->hit = true;
		iret->t = bottomHT;
		iret->normal = bottom.n;
		iret->normal.normalize();
	}
	double topHT = GeometryUtils::planeRay(top, _ray);
	Pt3 E2 = _ray.at(topHT);
	Vec3 distT = E2 - P1;
	if (distT * distT <= radius * radius + EPS && topHT < iret->t && topHT > 0) {
		iret->hit = true;
		iret->t = topHT;
		iret->normal = top.n;
		iret->normal.normalize();
	}

	if (distT * distT <= radius * radius + EPS && distB * distB <= radius * radius + EPS) {
		return;
	}

	Ray ray(E1, E2 - E1);
	double dist = GeometryUtils::rayRayDist(ray, axisRay);
	if (dist - radius > EPS) {
		return;
	}
	else if (abs(ray.dir * A - mag(ray.dir)) <= EPS) {
		return;
	}
	Vec3 dirP = ray.dir;
	dirP = dirP - (dirP * A) * A;
	Vec3 dirP_unit = dirP / mag(dirP);
	Vec3 w = (E1 - P0) - ((E1 - P0) * dirP_unit) * dirP_unit;
	Ray rayP(E1, dirP);
	double x = sqrt(radius * radius - w * w);
	Pt3 Q = P0 + w;
	Pt3 D1 = Q - x * dirP_unit;
	Pt3 D2 = Q + x * dirP_unit;
	double t1 = rayP.param(D1);
	double t2 = rayP.param(D2);
	if (t1 > -EPS && t1 < EPS + 1) {
		Pt3 ip = ray.at(t1);
		t1 = _ray.param(ip);
		if (t1 < iret->t && t1 > 0) {
			iret->hit = true;
			iret->t = t1;
			iret->normal = (ip - P0) - ((ip - P0) * A) * A; iret->normal.normalize();
		}
	}
	if (t2 > -EPS && t2 < EPS + 1) {
		Pt3 ip = ray.at(t2);
		t2 = _ray.param(ip);
		if (t2 < iret->t && t2 > 0) {
			iret->hit = true;
			iret->t = t2;
			iret->normal = (ip - P0) - ((ip - P0) * A) * A; iret->normal.normalize();
		}
	}
}

void Intersector::visit(Cylinder* op, void* ret) {
	IsectData* iret = (IsectData*)ret;

	Mat4 invMat = op->getInverseMat();
	Pt3 newPoint = _ray.p * invMat;
	// NOTE: Do not normalize this vector, or hit time will be wrong
	Vec3 newDir = _ray.dir * invMat;
	Ray ray(newPoint, newDir);

	intersectUnitCylinder(ray, iret);

	if (iret->hit) {
		Mat4 m_u = Transform::rotOnly(op->getForwardMat());
		iret->normal = iret->normal * transpose(!m_u);
		iret->normal.normalize();
	}
}

void Cone::updateTransform() {
	Vec3 axis = cross(Vec3(0, 0, 1, 0), _a);
	double sT = mag(axis);
	double cT = Vec3(0, 0, 1, 0) * _a;
	axis.normalize();
	_mat = axisRotationMatrix(cT, sT, axis)
		* Transform::scale(ORIGIN, _u, _su)
		* Transform::scale(ORIGIN, _v, _sv)
		* Transform::scale(ORIGIN, _a, _h)
		* Transform::translate(_q);
	_imat = !_mat;
	Geometry::updateTransform();
}

void intersectUnitCone(const Ray& _ray, IsectData* iret) {
	iret->hit = false; // no collision
	iret->t = DINF;
	iret->inside = false;
	double radius = 1;
	Pt3 v = Pt3(0, 0, 1);
	Pt3 _Q = ORIGIN;
	Vec3 A = _Q - v;
	double h = mag(A); A.normalize();
	Pt3 q = v + A;
	Pt3 p = _ray.p;
	Vec3 U = _ray.dir;
	Pt3 x0 = -p + ((q - p) * A) * v;
	double m0 = (v - p) * A;
	Pt3 x1 = ((v - q) * A) * U - (U * A) * v;
	double m1 = -U * A;
	Ray ray_trans;
	double hT = GeometryUtils::planeRay(Plane(_Q, A), _ray);
	double dist = mag(_ray.at(hT) - _Q);
	if (dist <= radius + EPS) {
		if (hT > 0 && hT < iret->t) {
			iret->hit = true;
			iret->inside = false;
			iret->t = hT;
			iret->normal = A; iret->normal.normalize();
		}
	}
	if (abs(m0) <= EPS && abs(m1) <= EPS) {
		return;
	}
	else if (abs(m1) <= EPS) {
		ray_trans = Ray(x0 / m0, x1);
	}
	else if (abs(m0) <= EPS) {
		ray_trans = Ray(x1 / m1, x0);
	}
	else {
		ray_trans = Ray(x0 / m0, x0 / m0 - x1 / m1);
	}
	Vec3 dirP = ray_trans.dir / mag(ray_trans.dir);
	Vec3 w = (ray_trans.p - q) - ((ray_trans.p - q) * dirP) * dirP;
	Pt3 tmp_q = q + w;
	double tmp_r = radius / h;
	double x = tmp_r * tmp_r - w * w;
	if (x < 0)return;
	Pt3 D2 = tmp_q + sqrt(x) * dirP;
	Pt3 D1 = tmp_q - sqrt(x) * dirP;
	double t1 = ((m0 * D1 - x0) * (x1 - m1 * D1)) / ((x1 - m1 * D1) * (x1 - m1 * D1));
	double t2 = ((m0 * D2 - x0) * (x1 - m1 * D2)) / ((x1 - m1 * D2) * (x1 - m1 * D2));
	if (t1 > 0 && iret->t > t1) {
		Pt3 hp = _ray.at(t1);
		double cT = (hp - v) * A / mag(hp - v);
		if (mag(hp - v) * cT < h && (hp - v) * A>0) {
			iret->inside = false;
			iret->hit = true;
			iret->t = t1;
			iret->normal = cT * (hp - v) - A * mag(hp - v); iret->normal.normalize();
		}
	}
	if (t2 > 0 && t2 < iret->t) {
		Pt3 hp = _ray.at(t2);
		double cT = (hp - v) * A / mag(hp - v);
		if (mag(hp - v) * cT < h && (hp - v) * A>0) {
			iret->inside = true;
			iret->hit = true;
			iret->t = t2;
			iret->normal = cT * (hp - v) - A * mag(hp - v); iret->normal.normalize();
		}
	}

	iret->inside = t1 < 0 || t2 < 0;
}

void Intersector::visit(Cone* op, void* ret) {
	IsectData* iret = (IsectData*)ret;

	Mat4 invMat = op->getInverseMat();
	Pt3 newPoint = _ray.p * invMat;
	// NOTE: Do not normalize this vector, or hit time will be wrong
	Vec3 newDir = _ray.dir * invMat;
	Ray ray(newPoint, newDir);

	intersectUnitCone(ray, iret);

	if (iret->hit) {
		Mat4 m_u = Transform::rotOnly(op->getForwardMat());
		iret->normal = iret->normal * transpose(!m_u);
		iret->normal.normalize();
	}
}
