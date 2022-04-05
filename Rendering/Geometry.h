#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "Common/Matrix.h"
#include "Rendering/Scene.h"
#include <iostream>

#define EPS 1e-6

class Ray {
public:
	Pt3 p;
	Vec3 dir;
	Pt3 at(double t) const { return p + (dir * t); }

	double param(const Pt3& q) const {
		return (q - p) * dir / (dir * dir);
	}

	Ray() : p(0, 0, 0), dir(1, 0, 0, 0) {};

	Ray(const Pt3& p, const Vec3& d) {
		this->p = p;
		this->dir = d;
	}
};

class Plane {
public:
	Pt3 p;
	Vec3 n;

	Plane() : p(0, 0, 0), n(1, 0, 0, 0) {}
	Plane(const Vec3& p, const Vec3& n) {
		this->p = p;
		this->n = n;
	}
};

class Sphere;
class Ellipsoid;
class Box;
class Cylinder;
class Cone;
class Operator;

class GeometryVisitor {
public:
	virtual void visit(Sphere* sphere, void* ret) = 0;
	virtual void visit(Ellipsoid* op, void* ret) = 0;
	virtual void visit(Box* op, void* ret) = 0;
	virtual void visit(Cylinder* op, void* ret) = 0;
	virtual void visit(Cone* op, void* ret) = 0;
	virtual void visit(Operator* op, void* ret) = 0;
};

class Geometry : public SceneObject {
protected:
	float _glmat[16];
	float _glimat[16];
	Mat4 _mat;  // Each geometry comes with an affine transformation on the shape
	Mat4 _imat; // Also stores the inverse, so when updating the affine transformation, you also need to update the inverse
public:
	Geometry() {}

	float* getGLInverseMat() { return _glimat; }
	float* getGLForwardMat() { return _glmat; }
	Mat4& getInverseMat() { return _imat; }
	Mat4& getForwardMat() { return _mat; }

	// Updates the library transformation matrices
	virtual void updateTransform() {
		for (int j = 0; j < 16; j++) {
			_glmat[j] = _mat[j >> 2][j & 3];
			_glimat[j] = _imat[j >> 2][j & 3];
		}
	}

	virtual void accept(GeometryVisitor* visitor, void* ret) = 0;
	virtual void accept(SceneObjectVisitor* visitor, void* ret) {}
};

class GeometryUtils {
public:
	static double pointRayClosest(const Pt3& pt, const Ray& ray);
	static double pointRayDist(const Pt3& pt, const Ray& ray);
	static double rayRayDist(const Ray& r1, const Ray& r2);
	static double lineSegRayDist(const Pt3& p0, const Pt3& p1, const Ray& r);
	static double planeRay(const Plane& pl, const Ray& r);
	static double planeRayDeg(const Plane& pl, const Vec3& xa, const Ray& r);
	static Vec3 reflect(const Vec3& v, const Vec3& n);
	static Vec3 refract(const Vec3& v, const Vec3& n, double n1, double n2);
};

#endif