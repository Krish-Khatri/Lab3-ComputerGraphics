#ifndef SHADE_AND_SHAPES_H
#define SHADE_AND_SHAPES_H

#include "Rendering/Operator.h"
#include "Rendering/Geometry.h"
#include "Rendering/Scene.h"

#include "Common/Common.h"
#include "Common/Matrix.h"
#include "Rendering/Scene.h"

class Material : public SceneObject {
protected:
	Color _ambient, _diffuse, _specular;
	double _specExp, _reflective, _transp, _refractInd;
public:
	inline Material() : _specExp(1), _reflective(0), _transp(0), _refractInd(1) {};

	inline const Color& getAmbient() const { return _ambient; }
	inline const Color& getDiffuse() const { return _diffuse; }
	inline const Color& getSpecular() const { return _specular; }
	inline double getSpecExponent() const { return _specExp; }
	inline double getReflective() const { return _reflective; }
	inline double getTransparency() const { return _transp; }
	inline double getRefractIndex() const { return _refractInd; }

	inline void setAmbient(const Color& amb) { _ambient = amb; }
	inline void setDiffuse(const Color& diff) { _diffuse = diff; }
	inline void setSpecular(const Color& spec) { _specular = spec; }
	inline void setSpecExponent(double s) { _specExp = s; }
	inline void setReflective(double r) { _reflective = r; }
	inline void setTransparency(double t) { _transp = t; }
	inline void setRefractIndex(double r) { _refractInd = r; }
	virtual void accept(SceneObjectVisitor* visitor, void* ret) { visitor->visit(this, ret); }
};

class Light : public SceneObject {
protected:
	unsigned int _id;
	Pt3 _pos;
	Color _color;
	Color _ambient; // Keeping ambient here makes it easier to code than using a global ambient
public:
	inline Light() : _color(Color(1, 1, 1)) {} // Default: white
	inline Light(const Pt3& pos, const Color& color) : _pos(pos), _color(color) {}

	inline const Pt3& getPos() const { return _pos; }
	inline const Color& getColor() const { return _color; }
	inline const Color& getAmbient() const { return _ambient; }
	inline unsigned int getId() { return _id; }

	inline void setPos(const Pt3& p) { _pos = p; }
	inline void setColor(const Color& c) { _color = c; }
	inline void setAmbient(const Color& c) { _ambient = c; }
	inline void setId(unsigned int id) { _id = id; }

	virtual void accept(SceneObjectVisitor* visitor, void* ret) { visitor->visit(this, ret); }
};


//========================================================================
// Shapes
//
// TODO: Finish the definitions for some of the shapes below
//========================================================================

class Sphere : public Geometry, public Operand {
protected:
	Pt3 _center;
	double _rad;
public:
	Sphere(const Pt3& c = ORIGIN, double r = 1) : _center(c), _rad(r) {}

	// set/get methods
	void setRadius(double r) { _rad = r; }
	double getRadius() { return _rad; }
	void setCenter(const Pt3& c) { _center = c; }

	Pt3 getCenter() { return _center; }
	void translate(const Vec3& trans);
	void rotate(double d, int axis);
	void updateTransform();

	virtual void accept(GeometryVisitor* visitor, void* ret) { visitor->visit(this, ret); }
	virtual void accept(SceneObjectVisitor* visitor, void* ret) { visitor->visit(this, ret); }
};

class Box : public Geometry, public Operand {
protected:
	// NOTE: The corner and center for boxes are defined for your convenience but may be redundant.
	// You are free to implement the box or any other shape however you want. Just make sure that it works.

	Pt3 _corner;
	Pt3 _center;
	Vec3 _lengthv;
	Vec3 _widthv;
	Vec3 _heightv;

	double _length;
	double _width;
	double _height;

public:
	Box(const Pt3& corner = ORIGIN,
		const Vec3& lv = Vec3(1, 0, 0, 0),
		const Vec3& wv = Vec3(0, 1, 0, 0),
		const Vec3& hv = Vec3(0, 0, 1, 0),
		double l = 1, double w = 1, double h = 1) {
		_corner = corner;
		_lengthv = lv;
		_widthv = wv;
		_heightv = hv;
		_length = l;
		_width = w;
		_height = h;
		updateTransform();
	}

	// set/get methods
	Vec3 getLengthVec() { return _lengthv; }
	Vec3 getWidthVec() { return _widthv; }
	Vec3 getHeightVec() { return _heightv; }
	Pt3 getCorner() { return _corner; }

	void setLengthVec(const Vec3& v) { _lengthv = v; }
	void setWidthVec(const Vec3& v) { _widthv = v; }
	void setHeightVec(const Vec3& v) { _heightv = v; }
	void setCorner(const Pt3& v) { _corner = v; }

	double getLength() { return _length; }
	double getWidth() { return _width; }
	double getHeight() { return _height; }

	void setLength(double l) { _length = l; }
	void setWidth(double w) { _width = w; }
	void setHeight(double h) { _height = h; }

	Pt3 getCenter() { return _center; }
	void translate(const Vec3& trans);
	void rotate(double d, int axis);
	void updateTransform();

	virtual void accept(GeometryVisitor* visitor, void* ret) { visitor->visit(this, ret); }
	virtual void accept(SceneObjectVisitor* visitor, void* ret) { visitor->visit(this, ret); }
};

// TODO: finish the definition of the following classes
class Ellipsoid : public Geometry, public Operand {
protected:
	Pt3 _center;
	double _su, _sv, _sw;
	Vec3 _u, _v, _w;
public:
	Ellipsoid(const Pt3& center = ORIGIN,
		const Vec3& u = Vec3(1, 0, 0, 0),
		const Vec3& v = Vec3(0, 1, 0, 0),
		const Vec3& w = Vec3(0, 0, 1, 0),
		double a = 1, double b = 1, double c = 1) {
		_center = center;
		_su = a;
		_sv = b;
		_sw = c;
		_u = u;
		_v = v;
		_w = w;
		updateTransform();
	}
	Pt3 getCenter() { return _center; }
	void setCenter(const Pt3& c) { _center = c; }
	double getA() { return _su; }
	void setA(double a) { _su = a; }
	double getB() { return _sv; }
	void setB(double b) { _sv = b; }
	double getC() { return _sw; }
	void setC(double c) { _sw = c; }
	Vec3 getUAxis() { return _u; }
	void setUAxis(const Vec3& vec) { _u = vec; _u.normalize(); }
	Vec3 getVAxis() { return _v; }
	void setVAxis(const Vec3& vec) { _v = vec; _v.normalize(); }
	Vec3 getWAxis() { return _w; }
	void setWAxis(const Vec3& vec) { _w = vec; _w.normalize(); }
	void translate(const Vec3& trans);
	void rotate(double d, int axis);
	void updateTransform();

	virtual void accept(GeometryVisitor* visitor, void* ret) { visitor->visit(this, ret); }
	virtual void accept(SceneObjectVisitor* visitor, void* ret) { visitor->visit(this, ret); }
};

class Cylinder : public Geometry, public Operand {
protected:
	Pt3 _q;
	Vec3 _u, _v, _a;
	double _su, _sv, _h;
public:
	Cylinder(const Pt3& q = ORIGIN,
		const Pt3& u = Vec3(1, 0, 0, 0),
		const Pt3& v = Vec3(0, 1, 0, 0),
		const Vec3& a = Vec3(0, 0, 1, 0),
		double su = 1, double sv = 1, double h = 1)
		: _q(q), _a(a), _u(u), _v(v), _su(su), _sv(sv), _h(h)
	{
		updateTransform();
	}

	Pt3 getCenter() { return _q + _a * (_h / 2); }

	Pt3 getQ() { return Pt3(_q); }
	void setQ(const Pt3& a) { _q = a; }

	Vec3 getUAxis() { return _u; }
	void setUAxis(const Vec3& vec) { _u = vec; _u.normalize(); }

	Vec3 getVAxis() { return _v; }
	void setVAxis(const Vec3& vec) { _v = vec; _v.normalize(); }

	Vec3 getA() { return Pt3(_a); }
	void setA(const Vec3& a) { _a = a; _a.normalize(); }

	double getULength() { return _su; }
	void setULength(double a) { _su = a; }

	double getVLength() { return _sv; }
	void setVLength(double b) { _sv = b; }

	double getHeight() { return _h; }
	void setHeight(double b) { _h = b; }

	void translate(const Vec3& trans);
	void rotate(double d, int axis);
	void updateTransform();

	virtual void accept(GeometryVisitor* visitor, void* ret) { visitor->visit(this, ret); }
	virtual void accept(SceneObjectVisitor* visitor, void* ret) { visitor->visit(this, ret); }
};

class Cone : public Geometry, public Operand
{
protected:
	Pt3 _q;
	Vec3 _u, _v, _a;
	double _su, _sv, _h;
public:
	Cone(const Pt3& q = ORIGIN,
		const Pt3& u = Vec3(1, 0, 0, 0),
		const Pt3& v = Vec3(0, 1, 0, 0),
		const Vec3& a = Vec3(0, 0, 1, 0),
		double su = 1, double sv = 1, double h = 1)
		: _q(q), _a(a), _u(u), _v(v), _su(su), _sv(sv), _h(h) {
		updateTransform();
	}

	Pt3 getCenter() {
		// COM
		return _q + _a * (_h / 2);
	}

	Pt3 getQ() { return Pt3(_q); }
	void setQ(const Pt3& a) { _q = a; }

	Vec3 getUAxis() { return _u; }

	void setUAxis(const Vec3& vec) {
		_u = vec;
		_u.normalize();
	}

	Vec3 getVAxis() { return _v; }

	void setVAxis(const Vec3& vec) {
		_v = vec;
		_v.normalize();
	}

	Vec3 getA() { return Pt3(_a); }

	void setA(const Vec3& a) {
		_a = a;
		_a.normalize();
	}

	double getULength() { return _su; }
	void setULength(double a) { _su = a; }

	double getVLength() { return _sv; }
	void setVLength(double b) { _sv = b; }

	double getHeight() { return _h; }
	void setHeight(double b) { _h = b; }

	void translate(const Vec3& trans);
	void rotate(double d, int axis);
	void updateTransform();

	virtual void accept(GeometryVisitor* visitor, void* ret) { visitor->visit(this, ret); }
	virtual void accept(SceneObjectVisitor* visitor, void* ret) { visitor->visit(this, ret); }
};


//========================================================================
// END OF Shapes
//========================================================================


struct IsectData {
	// NOTE: You are free to modify this struct
	bool hit;
	bool inside;		// Did intersection happen inside the object?
	double t;
	Vec3 normal;
	IsectData() : hit(false), t(DINF), inside(false) {}
};

struct IsectAxisData {
	bool hit;
	int axis;
	IsectAxisData() : hit(false) {}
};

class Intersector : public GeometryVisitor {
protected:
	Ray _ray;
public:
	Intersector() : _ray(Pt3(0, 0, 0), Vec3(1, 0, 0, 0)) {}
	Intersector(const Ray& r) : _ray(r) {}
	void setRay(const Ray& r) { _ray = r; }
	virtual void visit(Sphere* sphere, void* ret);
	virtual void visit(Box* op, void* ret);
	virtual void visit(Ellipsoid* op, void* ret);
	virtual void visit(Cylinder* op, void* ret);
	virtual void visit(Cone* op, void* ret);
	virtual void visit(Operator* op, void* ret);
};

#endif