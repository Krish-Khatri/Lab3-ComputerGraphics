#pragma once

#include "Common/Matrix.h"

class Transform {
public:
	static Mat4 translate(const Vec3& w);
	static Mat4 scale(const Pt3& q, const Vec3& u, double s);
	static Mat4 rotOnly(const Mat4& m);
};