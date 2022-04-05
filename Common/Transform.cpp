#include "Common/Transform.h"

typedef Matrix<double, 3> Mat3;

// Returns the 3x3 matrix u^T * v
Mat3 transposeProd(const Vec3& u, const Vec3& v) {
	Mat3 res;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			res[i][j] = u[i] * v[j];
		}
	}
	return res;
}

Mat4 affine(const Mat3& rot, const Vec3& trans) {
	Mat4 res;

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			res[i][j] = rot[i][j];
		}
	}

	for (int i = 0; i < 3; i++) {
		res[3][i] = trans[i];
	}

	return res;
}

Mat4 Transform::translate(const Vec3& w) {
	return affine(Mat3::I, w);
}

Mat4 Transform::scale(const Pt3& q, const Vec3& u, double s) {

	return affine(Mat3::I + (s - 1) * transposeProd(u, u),
		(1 - s) * (q * u) * u);
}

Mat4 Transform::rotOnly(const Mat4& m) {
	Mat4 res(m);

	res[3][0] = res[3][1] = res[3][2] = 0;

	return res;
}
