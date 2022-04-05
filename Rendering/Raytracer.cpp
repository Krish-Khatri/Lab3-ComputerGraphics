#include "Rendering/Raytracer.h"
#include "Rendering/Shading.h"
#include <FL/glu.h>
#include "Common/Common.h"

Raytracer::Raytracer() {
	_pixels = NULL;
}

void Raytracer::drawInit(GLdouble modelview[16], GLdouble proj[16], GLint view[4]) {
	_width = view[2];
	_height = view[3];

	if (_pixels) delete[] _pixels;
	_pixels = new float[_width * _height * 4];

	memcpy(_modelview, modelview, 16 * sizeof(modelview[0]));
	memcpy(_proj, proj, 16 * sizeof(proj[0]));
	memcpy(_view, view, 4 * sizeof(view[0]));

	Mat4 mv, pr;
	for (int j = 0; j < 16; j++) {
		mv[j / 4][j % 4] = modelview[j];
		pr[j / 4][j % 4] = proj[j];
	}

	_final = mv * pr;
	_invFinal = !_final;
	_last = 0;
}

Pt3 Raytracer::unproject(const Pt3& p) {
	Pt3	np = p;
	np[0] = (p[0] - _view[0]) / _view[2];
	np[1] = (p[1] - _view[1]) / _view[3];

	np[0] = np[0] * 2 - 1;
	np[1] = np[1] * 2 - 1;
	np[2] = np[2] * 2 - 1;

	Pt3 ret = np * _invFinal;
	ret[0] /= ret[3];
	ret[1] /= ret[3];
	ret[2] /= ret[3];
	ret[3] = 1;
	return move(ret);
}

bool Raytracer::draw(int step) {
	int size = _width * _height;
	int end = min(size, _last + step);

#pragma omp parallel for schedule(static, 1)
	for (int j = _last; j < end; j++) {
		int x = j % _width;
		int y = j / _width;

		Pt3 rst = unproject(Pt3(x, y, 0));
		Pt3 red = unproject(Pt3(x, y, 1));

		Ray r(rst, red - rst);
		r.dir.normalize();

		TraceResult res;
		if (_scene)
			res = trace(r, 0);

		res.color[3] = 1;
		int offset = j * 4;
		for (int i = 0; i < 4; i++)
			_pixels[offset + i] = res.color[i];
	}

	_last = end;
	return (_last >= size);
}

double visibility(Scene* scene, const Pt3& start, const Pt3& end) {
	double s = 1;

	Ray r(start, end - start);
	Intersector intersector(r);

	for (int j = 0; j < scene->getNumObjects(); j++) {
		IsectData data;
		data.hit = false;

		Geometry* geom = scene->getObject(j);
		Material* mat = scene->getMaterial(geom);
		geom->accept(&intersector, &data);

		if (data.hit && data.t > EPS && data.t < 1 && !data.inside) {
			// Object is in the way
			s *= mat->getTransparency();
		}

		if (s < EPS) {
			return 0;
		}
	}

	return s;
}

TraceResult Raytracer::trace(const Ray& ray, int depth) {
	// TODO: perform recursive tracing here
	// TODO: perform proper illumination, shadowing, etc...
	TraceResult res;

	// Default background color if missing all objects (default)
	res.color = Color(0, 0, 0);

	if (depth > MAX_DEPTH) {
		return res;
	}

	const Color ambientI = _scene->getLight(0)->getAmbient();

	Intersector intersector(ray);

	double t = DINF;
	Vec3 normal;
	Material* material = NULL;
	bool inside = false;

	for (int j = 0; j < _scene->getNumObjects(); j++) {
		IsectData data;
		data.hit = false;

		Geometry* geom = _scene->getObject(j);
		Material* mat = _scene->getMaterial(geom);
		geom->accept(&intersector, &data);

		if (data.hit) {
			if (t > data.t) {
				t = data.t;
				inside = data.inside;
				normal = (inside ? -1 : 1) * data.normal;
				material = mat;
			}
		}
	}

	if (t < DINF) { // hit
		// NOTE: Modify as needed
		const Pt3 P = ray.at(t);

		// Do lighting and reflection only if we are in air
		if (!inside) {
			res.color += ambientI.scale(material->getAmbient());

			// Lighting
			for (int i = 0; i < _scene->getNumLights(); ++i) {
				const Light* light = _scene->getLight(i);
				const Pt3 lightPos = light->getPos();

				double shadow = visibility(_scene, P, lightPos);

				if (shadow <= EPS) {
					continue;
				}

				Vec3 L = lightPos - P;
				L.normalize();

				double diffuse = L * normal;

				if (diffuse > 0) {
					Vec3 nR = GeometryUtils::reflect(L, normal);
					nR.normalize();

					// nR = -R, ray.dir = -V
					double specular = pow(max(0, nR * ray.dir), material->getSpecExponent());
					res.color += shadow * (
						diffuse * light->getColor().scale(material->getDiffuse()) +
						specular * light->getColor().scale(material->getSpecular())
						);
				}
			}

			// Reflection
			double reflective = material->getReflective();

			if (reflective > 0) {
				Vec3 reflectDir = GeometryUtils::reflect(ray.dir, normal);
				Ray reflected(P + EPS * reflectDir, reflectDir);
				TraceResult result = this->trace(reflected, depth + 1);
				res.color += reflective * result.color;
			}
		}

		// Refraction
		double transparency = material->getTransparency();

		if (transparency > 0) {
			// TODO: replace 1 with refractive index of outside material
			// TODO: add inside support for all shapes
			Vec3 refractDir;
			if (!inside) {
				refractDir = GeometryUtils::refract(ray.dir, normal,
					1, material->getRefractIndex());
			}
			else {
				refractDir = GeometryUtils::refract(ray.dir, normal,
					material->getRefractIndex(), 1);
			}

			if (mag2(refractDir) > EPS) {
				Ray refracted(P + EPS * refractDir, refractDir);
				TraceResult result = this->trace(refracted, depth + 1);
				res.color += transparency * result.color;
			}
			else if (inside) {
				// Total internal reflection
				Vec3 reflectDir = GeometryUtils::reflect(ray.dir, normal);
				Ray reflected(P + EPS * reflectDir, reflectDir);
				TraceResult result = this->trace(reflected, depth + 1);
				res.color += result.color;
			}
		}
	}

	res.color.clamp01();

	return res;
}
