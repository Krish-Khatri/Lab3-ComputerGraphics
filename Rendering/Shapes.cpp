#include "Common/Common.h" 
#include "Rendering/Shapes.h" 


void Sphere::translate(const Vec3& trans){
	_center+=trans;
}

// rotate a sphere has no effect
void Sphere::rotate(float r, int axis){}

void Box::translate(const Vec3& trans){
	_corner+=trans;
	updateTransform(); 
}

void Box::rotate(float d, int axis){ 
	// the rotation code for the other objects should look very similar.
	// first construct a rotation about x/y/z axis

	Mat4 m; 
	switch(axis){
		case OP_XAXIS: 
			m[1][1] = cos(d); 
			m[1][2] = -sin(d); 
			m[2][1] = sin(d); 
			m[2][2] = cos(d); 
			break; 
		case OP_YAXIS: 
			m[0][0] = cos(d); 
			m[0][2] = sin(d); 
			m[2][0] = -sin(d); 
			m[2][2] = cos(d); 
			break; 
		case OP_ZAXIS: 
			m[0][0] = cos(d); 
			m[0][1] = -sin(d); 
			m[1][0] = sin(d); 
			m[1][1] = cos(d); 
			break; 
		default: 
			break; 
	}

	Vec3 cv = _corner-_center; 

	// rotate the three linearly independent vectors that determine the shape
	_lengthv = _lengthv*m;
	_widthv = _widthv*m;
	_heightv = _heightv*m;

	// rotation is about the center, so points that are not the center will also be rotated; 
	_corner = _center + (cv*m); 

	updateTransform(); 
}
// TODO: fill in these functions
void Cylinder::translate(const Vec3& trans) {}
void Cylinder::rotate(float d, int axis) {}
void Cone::translate(const Vec3& trans) {}
void Cone::rotate(float d, int axis) {}
void Ellipsoid::translate(const Vec3& trans) {}
void Ellipsoid::rotate(float d, int axis) {}

// definition of the sphere can be pretty sparse.  You don't need to define the transform associated with a sphere.
void Sphere::updateTransform() {}

void Intersector::visit(Sphere* sphere, void* ret){
	// TODO: rewrite this as described in the textbook
	IsectData* iret = (IsectData*) ret; 

	Pt3 center = sphere->getCenter(); 
	float radius = sphere->getRadius(); 

	float closest = GeometryUtils::pointRayClosest(center,_ray); 
    Vec3 C2P = _ray.at(closest) - center; 
    float r2 = radius*radius; 
	float d2 = C2P*C2P; 

	if(d2 > r2 + EPS){
        iret->hit = false; 
		iret->t = 0; 
	}
	else{
		float h = sqrt(r2 - d2); 
		iret->t = closest-h; 
		iret->hit = true; 
		iret->normal = (_ray.at(iret->t) - center); 
		iret->normal.normalize(); 
	}
}


// The updateTransform functions should properly set up the transformation of this geometry.
// The transformation takes a unit shape into the shape described by the parameters.
// This function also computes the inverse of the transformation.
void Box::updateTransform(){
	Vec3 ncenter = _corner; 
	ncenter += _length/2 * _lengthv; 
	ncenter += _width/2 * _widthv; 
	ncenter += _height/2 * _heightv; 
	_center = ncenter; 

	_mat[0][0] = _lengthv[0]*_length;
	_mat[0][1] = _lengthv[1]*_length;
	_mat[0][2] = _lengthv[2]*_length;
	_mat[0][3] = 0; 
	_mat[1][0] = _widthv[0]*_width;
	_mat[1][1] = _widthv[1]*_width;
	_mat[1][2] = _widthv[2]*_width;
	_mat[1][3] = 0; 
	_mat[2][0] = _heightv[0]*_height;
	_mat[2][1] = _heightv[1]*_height;
	_mat[2][2] = _heightv[2]*_height;
	_mat[2][3] = 0; 
	_mat[3][0] = _corner[0];
	_mat[3][1] = _corner[1];
	_mat[3][2] = _corner[2];
	_mat[3][3] = 1;

	_imat = !_mat; 

	Geometry::updateTransform(); // must call this at the end
}

// a box has six faces, which are basically six planes with rectangular boundaries.
// TODO: still need to return the normal at the intersection// still need to return the normal at the intersection
void Intersector::visit(Box* op, void* ret){

	IsectData* iret = (IsectData*) ret; 

	//if(true) {iret->hit=false; return;}

	Pt3 np = _ray.p*op->getInverseMat(); 
	Pt3 ppp = np*(1/np[3]); 
	Vec3 vvv = _ray.dir*op->getInverseMat(); 

	Ray ray(ppp,vvv); 

	Plane p(Pt3(.5,.5,0),Vec3(0,0,1,0)); 
	float d = GeometryUtils::planeRay(p,ray); 
	Pt3 pp = ray.at(d); 

	iret->hit = false; 
	iret->t = FINF32; 

	if(pp[0]>=0&&pp[0]<=1&&pp[1]>=0&&pp[1]<=1){
		iret->t = min(iret->t,d); 
		iret->hit = true; 
	}

	p = Plane(Pt3(.5,.5,1),Vec3(0,0,-1,0)); 
	d = GeometryUtils::planeRay(p,ray); 
	pp = ray.at(d); 

	if(pp[0]>=0&&pp[0]<=1&&pp[1]>=0&&pp[1]<=1){
		iret->t = min(iret->t,d); 
		iret->hit = true; 
	}

	p = Plane(Pt3(.5,0,.5),Vec3(0,1,0,0)); 
	d = GeometryUtils::planeRay(p,ray); 
	pp = ray.at(d); 

	if(pp[0]>=0&&pp[0]<=1&&pp[2]>=0&&pp[2]<=1){
		iret->t = min(iret->t,d); 
		iret->hit = true; 
	}

	p = Plane(Pt3(.5,1,.5),Vec3(0,-1,0,0)); 
	d = GeometryUtils::planeRay(p,ray); 
	pp = ray.at(d); 

	if(pp[0]>=0&&pp[0]<=1&&pp[2]>=0&&pp[2]<=1){
		iret->t = min(iret->t,d); 
		iret->hit = true; 
	}

	p = Plane(Pt3(0,.5,.5),Vec3(1,0,0,0)); 
	d = GeometryUtils::planeRay(p,ray); 
	pp = ray.at(d); 

	if(pp[1]>=0&&pp[1]<=1&&pp[2]>=0&&pp[2]<=1){
		iret->t = min(iret->t,d); 
		iret->hit = true; 
	}

	p = Plane(Pt3(1,.5,.5),Vec3(-1,0,0,0)); 
	d = GeometryUtils::planeRay(p,ray); 
	pp = ray.at(d); 

	if(pp[1]>=0&&pp[1]<=1&&pp[2]>=0&&pp[2]<=1){
		iret->t = min(iret->t,d); 
		iret->hit = true; 
	}
}

// TODO: need to fill in the following functions
void Ellipsoid::updateTransform() {}
void Intersector::visit(Ellipsoid* op, void* ret){}

void Cylinder::updateTransform() {}
void Intersector::visit(Cylinder* op, void* ret){}

void Cone::updateTransform() {}
void Intersector::visit(Cone* op, void* ret){}

// The operator is the widget that allows you to translate and rotate a geometric object
// It is colored as red/green/blue.  When one of the axis is highlighted, it is yellow.
void Intersector::visit(Operator* op, void* ret){
	IsectAxisData* iret = (IsectAxisData*) ret; 
	Pt3 center = op->getPrimaryOp()->getCenter(); 

	iret->hit = false; 

	if(op->getState()==OP_TRANSLATE){
		Ray rx(center,op->getDirX()); 
		Ray ry(center,op->getDirY()); 
		Ray rz(center,op->getDirZ()); 

		float dx = GeometryUtils::pointRayDist(rx.p+rx.dir*.5,_ray); 
		float dy = GeometryUtils::pointRayDist(ry.p+ry.dir*.5,_ray); 
		float dz = GeometryUtils::pointRayDist(rz.p+rz.dir*.5,_ray); 

		float tol = .12f; 
		if(dx>tol && dy>tol && dz>tol) 
			iret->hit = false; 
		else{
			iret->hit = true; 
			if(dx<dy){
				if(dx<dz) iret->axis = OP_XAXIS; 
				else      iret->axis = OP_ZAXIS; 
			}
			else{
				if(dy<dz) iret->axis = OP_YAXIS; 
				else      iret->axis = OP_ZAXIS; 
			}
		}
	}
	else if(op->getState()==OP_ROTATE){
		Plane xy(center,Vec3(0,0,1,0)); 
		Plane xz(center,Vec3(0,1,0,0)); 
		Plane yz(center,Vec3(1,0,0,0)); 

		float dxy = GeometryUtils::planeRay(xy,_ray); 
		float dxz = GeometryUtils::planeRay(xz,_ray); 
		float dyz = GeometryUtils::planeRay(yz,_ray); 

		float txy = abs(mag(_ray.at(dxy)-center)-OP_STEP); 
		float txz = abs(mag(_ray.at(dxz)-center)-OP_STEP); 
		float tyz = abs(mag(_ray.at(dyz)-center)-OP_STEP); 

		float tol = .05f; 
		if(txy>tol && txz>tol && tyz>tol) 
			iret->hit = false; 
		else{
			iret->hit = true; 
			if(txy<txz){
				if(txy<tyz) iret->axis = OP_ZAXIS; 
				else      iret->axis = OP_XAXIS; 
			}
			else{
				if(txz<tyz) iret->axis = OP_YAXIS; 
				else      iret->axis = OP_XAXIS; 
			}
		}
	}
}
