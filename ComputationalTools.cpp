#include "stdafx.h"

#include <iostream>
#include <fstream>
using namespace std;

#include "MeshObject.h"

void rotateSig(signature* sig){
	signature temp;

	temp = sig[1];
	sig[1] = sig[2];
	sig[2] = temp;

	temp = sig[7];
	sig[7] = sig[5];
	sig[5] = temp;

	temp = sig[3];
	sig[3] = sig[4];
	sig[4] = temp;

	temp = sig[6];
	sig[6] = sig[8];
	sig[8] = temp;

	float tempD;

	tempD = sig[0].dis58;
	sig[0].dis58 = sig[0].dis67;
	sig[0].dis67 = tempD;

	tempD = sig[0].dis56;
	sig[0].dis56 = sig[0].dis78;
	sig[0].dis78 = tempD;
}

Eigen::Quaternionf constructQuaternion(signature* sig1,signature* sig2){
	Eigen::Matrix3f a,b1,b2;
	a << sig1->d1.x(),sig1->d2.x(),sig1->normal.x(),
		sig1->d1.y(),sig1->d2.y(),sig1->normal.y(),
		sig1->d1.z(),sig1->d2.z(),sig1->normal.z();

	b1 << sig2->d1.x(),sig2->d2.x(),sig2->normal.x(),
		sig2->d1.y(),sig2->d2.y(),sig2->normal.y(),
		sig2->d1.z(),sig2->d2.z(),sig2->normal.z();
	b2 << -sig2->d1.x(),-sig2->d2.x(),sig2->normal.x(),
		-sig2->d1.y(),-sig2->d2.y(),sig2->normal.y(),
		-sig2->d1.z(),-sig2->d2.z(),sig2->normal.z();

	Eigen::Quaternionf r1(b1.transpose()*a);
	Eigen::Quaternionf r2(b2.transpose()*a);

	// 	cout<<r1.w()<<' '<<r1.x()<<' '<<r1.y()<<' '<<r1.z()<<endl;
	// 	cout<<r2.w()<<' '<<r2.x()<<' '<<r2.y()<<' '<<r2.z()<<endl;
	// 
	// 	int i;
	// 	cin>>i;

	if (r2.w()>r1.w())
		return r2;
	else
		return r1;
}

threeTuple getJetColor(float value) {
	float fourValue = 4 * value;

	float red   = min(fourValue - 1.5, -fourValue + 4.5);
	float green = min(fourValue - 0.5, -fourValue + 3.5);
	float blue  = min(fourValue + 0.5, -fourValue + 2.5);

	red = max_zenyo(0,min_zenyo(1,red));
	green = max_zenyo(0,min_zenyo(1,green));
	blue = max_zenyo(0,min_zenyo(1,blue));

	threeTuple color;
	color.x = red;
	color.y = green;
	color.z = blue;

	return color;
}

float computeLength(PolyhedralSurf::Halfedge_const_handle e){
	float dx = e->vertex()->point().x() - e->opposite()->vertex()->point().x();
	float dy = e->vertex()->point().y() - e->opposite()->vertex()->point().y();
	float dz = e->vertex()->point().z() - e->opposite()->vertex()->point().z();

	return sqrt(dx*dx+dy*dy+dz*dz);
}

float computeArea(PolyhedralSurf::Facet_const_handle f){
	if (f == NULL) return 0;

	PolyhedralSurf::Halfedge_const_handle h1,h2,h3;
	h1 = f->halfedge();
	h2 = h1->next();
	h3 = h2->next();

	float len1 = computeLength(h1);
	float len2 = computeLength(h2);
	float len3 = computeLength(h3);

	float s = (len1 + len2 + len3)/2;

	return sqrt(s*(s-len1)*(s-len2)*(s-len3));
}

void writeOFF(BasicMesh* mesh, string filename, int idx){

	ofstream file;
	file.open(filename.c_str());

	file<<"OFF"<<endl;
	file<<mesh->P.size_of_vertices()<<' '<<mesh->P.size_of_facets()<<" 0"<<endl;

	PolyhedralSurf::Vertex_const_iterator vb = mesh->P.vertices_begin();
	PolyhedralSurf::Vertex_const_iterator ve = mesh->P.vertices_end();
	PolyhedralSurf::Vertex_const_handle vh;

	std::map<Vertex_const_handle, float>::iterator iter;
	float color;

	for (int i = 0; vb != ve; vb++,i++){
		vh = vb;

		file<<vh->point().x() + Surface[idx]->u[i*3] <<' '<<vh->point().y() + Surface[idx]->u[i*3+1] <<' '<<vh->point().z() + Surface[idx]->u[i*3+2] <<endl;
	}

	for (int i = 0; i < mesh->faceNum;i++){
		file<<"3 "<<mesh->faceList[i].index[1]<<' '<<mesh->faceList[i].index[2]<<' '<<mesh->faceList[i].index[3]<<endl;
	}
	file.close();
}

void readOFF(string filename,double* u){
	ifstream fin;
	fin.open(filename.c_str());

	string newLine;
	std::getline(fin,newLine);

	int extra;
	int vn,fn;
	float* x;
	float* y;
	float* z;
	int* fa,*fb,*fc;

	fin>>vn>>fn>>extra;

	x = new float[vn];
	y = new float[vn];
	z = new float[vn];

	for (int i = 0; i < vn; i++){
		fin>>x[i]>>y[i]>>z[i];
		u[i*3] = x[i];
		u[i*3+1] = y[i];
		u[i*3+2] = z[i];
	}
	fa = new int[fn];
	fb = new int[fn];
	fc = new int[fn];

	for (int i = 0; i < fn; i++)
		fin>>extra>>fa[i]>>fb[i]>>fc[i];

	fin.close();
}

Vector_3 projectVertex(Vector_3 d1,Vector_3 d2,Vector_3 edge){
	Vector_3 normal = cross_product(d1,d2);
	normal = normal / sqrt(normal.squared_length());

	float signN = edge*normal;
	Vector_3 inNormal = (edge*normal) * normal; 
	Vector_3 projected = edge - inNormal;

	projected = projected / sqrt(projected.squared_length());
	return projected;
}

Vector_3 rotate(Vector_3 dir,Vector_3 axis,double rotation){
	float nx,ny,nz;
	float x,y,z;

	x = axis.x();
	y = axis.y();
	z = axis.z();

	float ox,oy,oz;
	ox = dir.x();
	oy = dir.y();
	oz = dir.z();

	nx = (x*x+(1-x*x)*cos(rotation))*ox + (x*y*(1-cos(rotation))-z*sin(rotation))*oy + (x*z*(1-cos(rotation))+y*sin(rotation))*oz;
	ny = (y*x*(1-cos(rotation))+z*sin(rotation))*ox + (y*y+(1-y*y)*cos(rotation))*oy + (y*z*(1-cos(rotation))-x*sin(rotation))*oz;
	nz = (z*x*(1-cos(rotation))-y*sin(rotation))*ox + (z*y*(1-cos(rotation))+x*sin(rotation))*oy + (z*z+(1-z*z)*cos(rotation))*oz;

	return Vector_3(nx,ny,nz);
}

double computeDeterminant(Vertex_const_handle v1,Vertex_const_handle v2,Vertex_const_handle v3,Vertex_const_handle v4){
	double det234 = v2->point().x() * v3->point().y() * v4->point().z() + 
					v2->point().y() * v3->point().z() * v4->point().x() +
					v2->point().z() * v3->point().x() * v4->point().y() -
					v2->point().z() * v3->point().y() * v4->point().x() -
					v2->point().y() * v3->point().x() * v4->point().z() -
					v2->point().x() * v3->point().z() * v4->point().y();

	double det134 = 1 * v3->point().y() * v4->point().z() + 
					v2->point().y() * v3->point().z() * 1 +
					v2->point().z() * 1 * v4->point().y() -
					v2->point().z() * v3->point().y() * 1 -
					v2->point().y() * 1 * v4->point().z() -
					1 * v3->point().z() * v4->point().y();

	double det124 = 1 * v3->point().x() * v4->point().z() +
					v2->point().x() * v3->point().z() * 1 + 
					v2->point().z() * 1 * v4->point().x() -
					v2->point().z() * v3->point().x() * 1 -
					v2->point().x() * 1 * v4->point().z() -
					1 * v3->point().z() * v4->point().x();

	double det123 = v2->point().x() * v3->point().y() * 1 + 
					v2->point().y() * 1 * v4->point().x() +
					1 * v3->point().x() * v4->point().y() -
					1 * v3->point().y() * v4->point().x() -
					v2->point().y() * v3->point().x() * 1 -
					v2->point().x() * 1 * v4->point().y();

	return 1*det234 - v1->point().x()*det134 +
					  v1->point().y()*det124 -
					  v1->point().z()*det123;


}

double* computeDeterminant(Point_3 v1,Point_3 v2,Point_3 v3,Point_3 v4){
	double* det = new double[5];
	
	double det234 = v2.x() * v3.y() * v4.z() + 
					v2.y() * v3.z() * v4.x() +
					v2.z() * v3.x() * v4.y() -
					v2.z() * v3.y() * v4.x() -
					v2.y() * v3.x() * v4.z() -
					v2.x() * v3.z() * v4.y();

	double det134 = 1 * v3.y() * v4.z() + 
					v2.y() * v3.z() * 1 +
					v2.z() * 1 * v4.y() -
					v2.z() * v3.y() * 1 -
					v2.y() * 1 * v4.z() -
					1 * v3.z() * v4.y();

	double det124 = 1 * v3.x() * v4.z() +
					v2.x() * v3.z() * 1 + 
					v2.z() * 1 * v4.x() -
					v2.z() * v3.x() * 1 -
					v2.x() * 1 * v4.z() -
					1 * v3.z() * v4.x();

	double det123 = v2.x() * v3.y() * 1 + 
					v2.y() * 1 * v4.x() +
					1 * v3.x() * v4.y() -
					1 * v3.y() * v4.x() -
					v2.y() * v3.x() * 1 -
					v2.x() * 1 * v4.y();
	
	det[1] = det234;
	det[2] = det134;
	det[3] = det124;
	det[4] = det123;
	det[0] = 1*det[1]- v1.x()*det[2] + v1.y()*det[3] - v1.z()*det[4];

	return det;
}

float computeEuclideanDis(Point_3 a,Point_3 b){
	float dis = sqrt((a.x()-b.x())*(a.x()-b.x()) + (a.y()-b.y())*(a.y()-b.y()) + (a.z()-b.z())*(a.z()-b.z()));
	return dis;
}

float computeEuclideanDis(threeTuple a ,threeTuple b){
	float dis = sqrt((a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y) + (a.z-b.z)*(a.z-b.z));
	return dis;
}

float computeAngle(Point_3 V1,Point_3 V2,Point_3 V3,Point_3 V4){
	double theta1;

	Vector_3 b1(V1,V2);
	Vector_3 b2(V2,V3);
	Vector_3 s1(V1,V4);

	Vector_3 nm = cross_product(b1,b2);
	Vector_3 n1 = cross_product(s1,b1);

	double sign1 = s1*nm;
	double x1 = n1*nm / (sqrt(n1.squared_length()+EPS)*sqrt(nm.squared_length()+EPS));
	if (sign1 > 0) x1 = (2 - x1);

	if (x1 <= 1) theta1 = (-0.69813170079773212 * x1 * x1 - 0.87266462599716477) * x1 + 1.5707963267948966;
	else	     theta1 = (-0.69813170079773212 * (x1-2) * (x1-2) - 0.87266462599716477) * (x1-2) - 1.5707963267948966;

	return theta1;
}

double computeAngle(Vector_3 v1,Vector_3 v2){
	double costheta = v1*v2/(sqrt(v1.squared_length())*sqrt(v2.squared_length()));
	return acos(costheta);
}

Vector_3 closesPointOnTriangle( const Vector_3 *triangle, const Vector_3 &sourcePosition )
{
	Vector_3 edge0 = triangle[1] - triangle[0];
	Vector_3 edge1 = triangle[2] - triangle[0];
	Vector_3 v0 = triangle[0] - sourcePosition;

	float a = edge0 * edge0;
	float b = edge0 * edge1;
	float c = edge1 * edge1;
	float d = edge0 * v0;
	float e = edge1 * v0;

	float det = a*c - b*b;
	float s   = b*e - c*d;
	float t   = b*d - a*e;

	if ( s + t < det )
	{
		if ( s < 0.f )
		{
			if ( t < 0.f )
			{
				if ( d < 0.f )
				{
					s = clamp_zenyo( -d/a, 0.f, 1.f );
					t = 0.f;
				}
				else
				{
					s = 0.f;
					t = clamp_zenyo( -e/c, 0.f, 1.f );
				}
			}
			else
			{
				s = 0.f;
				t = clamp_zenyo( -e/c, 0.f, 1.f );
			}
		}
		else if ( t < 0.f )
		{
			s = clamp_zenyo( -d/a, 0.f, 1.f );
			t = 0.f;
		}
		else
		{
			float invDet = 1.f / det;
			s *= invDet;
			t *= invDet;
		}
	}
	else
	{
		if ( s < 0.f )
		{
			float tmp0 = b+d;
			float tmp1 = c+e;
			if ( tmp1 > tmp0 )
			{
				float numer = tmp1 - tmp0;
				float denom = a-2*b+c;
				s = clamp_zenyo( numer/denom, 0.f, 1.f );
				t = 1-s;
			}
			else
			{
				t = clamp_zenyo( -e/c, 0.f, 1.f );
				s = 0.f;
			}
		}
		else if ( t < 0.f )
		{
			if ( a+d > b+e )
			{
				float numer = c+e-b-d;
				float denom = a-2*b+c;
				s = clamp_zenyo( numer/denom, 0.f, 1.f );
				t = 1-s;
			}
			else
			{
				s = clamp_zenyo( -e/c, 0.f, 1.f );
				t = 0.f;
			}
		}
		else
		{
			float numer = c+e-b-d;
			float denom = a-2*b+c;
			s = clamp_zenyo( numer/denom, 0.f, 1.f );
			t = 1.f - s;
		}
	}

	return triangle[0] + s * edge0 + t * edge1;
}
