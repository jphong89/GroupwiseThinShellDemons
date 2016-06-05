#pragma once
#include "MeshObject.h"
#include <gl\glut.h>
#include "pthread.h"

extern BasicMesh** Surface;
extern int surfaceNum;
extern string* fileName;

extern float camera_Scale;
extern float camera_Up;
extern float camera_Right;

extern float rotate_X;
extern float rotate_Y;
extern float rotate_Z;

extern bool switch1;
extern bool switch2;
extern bool switch3;
extern bool switch4;

/*debug*/
extern double facetBending;
extern double facetStretching;
extern double distantLink;
extern LARGE_INTEGER t1;
extern LARGE_INTEGER t2;
extern LARGE_INTEGER tc;
extern double thetaDif1;
extern double thetaDif2;
extern double thetaDif3;
extern float geodesicDistance;
extern pthread_mutex_t lock;

// void drawSurface(BasicMesh* mesh,int deform);
// void drawExtra(BasicMesh*);
// void display();
// void reshape(int w,int h);
// void init();
// int beginGlut(int argc, char* argv[]);
// void ProcessSpecialKeyboead(int key, int x, int y);
// void ProcessKeyboard(unsigned char key,int x,int y);
// void Mouse(int button, int state, int x, int y);
// void paintColor(BasicMesh*,threeTuple,threeTuple,threeTuple);
// PolyhedralSurf::Vertex_const_handle mouseEvent(BasicMesh*,threeTuple);
// void drawSphere(GLfloat xx, GLfloat yy, GLfloat zz, GLfloat radius, GLfloat M, GLfloat N);

//CPU
double evaluateError(double bestError);
void startOptimization();
void initialDVF();

static int progress(
	void *instance,
	const lbfgsfloatval_t *u,
	const lbfgsfloatval_t *g,
	const lbfgsfloatval_t fx,
	const lbfgsfloatval_t xnorm,
	const lbfgsfloatval_t gnorm,
	const lbfgsfloatval_t step,
	int n,
	int k,
	int ls
	);

static lbfgsfloatval_t evaluate(
	void *instance,
	const lbfgsfloatval_t *u,
	lbfgsfloatval_t *g,
	const int n,
	const lbfgsfloatval_t step
	);

lbfgsfloatval_t penalizeData(const lbfgsfloatval_t *u, lbfgsfloatval_t *g, int idx);
lbfgsfloatval_t penalizeStretch(const lbfgsfloatval_t *u, lbfgsfloatval_t *g, int idx);
lbfgsfloatval_t penalizeLink(const lbfgsfloatval_t *u, lbfgsfloatval_t *g, int idx);
lbfgsfloatval_t penalizeBendQuadratic(const lbfgsfloatval_t *tu, lbfgsfloatval_t *g, int idx);
lbfgsfloatval_t penalizeBendAngleLP(const lbfgsfloatval_t *u, lbfgsfloatval_t *g, int idx);
lbfgsfloatval_t penalizeBend(const lbfgsfloatval_t *u, lbfgsfloatval_t *g, int idx);
lbfgsfloatval_t penalizeLandmark(const lbfgsfloatval_t *u, lbfgsfloatval_t *g, int idx);
lbfgsfloatval_t penalizeStretchQuadratic(const lbfgsfloatval_t *u, lbfgsfloatval_t *g, int surfaceIdx);

double computeStretchGradient(int n,int nd,int e,int ed);
double computeAngleGradient(int i,int j,Point_3 V1,Point_3 V2,Point_3 V3,Point_3 V4);
double computeGradient(int n,int nd,int e,int ed);
double computeDetGradient(int i,int j,Point_3 v1,Point_3 v2,Point_3 v3,Point_3 v4,double* det);
inline threeTuple computeGradient(int v,int i,int j);
inline threeTuple computeCrossProductGradient(Vector_3 b1,Vector_3 b2,threeTuple db1,threeTuple db2);
inline double computeDotProductGradient(Vector_3 n1,Vector_3 n2,threeTuple dn1,threeTuple dn2);
