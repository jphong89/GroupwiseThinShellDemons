#pragma once

#include <CGAL/Cartesian.h>
#include <CGAL/Ridges.h>
#include <CGAL/Umbilics.h>
#include <CGAL/Monge_via_jet_fitting.h>
#include <fstream>
#include <cassert>
#include <set>
#include "omp.h"
//this is an enriched Polyhedron with facets' normal
#include "PolyhedralSurf.h"
#include "PolyhedralSurf_rings.h"
#include "lbfgs.h"

#include <Eigen/Geometry>
#include <Eigen/LU>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

/* code borrowed from CGAL */
typedef PolyhedralSurf::Traits          Kernel;
typedef Kernel::FT                      FT;
typedef Kernel::Point_3                 Point_3;
typedef Kernel::Vector_3                Vector_3;
typedef Kernel::Vector_2                Vector_2;

typedef PolyhedralSurf::Vertex_const_handle   Vertex_const_handle;
typedef PolyhedralSurf::Vertex_const_iterator Vertex_const_iterator;
typedef PolyhedralSurf::Facet_const_handle Facet_const_handle;

typedef T_PolyhedralSurf_rings<PolyhedralSurf> Poly_rings;
typedef CGAL::Monge_via_jet_fitting<Kernel>    Monge_via_jet_fitting;
typedef Monge_via_jet_fitting::Monge_form      Monge_form;

typedef CGAL::Vertex2Data_Property_Map_with_std_map<PolyhedralSurf> Vertex2Data_Property_Map_with_std_map;
typedef Vertex2Data_Property_Map_with_std_map::Vertex2FT_map Vertex2FT_map;
typedef Vertex2Data_Property_Map_with_std_map::Vertex2Vector_map Vertex2Vector_map;
typedef Vertex2Data_Property_Map_with_std_map::Vertex2FT_property_map Vertex2FT_property_map;
typedef Vertex2Data_Property_Map_with_std_map::Vertex2Vector_property_map Vertex2Vector_property_map;

//RIDGES
typedef CGAL::Ridge_line<PolyhedralSurf> Ridge_line;
typedef CGAL::Ridge_approximation < PolyhedralSurf,
        Vertex2FT_property_map,
        Vertex2Vector_property_map > Ridge_approximation;
//UMBILICS
typedef CGAL::Umbilic<PolyhedralSurf> Umbilic;
typedef CGAL::Umbilic_approximation < PolyhedralSurf,
        Vertex2FT_property_map,
        Vertex2Vector_property_map > Umbilic_approximation;

/* MY OWN DATA STRUCTURE, WHICH IS QUITE USELESS NOW */
struct threeTuple 
{
    float x;
    float y;
    float z;

    threeTuple(){
        x = 0;
        y = 0;
        z = 0;
    }

    threeTuple(float a,float b,float c){
        x = a;
        y = b;
        z = c;
    }

    threeTuple(Point_3 input){
        x = input.x();
        y = input.y();
        z = input.z();
        next = NULL;
    }

    threeTuple(Vector_3 input){
        x = input.x();
        y = input.y();
        z = input.z();
        next = NULL;
    }
    threeTuple* next;
};

struct facet{
    threeTuple p1,p2,p3; // VERTEX INFORMATION

    threeTuple color1,color2,color3; // VERTEX COLOR

    threeTuple d1,d2; //PRINCIPAL DIRECTION

    Vertex_const_handle v1,v2,v3;
    Vertex_const_handle nb1,nb2,nb3; // 3 NEIGHBOUR VERTICES
    double theta1,theta2,theta3;	// 3 ANGLES (BETWEEN NEIGHBOUR FACETS)
    double sideArea1,sideArea2,sideArea3;	//AREAS OF NEIGHBOUR FACETS
    double l1,l2,l3; 

    double SO1[2][2]; //  
    double SO2[2][2]; //
    double SO3[2][2]; //

    Eigen::Matrix3d basis; 
    int index[4];
    int nbIdx1,nbIdx2,nbIdx3;
    int faceIdx1,faceIdx2,faceIdx3;

    Vector_2 local_v1,local_v2,local_v3;
    Vector_2 t1,t2,t3;

    double inverse[2][2]; //INVERSE OF STRETCHING TENSOR
    double area; //FACET AREA

    bool isBorder; // BORDER VERTEX?

    threeTuple normal;
};

struct edge{
    Vertex_const_handle v1,v2; //from v1 to v2
    int index1,index2;

    Vertex_const_handle vl,vr;
    int indexl,indexr;

    int index[5];

    PolyhedralSurf::Halfedge_const_handle eh;

    double angle;

    double stretchWeight;
    double bendWeight;

    bool isBoundary;
};

struct link{
    int index1,index2;

    double length;
};

struct signature{
    float k1,k2;
    float C,S;
    float K,H;

    Vector_3 normal;
    Vector_3 point;
    Vector_3 d1,d2;

    Eigen::Quaternionf rotation;

    float deltaS,deltaN;

    float deltaN12,deltaN34;
    float dis12,dis34,dis56,dis78,dis58,dis67;

    threeTuple texture;

    signature(){
        k1 = 0;
        k2 = 0;
        C = 0;
        S = 0;
        K = 0;
        H = 0;

        deltaS = 0;
        deltaN = 0;

        deltaN12 = 0; //normal between v1 v2
        deltaN34 = 0; //normal between v3 v4

        dis12 = 0;	//euc dis between v1 v2
        dis34 = 0;
        dis56 = 0;
        dis78 = 0;
        dis58 = 0;
        dis67 = 0;

        texture.x = 0;
        texture.y = 0;
        texture.z = 0;

        normal = Vector_3(1,0,0);
        point = Vector_3(0,0,0);
        d1 = Vector_3(0,1,0);
        d2 = Vector_3(0,0,1);

        rotation = Eigen::Quaternionf(0,0,0,0);
    }

    const signature& operator=( const signature& other )
    {
        if ( this == &other )
        {
            return *this;
        }
        k1 = other.k1;
        k2 = other.k2;
        C = other.C;
        S = other.S;
        K = other.K;
        H = other.H;
        point = other.point;
        normal = other.normal;
        d1 = other.d1;
        d2 = other.d2;
        rotation = other.rotation;

        deltaS = other.deltaS;
        deltaN = other.deltaN;

        texture = other.texture;

        return other;
    }
};

struct node{
    Vertex_const_handle vh;
    double dis;

    friend bool operator < (node a, node b)
    {
        return a.dis > b.dis;
    }

    friend bool operator > (node a, node b)
    {
        return a.dis < b.dis;
    }

    node(){}

    node(Vertex_const_handle v,double d){
        vh = v;
        dis = d;
    }
};

struct vertex{
    int idx;
    int neighbourNum;
    //Point_3* pt;
    double* weight;
    int* nbIdx;
    double area;
};

/* node in a queue , for qsort purpose */
struct qnode{
    float value;
    int i,j;
};

class BasicMesh
{
    public:
        // cgal required para's
        Vertex2FT_map vertex2k1_map, vertex2k2_map,
                      vertex2b0_map, vertex2b3_map,
                      vertex2P1_map, vertex2P2_map;
        Vertex2Vector_map vertex2d1_map, vertex2d2_map;

        // default fct parameter values and global variables
        unsigned int d_fitting;
        unsigned int d_monge;

        unsigned int nb_rings;//seek min # of rings to get the required #pts
        unsigned int nb_points_to_use;//

        CGAL::Ridge_order tag_order;
        unsigned int int_tag;

        double umb_size;
        unsigned int min_nb_points;

        /* data information */
        string if_name;
        int idx;

        Vertex_const_handle centre;
        Vertex_const_handle corresponding;

    public:
        PolyhedralSurf P; //the surface, core structure

    public:
        /* geometric feature construction*/
        int ComputeMeshProperty();
        void translateMesh(Vertex2FT_property_map vertex2k1_pm,Vertex2FT_property_map vertex2k2_pm); // compute facet information
        void ComputeStretchInverse(Vertex_const_handle v1,Vertex_const_handle v2,Vertex_const_handle v3,facet* face);
        void ComputeShapeOperator(facet* face);

        Vertex_const_handle marchOnSurface(PolyhedralSurf::Halfedge_const_handle,PolyhedralSurf::Facet_const_handle,Vector_3,float,float,float,int,threeTuple*); // geodesic marching

        void gather_fitting_points(Vertex_const_handle v,std::vector<Point_3> &in_points,Poly_rings& poly_rings);
        void compute_differential_quantities(PolyhedralSurf& P, Poly_rings& poly_rings);
        void deleteMarchList();

        signature* findSignature(int dirOpt,Vertex_const_handle vh,double rotation = 0);
        signature* constructSignature(Vertex_const_handle vh,double rotation = 0);
        Vertex_const_handle findVertexHandle(threeTuple a);
        void findSignatureAll();

    public:
        BasicMesh(string,int);
        BasicMesh(){}

        /* Prepare registration data structure */
        void constructVertexList();
        void constructEdgeList();
        void constructLink();
        void initialDeformation();
        void summerizeForce();
        void constructLandmark(string);

        /*REGISTRATION matching pipeline*/
        void findCorrespondenceBothWay(BasicMesh*,double); //compute feature similarity based on euc dis & curvatures3
        void findMatch2(BasicMesh* secondMesh, int surfaceIdx); //from static to moving
        void outputAffinity(BasicMesh* secondMesh, int surfaceIdx); // for test purpose
        void outputForce();
        void outputFeature();

        /* surface information */
        int faceNum;
        int vertexNum;
        int edgeNum;
        int dirNum;
        int *landmarkNum;
        int **landmarks1;
        int **landmarks2;

        /* registration data structure */
        Vertex_const_handle matched;
        facet* faceList;
        edge* edgeList;
        vertex* vertexList;

        Vector_3** bestMatch;
        float** matchWeight;
        double* u;
        Vector_3* targetPos;
        double* overallWeight;
        float* affinity;

        /* supporting data structure */
        threeTuple* marchList[8];
        std::map<Vertex_const_handle, signature*> signatureMap; 
        std::map<Vertex_const_handle, int> indexMap;
        std::map<Vertex_const_handle, float> meshColor;
        Vertex_const_handle *vertexIndex;
        std::set<PolyhedralSurf::Facet_const_handle> vertexBool;
};

/* helper functions */
Vector_3 projectVertex(Vector_3 d1,Vector_3 d2,Vector_3 edge);
float compareSignature(signature* sig1,signature* sig2);
Vector_3 rotate(Vector_3 dir,Vector_3 axis,double rotation);
threeTuple getJetColor(float);
float computeEuclideanDis(Point_3,Point_3);
float computeEuclideanDis(threeTuple,threeTuple);
void rotateSig(signature*);
Eigen::Quaternionf constructQuaternion(signature*,signature*);
void writeOFF(BasicMesh* mesh, string filename, int idx);
void readOFF(string filename,double* u);
float computeArea(PolyhedralSurf::Facet_const_handle f);
float computeLength(PolyhedralSurf::Halfedge_const_handle e);
double computeAngle(Vector_3 v1,Vector_3 v2);
double computeDeterminant(Vertex_const_handle v1,Vertex_const_handle v2,Vertex_const_handle v3,Vertex_const_handle v4);
double* computeDeterminant(Point_3 v1,Point_3 v2,Point_3 v3,Point_3 v4);
float computeAngle(Point_3 V1,Point_3 V2,Point_3 V3,Point_3 V4);
double computeAngle(Vector_3 v1,Vector_3 v2);
Vector_3 closesPointOnTriangle( const Vector_3 *triangle, const Vector_3 &sourcePosition );

extern float geodesicDistance;

extern BasicMesh** Surface;
extern int surfaceNum;

extern float minx,miny,minz;
extern float maxx,maxy,maxz;
extern float max_all;
