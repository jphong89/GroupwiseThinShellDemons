// Demons.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "MeshObject.h"
#include "GlutManager.h"
#include <time.h>

#include <iostream>
#include <fstream>
using namespace std;

/* registration params */
double YOUNG = 2;
double POISSON = 0.05;
double REGWEIGHT = 40; //stretching weight
double BENDWEIGHT = 400;  //bending weight
double DISTHRESHOLD = 12.0; // we don't want to match things too far away
double EUCLWEIGHT = 1;
int MESHLABOPTION = 1;
double LANDMARKWEIGHT = 1;
int TEXTURE_RANGE = 2;
int GEOMETRIC_RANGE = 20; // no more than 20 surfaces allowed

float geodesicDistance = GEODECISDIS;

/* surface lists */
BasicMesh** Surface;
int surfaceNum;
string* fileName;
string* landmarkFile;

/* registration variables */
lbfgsfloatval_t* u;
double bendweight;
double regweight;

/*debug*/
double facetBending = 0;
double facetStretching = 0;
double distantLink = 0;
double thetaDif1,thetaDif2,thetaDif3;
double facetTrace,facetDet;

/* camera variables */
float camera_Scale = 1;
float camera_Up = 0;
float camera_Right = 0;

float rotate_X = 0;
float rotate_Y = 0;
float rotate_Z = 0;

bool switch1 = false;
bool switch2 = false;
bool switch3 = false;
bool switch4 = false;

float minx,miny,minz;
float maxx,maxy,maxz;
float max_all;

//pthread_mutex_t lock;

void readFileList(string list){
    ifstream fin;

    fin.open(list.c_str());

    fileName = new string[surfaceNum];

    for (int i = 0; i < surfaceNum; i++)
        fin>>fileName[i];

    for (int i = 0; i < surfaceNum; i++){
        cout<<i<<" file name: "<<fileName[i];
        cout<<endl;
    }
}

void readMatchingList(string list){
    ifstream fin;

    fin.open(list.c_str());

    landmarkFile = new string[surfaceNum];

    for (int i = 0; i < surfaceNum; i++)
        fin>>landmarkFile[i];

    for (int i = 0; i < surfaceNum; i++){
        cout<<i<<" landmark name: "<<landmarkFile[i];
        cout<<endl;
    }
}

int main(int argc, char* argv[])
{
    char* landmark1; //landmark1
    char* landmark2; //landmark2

    //format: demons.exe surface1.off surface2.off weight_output.txt landmark1.txt landmark2.txt
    if (argc > 2){
        surfaceNum = atoi(argv[1]);

        string fileList;
        fileList = argv[2];
        readFileList(fileList);

        if (argc > 3) {
            REGWEIGHT = atof(argv[3]);
            BENDWEIGHT = atof(argv[4]);
        }

        if (argc > 5) DISTHRESHOLD = atof(argv[5]);
        if (argc > 6) EUCLWEIGHT = atof(argv[6]);

        if (argc > 7) {
            YOUNG = atof(argv[7]);
            POISSON = atof(argv[8]);
        }

        if (argc > 9) {
            MESHLABOPTION = atoi(argv[9]);
        }

        string matchingFile;
        if (argc > 10) {
            matchingFile = argv[10];
            readMatchingList(matchingFile);
            TEXTURE_RANGE = atoi(argv[11]);
        }
        cout<<"STRETCHING WEIGHT: "<<REGWEIGHT<<" -BENDING WEIGHT: "<<BENDWEIGHT<<" -EUCLIDEAN THRESHOLD: "
            <<DISTHRESHOLD<<" -GEOMETRIC FEATURE WEIGHT: "<<EUCLWEIGHT
            <<" -YOUNG: "<<YOUNG<<" -POSSION: "<<POISSON<<" -MESHLAB: "<<MESHLABOPTION<<endl;
    }else{
        cout<<"not enough parameters, refer ReadMe.txt";
    }

    omp_set_dynamic(0);
    omp_set_num_threads(8);
    int i,j;
    ////////////////////////////////////*read & translate surfaces */////////////////////////////////////////////
    Surface = new BasicMesh*[surfaceNum];

#pragma omp parallel for private(i)
    for (i = 0; i < surfaceNum; i++){
        cout<<"Begin Constructing Mesh: "<<fileName[i]<<endl;
        Surface[i] = new BasicMesh(fileName[i],i);
        Surface[i]->ComputeMeshProperty();
        Surface[i]->findSignatureAll();
        Surface[i]->constructEdgeList();
        Surface[i]->constructVertexList();
        Surface[i]->initialDeformation();
        if (argc > 10) Surface[i]->constructLandmark(landmarkFile[i]);

        //Surface[i]->outputFeature();
    }
    /////////////////////////////* 3compute correspondences for registration *///////////////////////////////////
    for (i = 0; i < surfaceNum; i++){
        cout<<"Computing attraction force: "<<fileName[i]<<endl;

        int startSurface = max_zenyo(i - GEOMETRIC_RANGE, 0);
        int endSurface = min_zenyo(i + GEOMETRIC_RANGE, surfaceNum - 1);
        for (int j = startSurface; j <= endSurface; j++)
            if (i != j) {
                Surface[i]->findCorrespondenceBothWay(Surface[j],EUCLWEIGHT);
                Surface[i]->findMatch2(Surface[j],j);
                /*Surface[i]->outputAffinity(Surface[j],j);*/
            }

        Surface[i]->summerizeForce();
        //Surface[i]->outputForce();
    }
    /////////////////////////////* rendering and optimization *////////////////////////////////////////////////////////////
    startOptimization(); //direct begin for script
    return 0;
}

