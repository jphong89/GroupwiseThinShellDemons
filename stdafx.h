// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

//#pragma once

#ifndef __STDAFX_H__
#define __STDAFX_H__

#include <stdio.h>
#include <cstdlib>
#define VAL(SA,ia,ja,na) SA[(ia)*(na)+(ja)]
const double PI = 3.141592654;
const double EPS = 1e-6;
#define max_zenyo(a,b)            (((a) > (b)) ? (a) : (b))
#define min_zenyo(a,b)            (((a) < (b)) ? (a) : (b))
#define clamp_zenyo(x,y,z)              (min_zenyo(max_zenyo((x),(y)),(z)))

#define DIR_NUM 8
#define GEODECISDIS 4.0 // inactive: length of the "spider leg"
#define STRETCH_MIU (YOUNG/2/(1+POISSON)) //inactive 
#define STRETCH_LAMBDA (YOUNG*POISSON/(1+POISSON)/(1-2*POISSON)) //inactive
#define LINKWEIGHT 200  //structural links weight


#define GPU_ON 0

extern double YOUNG;
extern double POISSON;
extern double REGWEIGHT; //stretching weight
extern double BENDWEIGHT;  //bending weight
extern double DISTHRESHOLD; // we don't want to match things too far away
extern double EUCLWEIGHT;
extern int MESHLABOPTION;
extern double LANDMARKWEIGHT;
extern int TEXTURE_RANGE;
extern int GEOMETRIC_RANGE;

// TODO: reference additional headers your program requires here
#endif // __STDAFX_H__
