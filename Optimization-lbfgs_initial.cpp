#include "stdafx.h"
#include "GlutManager.h"
#include <omp.h>
#include <iostream>
#include <fstream>
using namespace std;

double* tempG;

static int progress(
	void *instance,
	const lbfgsfloatval_t *tu,
	const lbfgsfloatval_t *g,
	const lbfgsfloatval_t fx,
	const lbfgsfloatval_t xnorm,
	const lbfgsfloatval_t gnorm,
	const lbfgsfloatval_t step,
	int n,
	int k,
	int ls
	)
{
	pthread_mutex_lock(&lock);
	if ((k % 100 == 0) || (k <=1)){
		memcpy(u,tu,sizeof(lbfgsfloatval_t)*affinityM*3);

		cout<<"Iteration: "<<k<<" -------------------------------"<<endl;
		cout<<"fx = "<<fx<<endl;
		cout<<"xnorm = "<<xnorm<<", gnorm = "<<gnorm<<", step = "<<step<<endl;
		cout<<endl;
		cout<<"Bending: "<<facetBending<<" Stretching: "<<facetStretching<<" Link: "<<distantLink<<endl;
		cout<<thetaDif1<<' '<<thetaDif2<<' '<<thetaDif3<<endl;

		QueryPerformanceFrequency(&tc);
		QueryPerformanceCounter(&t2);
		cout<<"Time Elapse:"<<(t2.QuadPart - t1.QuadPart)*1.0/tc.QuadPart * 1000<<endl;
		cout<<endl;
	}
	pthread_mutex_unlock(&lock);
	return 0;
}

void initialDVF(double* initialU){	//
	char filename2[50] = "temp8_occluded2.off";

	double* newConfig = new double[affinityM*3];
	readOFF(filename2,newConfig);
	
	PolyhedralSurf::Vertex_const_iterator vb = RCSurface->P.vertices_begin();
	PolyhedralSurf::Vertex_const_iterator ve = RCSurface->P.vertices_end();
	PolyhedralSurf::Vertex_const_handle vh;

	for (int i = 0; vb != ve; vb++,i++){
		vh = vb;

		initialU[i*3] = newConfig[i*3] - vh->point().x();
		initialU[i*3+1] = newConfig[i*3+1] - vh->point().y();
		initialU[i*3+2] = newConfig[i*3+2] - vh->point().z();
	}
}

double evaluateError(double bestError){
	/* output errer */
	float finalErr,finalVar,finalBoundErr,finalBoundVar,maxErr;
	computeStatisticBool(RCSurface,CTSurface,finalErr,finalBoundErr,finalVar,finalBoundVar,maxErr);

	if (finalBoundErr < bestError){
		bestError = finalBoundErr;

		ofstream fout;
		fout.open("error.txt");

		fout<<initialError<<' '<<initialVariance<<endl;
		fout<<initialBoundError<<' '<<initialBoundVar<<endl;


		fout<<finalErr<<' '<<finalVar<<endl;
		fout<<finalBoundErr<<' '<<finalBoundVar<<endl;

		fout.close();
	}

	return bestError;
	/* output errer */
}

void *startOptimization(void * arg){
	/* initialize optimization */
	char filename[50];

	int ret;
	lbfgsfloatval_t *tempU = new lbfgsfloatval_t[affinityM*3];
	tempG = new double[RCSurface->faceNum*9];
	memset(tempG,0,sizeof(double)*RCSurface->faceNum*9);
	memset(tempU,0,sizeof(lbfgsfloatval_t) * affinityM * 3);

	//omp_set_dynamic(0);
	//omp_set_num_threads(8);

	double bestError = 100000;
	regweight = REGWEIGHT;
	bendweight = BENDWEIGHT;

	for (int i = 0; i < 10; i++){
		cout<<"Velocity Iteration: "<<i<<endl;

		if (i > 0){
			cout<<"Reset surface..."<<endl;

			RCSurface->ComputeMeshProperty(filename);
			bestError = evaluateError(bestError);

			cout<<"Recompute geometric features..."<<endl;
			RCSurface->findSignatureAll();

			memset(tempU,0,sizeof(lbfgsfloatval_t) * affinityM * 3);
			memset(affinity,0,affinityM * affinityN * sizeof(double));
			memset(bestMatch,0,sizeof(int)*affinityM);
			memset(matchWeight,0,sizeof(float)*affinityM);
			
			cout<<"Recompute affinity map..."<<endl;

			/* prepare registration */
			RCSurface->constructEdgeList();
			//RCSurface->constructLink();
			RCSurface->constructVertexList();
			RCSurface->findCorrespondenceBothWay(CTSurface,EUCLWEIGHT);
			cout<<"Recompute affinity map..."<<endl;
			RCSurface->findMatch2(CTSurface);
			//RCSurface->findClosestPoint(CTSurface);
		}

		cout<<"Begin optimization..."<<endl;    /* actual registration (optimization) */

		int opIterNum = 80000;

		lbfgs_parameter_t param;
		ret = 0;

		clock_t start, finish;
		double duration;
		start = clock();
		QueryPerformanceCounter(&t1);

		lbfgsfloatval_t fx = 0;
		lbfgs_parameter_init(&param);
		param.max_iterations = opIterNum;
		param.epsilon = 1e-5f;
		ret = lbfgs(affinityM * 3, tempU, &fx, evaluate, progress, NULL, &param);

		finish = clock();
		duration = (double)(finish - start) / CLOCKS_PER_SEC;
		cout<<"function value: "<<fx<<"  Elapsed Time: "<<duration<<" s"<<endl;

		memcpy(u,tempU,sizeof(lbfgsfloatval_t)*affinityM*3);
		//output file
		filename[0] = '\0';
		strcat(filename,"temp");
		char iterNum[10];
		itoa(i,iterNum,10);
		strcat(filename,iterNum);
		strcat(filename,".off");
		writeOFF(RCSurface,filename);

		cout<<"L-BFGS optimization terminated with status code: "<<ret<<"fx: "<<fx<<endl;

		if ((MESHLABOPTION > 0) && (i % MESHLABOPTION == 0)) {
			char cmdLine[100];
			sprintf(cmdLine,"meshlabserver -i %s -o %s -s meshlabscript_demons.mlx\n",filename,filename);
			system(cmdLine);
		}
		
		regweight *= 0.9;
		bendweight *= 0.9;
	}

	//final iteration evaluation
	RCSurface->ComputeMeshProperty(filename);
	bestError = evaluateError(bestError);

	return NULL;
}

static lbfgsfloatval_t evaluate(
	void *instance,
	const lbfgsfloatval_t *u,
	lbfgsfloatval_t *g,
	const int n,
	const lbfgsfloatval_t step
	)
{
	pthread_mutex_lock(&lock);
	//cout<<"begin optimization"<<endl;
	double data = penalizeData(u,g);

	/* stretching */
	double stretching = penalizeStretch(u,g);
	//double stretching = 0;
	/* bending */
	double bending = penalizeBendQuadratic(u,g);
	//double bending = 0;
	double landmark = penalizaeLandmark(u,g);

 	//double linkStretching = penalizeLink(u,g);
// 	double linkStretching = 0;
// 	for (int i = 0; i < RCSurface->linkNum; i++){
// 		int idx1 = linkList[i].index1;
// 		int idx2 = linkList[i].index2;
// 
// 		Vector_3 deformVec(u[idx2*3] - u[idx1*3],u[idx2*3+1] - u[idx1*3+1],u[idx2*3+2] - u[idx1*3+2]);
// 
// 		g[idx1*3] += -2 * LINKWEIGHT * deformVec.x();
// 		g[idx1*3+1] += -2 * LINKWEIGHT * deformVec.y();
// 		g[idx1*3+2] += -2 * LINKWEIGHT * deformVec.z();
// 		g[idx2*3] += 2 * LINKWEIGHT * deformVec.x();
// 		g[idx2*3+1] += 2 * LINKWEIGHT * deformVec.y();
// 		g[idx2*3+2] += 2 * LINKWEIGHT * deformVec.z();
// 
// 		linkStretching += LINKWEIGHT * deformVec.squared_length();
// 	}
	pthread_mutex_unlock(&lock);
	return data + stretching +  bending +landmark /*+ linkStretching*/;
}

lbfgsfloatval_t penalizeData(const lbfgsfloatval_t *u, lbfgsfloatval_t *g){
	PolyhedralSurf::Vertex_const_iterator vb = RCSurface->P.vertices_begin();
	PolyhedralSurf::Vertex_const_iterator ve = RCSurface->P.vertices_end();
	PolyhedralSurf::Vertex_const_handle vh;

	lbfgsfloatval_t fx = 0.0;

	for (int i = 0; vb != ve; vb++,i++){
		vh = vb;

		if (occluded[i]){
			g[i*3] = 0;
			g[i*3+1] = 0;
			g[i*3+2] = 0;
			continue; // occluded region
		}

		Point_3 deformed = vh->point() + Vector_3(u[i*3],u[i*3+1],u[i*3+2]);
		Point_3 cor = CTSurface->vertexIndex[bestMatch[i]]->point();
		double weight = matchWeight[i];

		Vector_3 dis(deformed,cor);

		g[i*3] = -2 * dis.x() * weight;
		g[i*3+1] = -2 * dis.y() * weight;
		g[i*3+2] = -2 * dis.z() * weight;
		fx = fx + weight * dis.squared_length();
	}

	return fx;
}

lbfgsfloatval_t penalizeStretch(const lbfgsfloatval_t *u, lbfgsfloatval_t *g){
	double stretching = 0;

	if (regweight < 0){ // edge-based  stretching energy
		for (int i=0; i < RCSurface->edgeNum/2;i++){


			edge he = RCSurface->edgeList[i];

			double weight =  regweight * 2 * he.stretchWeight;
			int idx1 = he.index1;
			int idx2 = he.index2;

			Vector_3 deformVec(u[idx2*3] - u[idx1*3],u[idx2*3+1] - u[idx1*3+1],u[idx2*3+2] - u[idx1*3+2]);

			g[idx1*3] += -2 * deformVec.x() * weight;
			g[idx1*3+1] += -2 * deformVec.y() * weight;
			g[idx1*3+2] += -2 * deformVec.z() * weight;
			g[idx2*3] += 2 * deformVec.x() * weight;
			g[idx2*3+1] += 2 * deformVec.y() * weight;
			g[idx2*3+2] += 2 * deformVec.z() * weight;

			stretching += (deformVec.squared_length() * weight);
		}
	}else{ //triangle-based stretching energy
		facet* faceList = RCSurface->faceList;

		double Dtemp_v1_x,Dtemp_v1_y,Dtemp_v1_z;
		double Dtemp_v2_x,Dtemp_v2_y,Dtemp_v2_z;

		Point_3 newV1,newV2,newV3;
		Vector_3 temp_v1,temp_v2;
		double v11,v12,v21,v22,J11,J12,J21,J22,S11,S12,S21,S22,k;
		double Dv11,Dv12,Dv21,Dv22,DJ11,DJ12,DJ21,DJ22,DS11,DS12,DS21,DS22,dk;
		double weight,trace,trace_2,det,det_inv,Ddet;

		double* tempS  = new double[RCSurface->faceNum];
		memset(tempS,0,sizeof(double)*RCSurface->faceNum);

				cout<<"complete 1 single energy computing"<<endl;

		#pragma omp parallel for private(idx)
		for (int idx=0; idx < RCSurface->faceNum;idx++){
			facet f = faceList[idx];
			weight = f.area * regweight;

			newV1 = f.v1->point() + Vector_3(u[f.index[1]*3],u[f.index[1]*3+1],u[f.index[1]*3+2]);
			newV2 = f.v2->point() + Vector_3(u[f.index[2]*3],u[f.index[2]*3+1],u[f.index[2]*3+2]);
			newV3 = f.v3->point() + Vector_3(u[f.index[3]*3],u[f.index[3]*3+1],u[f.index[3]*3+2]);

			temp_v1 = Vector_3(newV1,newV2);
			temp_v2 = Vector_3(newV1,newV3);

			v11 = 0;
			v12 = sqrt(temp_v1.squared_length()+EPS);

			k = temp_v1*temp_v2;
			v22 = k / v12;

			double sqv21 = temp_v2.squared_length() - v22*v22;

			v21 = - sqrt(sqv21+EPS);

			J11 = 0*f.inverse[0][0] + v21*f.inverse[1][0];
			J12 = 0*f.inverse[0][1] + v21*f.inverse[1][1];
			J21 = v12*f.inverse[0][0] + v22*f.inverse[1][0];
			J22 = v12*f.inverse[0][1] + v22*f.inverse[1][1];

			S11 = J11*J11 + J21*J21; S12 = J11*J12+J21*J22; S21 = S12; S22 = J12*J12 + J22*J22;
			trace = (S11 + S22 - 2);
			
			//det = S11*S22 - S12*S21;
			//det_inv = 1 / det;
			trace_2 = 0.5*(S11-1)*(S11-1) + 4*S12*S21 + 0.5*(S22-1)*(S22-1);

			//stretching += weight * (STRETCH_MIU / 2 * trace + (STRETCH_LAMBDA - STRETCH_MIU*2)/8 * det + (STRETCH_LAMBDA + STRETCH_MIU*2)/8 * det_inv);
			//stretching += weight *YOUNG/2/(1-POISSON*POISSON)*((1-POISSON)*trace_2+POISSON*trace*trace);
			tempS[idx] = weight *YOUNG/2/(1-POISSON*POISSON)*((1-POISSON)*trace_2+POISSON*trace*trace);

			//cout<<"complete 1.1 single energy computing"<<endl;

			/* COMPUTE GRADIENT */
			for (int i = 1; i < 4; i++)
				for (int j = 1; j < 4; j++){
					Dtemp_v1_x = computeStretchGradient(1,1,i,j);Dtemp_v1_y = computeStretchGradient(1,2,i,j);Dtemp_v1_z = computeStretchGradient(1,3,i,j);
					Dtemp_v2_x = computeStretchGradient(2,1,i,j);Dtemp_v2_y = computeStretchGradient(2,2,i,j);Dtemp_v2_z = computeStretchGradient(2,3,i,j);

					dk = (Dtemp_v1_x*temp_v2.x()  +temp_v1.x()*Dtemp_v2_x) + (Dtemp_v1_y*temp_v2.y()  +temp_v1.y()*Dtemp_v2_y) + (Dtemp_v1_z*temp_v2.z()  +temp_v1.z()*Dtemp_v2_z);
					Dv12 = 0.5 * 1 / v12 * (2*temp_v1.x()*Dtemp_v1_x + 2*temp_v1.y()*Dtemp_v1_y + 2*temp_v1.z()*Dtemp_v1_z);
					Dv22 = (dk * v12 - Dv12 * k) / (v12*v12);
					Dv21 = 0.5 * 1 / v21 * (2*temp_v2.x()*Dtemp_v2_x + 2*temp_v2.y()*Dtemp_v2_y + 2*temp_v2.z()*Dtemp_v2_z - 2*v22*Dv22);

					DJ11 = Dv21*f.inverse[1][0];
					DJ12 = Dv21*f.inverse[1][1];
					DJ21 = Dv12*f.inverse[0][0] + Dv22*f.inverse[1][0];
					DJ22 = Dv12*f.inverse[0][1] + Dv22*f.inverse[1][1];

					DS11 = 2*J11*DJ11 + 2*J21*DJ21;
					DS22 = 2*J12*DJ12 + 2*J22*DJ22;
					DS12 = (DJ11*J12+J11*DJ12) + (DJ21*J22+J21*DJ22);
					DS21 = DS12;

					//Ddet = (DS11*S22+S11*DS22) - (DS12*S21+S12*DS21);

					//g[f.index[i]*3+j-1] += weight * STRETCH_MIU / 2 * (DS11 + DS22);
					//g[f.index[i]*3+j-1] += weight * (STRETCH_LAMBDA - STRETCH_MIU*2)/8 * Ddet;
					//g[f.index[i]*3+j-1] += weight * (STRETCH_LAMBDA + STRETCH_MIU*2)/8 * (-1/(det*det)*Ddet);
					
// 					g[f.index[i]*3+j-1] += weight *YOUNG/2/(1-POISSON*POISSON)*
// 						((1-POISSON)*(1*(S11-1)*DS11+4*S12*DS21+4*DS12*S21+1*(S22-1)*DS22)+
// 						 POISSON*2*trace*(DS11+DS22));
					tempG[idx*9+(i-1)*3+j-1] = weight *YOUNG/2/(1-POISSON*POISSON)*
						 					   ((1-POISSON)*(1*(S11-1)*DS11+4*S12*DS21+4*DS12*S21+1*(S22-1)*DS22)+
						    				   POISSON*2*trace*(DS11+DS22));
				}


		}

		cout<<"complete 2 single energy computing"<<endl;

		for (int idx=0; idx < RCSurface->faceNum;idx++){
			stretching += tempS[idx];
			for (int i = 1; i < 4; i++)
				for (int j = 1; j < 4; j++)
					g[faceList[idx].index[i]*3+j-1] += tempG[idx*9+(i-1)*3+j-1];
		}	
		delete[] tempS;
	}

	facetStretching = stretching;
	return stretching;
}

lbfgsfloatval_t penalizeBend(const lbfgsfloatval_t *u, lbfgsfloatval_t *g){
	double bending = 0;
	facet* faceList = RCSurface->faceList;

	Point_3 newV1,newV2,newV3,newNB1,newNB2,newNB3;
	double dv1x, dv1y, dv1z, dv2x, dv2y, dv2z, dv3x, dv3y, dv3z,
		   dnb1x,dnb1y,dnb1z,dnb2x,dnb2y,dnb2z,dnb3x,dnb3y,dnb3z;
	double Dtheta1,Dtheta2,Dtheta3,DSOD_00,DSOD_01,DSOD_10,DSOD_11;
	double weight,Dbend;

	for (int idx=0; idx < RCSurface->faceNum;idx++){
		facet f = faceList[idx];
		if (f.isBorder) continue;
		weight = bendweight;

		newV1 = f.v1->point() + Vector_3(u[f.index[1]*3],u[f.index[1]*3+1],u[f.index[1]*3+2]);
		newV2 = f.v2->point() + Vector_3(u[f.index[2]*3],u[f.index[2]*3+1],u[f.index[2]*3+2]);
		newV3 = f.v3->point() + Vector_3(u[f.index[3]*3],u[f.index[3]*3+1],u[f.index[3]*3+2]);
		newNB1 = f.nb1->point() + Vector_3(u[f.nbIdx1*3],u[f.nbIdx1*3+1],u[f.nbIdx1*3+2]);
		newNB2 = f.nb2->point() + Vector_3(u[f.nbIdx2*3],u[f.nbIdx2*3+1],u[f.nbIdx2*3+2]);
		newNB3 = f.nb3->point() + Vector_3(u[f.nbIdx3*3],u[f.nbIdx3*3+1],u[f.nbIdx3*3+2]);

		/* */
		double* det1 = computeDeterminant(newV1,newV2,newV3,newNB1);
		double* det2 = computeDeterminant(newV1,newV2,newV3,newNB2);
		double* det3 = computeDeterminant(newV1,newV2,newV3,newNB3);
		double theta1 = - det1[0] / ((f.area*f.sideArea1)/f.l1);
		double theta2 = - det2[0] / ((f.area*f.sideArea2)/f.l2);
		double theta3 = - det3[0] / ((f.area*f.sideArea3)/f.l3);
		/* */
		/* angle based 
		double theta1 = computeAngle(newV1,newV2,newV3,newNB1);
		double theta2 = computeAngle(newV2,newV3,newV1,newNB2);
		double theta3 = computeAngle(newV3,newV1,newV2,newNB3);
		*/

		double SOD_00 = (theta1-f.theta1) * f.SO1[0][0] + (theta2-f.theta2) * f.SO2[0][0] + (theta3-f.theta3) * f.SO3[0][0];
		double SOD_01 = (theta1-f.theta1) * f.SO1[0][1] + (theta2-f.theta2) * f.SO2[0][1] + (theta3-f.theta3) * f.SO3[0][1];
		double SOD_10 = (theta1-f.theta1) * f.SO1[1][0] + (theta2-f.theta2) * f.SO2[1][0] + (theta3-f.theta3) * f.SO3[1][0];
		double SOD_11 = (theta1-f.theta1) * f.SO1[1][1] + (theta2-f.theta2) * f.SO2[1][1] + (theta3-f.theta3) * f.SO3[1][1];

		bending += weight * (SOD_00*SOD_00 + SOD_01*SOD_01 + SOD_10*SOD_10 + SOD_11*SOD_11); 

// 		if ((f.index[1] == 2640)&& (f.index[2] == 2638) && (f.index[3] == 2570)){
// 			facetBending = weight * (SOD_00*SOD_00 + SOD_01*SOD_01 + SOD_10*SOD_10 + SOD_11*SOD_11);
// 			thetaDif1 = theta1 - f.theta1;
// 			thetaDif2 = theta2 - f.theta2;
// 			thetaDif3 = theta3 - f.theta3;
// 		}

		//COMPUTE GRADIENT

		for (int i = 1; i < 4; i++)
			for (int j = 1; j < 4; j++){
				/* */
				Dtheta1 = - computeDetGradient(i,j,newV1,newV2,newV3,newNB1,det1) / ((f.area*f.sideArea1)/f.l1);
				Dtheta2 = - computeDetGradient(i,j,newV1,newV2,newV3,newNB2,det2) / ((f.area*f.sideArea2)/f.l2);
				Dtheta3 = - computeDetGradient(i,j,newV1,newV2,newV3,newNB3,det3) / ((f.area*f.sideArea3)/f.l3);
				/* */
				/* angle based
				int i1 = i;
				int i2 = ((i+2)>3)? i-1 : i+2;
				int i3 = ((i+1)>3)? i-2 : i+1;

				cout<<i1<<i2<<i3;
				cin>>i1;

				Dtheta1 = computeAngleGradient(i1,j,newV1,newV2,newV3,newNB1);
				Dtheta2 = computeAngleGradient(i2,j,newV2,newV3,newV1,newNB2);
				Dtheta3 = computeAngleGradient(i3,j,newV3,newV1,newV2,newNB3);
				*/

				DSOD_00 = Dtheta1 * f.SO1[0][0] + Dtheta2 * f.SO2[0][0] + Dtheta3 * f.SO3[0][0];
				DSOD_01 = Dtheta1 * f.SO1[0][1] + Dtheta2 * f.SO2[0][1] + Dtheta3 * f.SO3[0][1];
				DSOD_10 = Dtheta1 * f.SO1[1][0] + Dtheta2 * f.SO2[1][0] + Dtheta3 * f.SO3[1][0];
				DSOD_11 = Dtheta1 * f.SO1[1][1] + Dtheta2 * f.SO2[1][1] + Dtheta3 * f.SO3[1][1];

				Dbend = weight * (2*SOD_00*DSOD_00 + 2*SOD_01*DSOD_01 + 2*SOD_10*DSOD_10 + 2*SOD_11*DSOD_11);

				g[f.index[i]*3+j-1] += Dbend;
			}

		for (int j = 1; j < 4; j++){
			Dtheta1 = - computeDetGradient(4,j,newV1,newV2,newV3,newNB1,det1) / ((f.area*f.sideArea1)/f.l1);
			Dtheta2 = 0;
			Dtheta3 = 0;

			DSOD_00 = Dtheta1 * f.SO1[0][0] + Dtheta2 * f.SO2[0][0] + Dtheta3 * f.SO3[0][0];
			DSOD_01 = Dtheta1 * f.SO1[0][1] + Dtheta2 * f.SO2[0][1] + Dtheta3 * f.SO3[0][1];
			DSOD_10 = Dtheta1 * f.SO1[1][0] + Dtheta2 * f.SO2[1][0] + Dtheta3 * f.SO3[1][0];
			DSOD_11 = Dtheta1 * f.SO1[1][1] + Dtheta2 * f.SO2[1][1] + Dtheta3 * f.SO3[1][1];

			Dbend = weight * (2*SOD_00*DSOD_00 + 2*SOD_01*DSOD_01 + 2*SOD_10*DSOD_10 + 2*SOD_11*DSOD_11);

			g[f.nbIdx1*3+j-1] += Dbend;
		}
		for (int j = 1; j < 4; j++){
			Dtheta1 = 0;
			Dtheta2 = - computeDetGradient(4,j,newV1,newV2,newV3,newNB2,det2) / ((f.area*f.sideArea2)/f.l2);
			Dtheta3 = 0;

			DSOD_00 = Dtheta1 * f.SO1[0][0] + Dtheta2 * f.SO2[0][0] + Dtheta3 * f.SO3[0][0];
			DSOD_01 = Dtheta1 * f.SO1[0][1] + Dtheta2 * f.SO2[0][1] + Dtheta3 * f.SO3[0][1];
			DSOD_10 = Dtheta1 * f.SO1[1][0] + Dtheta2 * f.SO2[1][0] + Dtheta3 * f.SO3[1][0];
			DSOD_11 = Dtheta1 * f.SO1[1][1] + Dtheta2 * f.SO2[1][1] + Dtheta3 * f.SO3[1][1];

			Dbend = weight * (2*SOD_00*DSOD_00 + 2*SOD_01*DSOD_01 + 2*SOD_10*DSOD_10 + 2*SOD_11*DSOD_11);

			g[f.nbIdx2*3+j-1] += Dbend;
		}
		for (int j = 1; j < 4; j++){
			Dtheta1 = 0;
			Dtheta2 = 0;
			Dtheta3 = - computeDetGradient(4,j,newV1,newV2,newV3,newNB3,det3) / ((f.area*f.sideArea3)/f.l3);
			
			DSOD_00 = Dtheta1 * f.SO1[0][0] + Dtheta2 * f.SO2[0][0] + Dtheta3 * f.SO3[0][0];
			DSOD_01 = Dtheta1 * f.SO1[0][1] + Dtheta2 * f.SO2[0][1] + Dtheta3 * f.SO3[0][1];
			DSOD_10 = Dtheta1 * f.SO1[1][0] + Dtheta2 * f.SO2[1][0] + Dtheta3 * f.SO3[1][0];
			DSOD_11 = Dtheta1 * f.SO1[1][1] + Dtheta2 * f.SO2[1][1] + Dtheta3 * f.SO3[1][1];

			Dbend = weight * (2*SOD_00*DSOD_00 + 2*SOD_01*DSOD_01 + 2*SOD_10*DSOD_10 + 2*SOD_11*DSOD_11);

			g[f.nbIdx3*3+j-1] += Dbend;
		}

		delete[] det1;
		delete[] det2;
		delete[] det3;
	}

	facetBending = bending;
	return bending;
}

lbfgsfloatval_t penalizeBendQuadratic(const lbfgsfloatval_t *tu, lbfgsfloatval_t *g){
	PolyhedralSurf::Vertex_const_iterator vb = RCSurface->P.vertices_begin();
	PolyhedralSurf::Vertex_const_iterator ve = RCSurface->P.vertices_end();
	PolyhedralSurf::Vertex_const_handle vh,nb;
	vertex* vl = RCSurface->vertexList;

	double fv = 0;
	int idx1,idx2;
	double laplacian,l1,l2,l3;
	for (int i = 0; vb != ve; vb++,i++){
		vh = vb;

		l1 = 0;l2 = 0;l3 = 0;

		for (int j = 0; j < vl[i].neighbourNum; j++){
			idx1 = vl[i].idx;
			idx2 = vl[i].nbIdx[j];

			l1 += vl[i].weight[j]*(tu[idx2*3] -   tu[idx1*3]);
			l2 += vl[i].weight[j]*(tu[idx2*3+1] - tu[idx1*3+1]);
			l3 += vl[i].weight[j]*(tu[idx2*3+2] - tu[idx1*3+2]);
		}

		double weight = bendweight / (vl[i].area*vl[i].area);
		fv += weight * (l1*l1 + l2*l2 + l3*l3);

		for (int j = 0; j < vl[i].neighbourNum; j++){
			idx1 = vl[i].idx;
			idx2 = vl[i].nbIdx[j];

			g[idx2*3]   += weight * 2 * l1 * vl[i].weight[j];
			g[idx2*3+1] += weight * 2 * l2 * vl[i].weight[j];
			g[idx2*3+2] += weight * 2 * l3 * vl[i].weight[j];

			g[idx1*3]   -= weight * 2 * l1 * vl[i].weight[j];
			g[idx1*3+1] -= weight * 2 * l2 * vl[i].weight[j];
			g[idx1*3+2] -= weight * 2 * l3 * vl[i].weight[j];
		}
	}

	facetBending = fv;
	return fv;
}

lbfgsfloatval_t penalizeLink(const lbfgsfloatval_t *u, lbfgsfloatval_t *g){
	float linkEnergy = 0;
	
	for (int i = 0; i < affinityM; i++){
		if (!attractor[i]) continue;

		int idx1 = i;
		int idx2 = linkTarget[i];

		Point_3 P1 = RCSurface->vertexIndex[idx1]->point() + Vector_3(u[idx1*3],u[idx1*3+1],u[idx1*3+2]);
		Point_3 P2 = RCSurface->vertexIndex[idx2]->point() + Vector_3(u[idx2*3],u[idx2*3+1],u[idx2*3+2]);
		Vector_3 dis(P1,P2);

		float linkLength = sqrt(dis.squared_length());
		if (linkLength < 0.00001) return 0;
		float a = 4;
		float linkWeight = 5000;

		if (linkLength < a){
			linkEnergy += linkWeight * (a*a/6) * (1-pow((1-linkLength*linkLength/a/a),3));

			float dB = linkWeight * linkLength * pow((1-linkLength*linkLength/a/a),2);
			g[idx1*3] += dB * (- 0.5* 1 / linkLength * dis.x());
			g[idx1*3+1] += dB * (- 0.5* 1 / linkLength * dis.y());
			g[idx1*3+2] += dB * (- 0.5* 1 / linkLength * dis.z());

			g[idx2*3] += dB * (0.5* 1 / linkLength * dis.x());
			g[idx2*3+1] += dB * (0.5* 1 / linkLength * dis.y());
			g[idx2*3+2] += dB * (0.5* 1 / linkLength * dis.z());
		}
		else{ 
			linkEnergy += linkWeight * (a*a/6);
		}
	}
	
	distantLink = linkEnergy;

	return linkEnergy;
}

lbfgsfloatval_t penalizaeLandmark(const lbfgsfloatval_t *u, lbfgsfloatval_t *g){
	PolyhedralSurf::Vertex_const_handle vh;
	double weight = LANDMARKWEIGHT;

	lbfgsfloatval_t fx = 0.0;
	for (int i = 0; i<RCSurface->landmarkNum; i++){
		int idx = RCSurface->landmark[i];
		vh = RCSurface->vertexIndex[idx];
		Point_3 deformed = vh->point() + Vector_3(u[idx*3],u[idx*3+1],u[idx*3+2]);
		Vector_3 dis(deformed,CTSurface->vertexIndex[CTSurface->landmark[i]]->point());

		fx = fx + weight*dis.squared_length();
		g[idx*3] += -2 * dis.x() * weight;
		g[idx*3+1] += -2 * dis.y() * weight;
		g[idx*3+2] += -2 * dis.z() * weight;

	}

	return fx;
}

double computeStretchGradient(int n,int m,int x,int y){
	if (x == 1){
		if (m != y) return 0; else return -1;
	}else if (x == 2){
		if (n == 2) return 0;
		if (m != y) return 0; else return 1;
	}if (x == 3){
		if (n == 1) return 0;
		if (m != y) return 0; else return 1;
	}
}

double computeAngleGradient(int i,int j,Point_3 V1,Point_3 V2,Point_3 V3,Point_3 V4){
	Vector_3 v1(V1,V4);
	Vector_3 v2(V1,V2);
	Vector_3 v3(V1,V3);

	Vector_3 n2 = cross_product(v2,v3);
	double sm2 = n2.squared_length();
	double sq2 = sqrt(sm2);
	double p1 = v1*n2;
	double dis1 = p1/sq2;

	double v1v2 = v1*v2;
	double sqv2 = sqrt(v2.squared_length());
	double p2 = v1v2/sqv2;
	double dis2 = sqrt(v1.squared_length()-p2*p2);

	/* gradient */

	double dn2x,dn2y,dn2z,dsm2,dsq2,dp1,ddis1,dtheta,dsv2,dsqv2,dv1v2,dp2,ddis2;
	double dv1x,dv1y,dv1z,dv2x,dv2y,dv2z,dv3x,dv3y,dv3z;

	dv1x = computeGradient(1,1,i,j);dv1y = computeGradient(1,2,i,j); dv1z = computeGradient(1,3,i,j);
	dv2x = computeGradient(2,1,i,j);dv2y = computeGradient(2,2,i,j); dv2z = computeGradient(2,3,i,j);
	dv3x = computeGradient(3,1,i,j);dv3y = computeGradient(3,2,i,j); dv3z = computeGradient(3,3,i,j);

	dn2x = (dv2y*v3.z()+v2.y()*dv3z)-(dv2z*v3.y()+v2.z()*dv3y);
	dn2y = (dv2z*v3.x()+v2.z()*dv3x)-(dv2x*v3.z()+v2.x()*dv3z);
	dn2z = (dv2x*v3.y()+v2.x()*dv3y)-(dv2y*v3.x()+v2.y()*dv3x);

	dsm2 = 2*n2.x()*dn2x + 2*n2.y()*dn2y + 2*n2.z()*dn2z;
	dsq2 = 0.5 * 1 / sq2 * dsm2;

	dp1 = (dv1x*n2.x()+v1.x()*dn2x) + (dv1y*n2.y()+v1.y()*dn2y) + (dv1z*n2.z()+v1.z()*dn2z);
	ddis1 = (dp1*sq2-p1*dsq2)/(sq2*sq2);
	
	//
	dsv2 = 2*v2.x()*dv2x + 2*v2.y()*dv2y + 2*v2.z()*dv2z;
	dsqv2 = 0.5 * 1 / sqv2 * dsv2;

	dv1v2 = (dv1x*v2.x()+v1.x()*dv2x) + (dv1y*v2.y()+v1.y()*dv2y) + (dv1z*v2.z()+v1.z()*dv2z);
	dp2 = (dv1v2*sqv2 - v1v2*dsqv2)/v2.squared_length();

	ddis2 = 0.5 * 1 / dis2 *(2*v1.x()*dv1x + 2*v1.y()*dv1y + 2*v1.z()*dv1z-2*p2*dp2);

	return (ddis1*dis2-dis1*ddis2)/(dis2*dis2);
}

double computeGradient(int n,int nd,int e,int ed){
	if (nd != ed) return 0;

	if (n == 1){
		if (e == 1) return -1;
		else if (e == 4) return 1;
		else return 0;
	}
	else if (n == 2){
		if (e == 1) return -1;
		else if (e == 2) return 1;
		else return 0;
	}
	else if (n == 3){
		if (e == 1) return -1;
		else if (e == 3) return 1;
		else return 0;
	}

	return 0;
}

double computeDetGradient(int i,int j,Point_3 v1,Point_3 v2,Point_3 v3,Point_3 v4,double* det){
	double D1,D2,D3,D4,Ddet;

	if (i == 1){
		switch (j)
		{
		case 1:
			Ddet = -det[2];
			break;
		case 2:
			Ddet = det[3];
			break;
		case 3:
			Ddet = -det[4];
		}
	}

	if (i == 2){
		switch (j)
		{
		case 1:
			D1 = v3.y() * v4.z() -v3.z() * v4.y(); 
			D2 = 0;
			D3 = v3.z() - v4.z();
			D4 = v3.y() - v4.y();
			Ddet = D1 - v1.x()*D2 + v1.y()*D3 - v1.z()*D4;
			break;
		case 2:
			D1 = v3.z() * v4.x() - v3.x() * v4.z();
			D2 = v3.z() - v4.z();
			D3 = 0;
			D4 = v4.x() - v3.x();
			Ddet = D1 - v1.x()*D2 + v1.y()*D3 - v1.z()*D4;
			break;

		case 3:
			D1 = v3.x() * v4.y() - v3.y() * v4.x();
			D2 = v4.y() - v3.y();
			D3 = v4.x() - v3.x();
			D4 = 0;
			Ddet = D1 - v1.x()*D2 + v1.y()*D3 - v1.z()*D4;
		}
	}

	if (i == 3){
		switch (j)
		{
		case 1:
			D1 = v2.z() * v4.y() - v2.y() * v4.z();
			D2 = 0;
			D3 = v4.z() - v2.z();
			D4 = v4.y() - v2.y();
			Ddet = D1 - v1.x()*D2 + v1.y()*D3 - v1.z()*D4;
			break;

		case 2:
			D1 = v2.x() * v4.z() - v2.z() * v4.x();
			D2 = v4.z() - v2.z();
			D3 = 0;
			D4 = v2.x() - v4.x();
			Ddet = D1 - v1.x()*D2 + v1.y()*D3 - v1.z()*D4;
			break;

		case 3:
			D1 = v2.y() * v4.x() - v2.x() * v4.y();
			D2 = v2.y() - v4.y();
			D3 = v2.x() - v4.x();
			D4 = 0;
			Ddet = D1 - v1.x()*D2 + v1.y()*D3 - v1.z()*D4;
		}
	}

	if (i == 4){
		switch (j)
		{
		case 1:
			D1 = v2.y() * v3.z() - v2.z() * v3.y();
			D2 = 0;
			D3 = v2.z() - v3.z();
			D4 = v2.y() - v3.y();
			Ddet = D1 - v1.x()*D2 + v1.y()*D3 - v1.z()*D4;
			break;

		case 2:
			D1 = v2.z() * v3.x() - v2.x() * v3.z();
			D2 = v2.z() - v3.z();
			D3 = 0;
			D4 = v3.x() - v2.x();
			Ddet = D1 - v1.x()*D2 + v1.y()*D3 - v1.z()*D4;
			break;

		case 3:
			D1 = v2.x() * v3.y() - v2.y() * v3.x();
			D2 = v3.y() - v2.y();
			D3 = v3.x() - v2.x();
			D4 = 0;
			Ddet = D1 - v1.x()*D2 + v1.y()*D3 - v1.z()*D4;
		}
	}

	return Ddet;
}
