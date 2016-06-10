#include "stdafx.h"
#include "GlutManager.h"
#include <omp.h>
#include <iostream>
#include <fstream>
using namespace std;

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
	if ((k % 100 == 0) || (k <=1)){
		int idx = *((int *)instance);
		memcpy(Surface[idx]->u,tu,sizeof(lbfgsfloatval_t)*Surface[idx]->vertexNum*3);

		cout<<"Surface: "<<idx<<" Iteration: "<<k<<" -------------------------------"<<endl;
		cout<<"fx = "<<fx<<endl;
		cout<<"xnorm = "<<xnorm<<", gnorm = "<<gnorm<<", step = "<<step<<endl;
		cout<<endl;
		cout<<"Bending: "<<facetBending<<" Stretching: "<<facetStretching<<" Link: "<<distantLink<<endl;
		cout<<thetaDif1<<' '<<thetaDif2<<' '<<thetaDif3<<endl;

	}
	return 0;
}

void startOptimization(){
	/* initialize optimization */

	for (int iter = 0; iter < 20; iter++){


		int i;
		#pragma omp parallel for private(i)
		for (i = 0; i < surfaceNum; i++){
			lbfgsfloatval_t *tempU = new lbfgsfloatval_t[Surface[i]->vertexNum*3];
			memset(tempU,0,sizeof(lbfgsfloatval_t) * Surface[i]->vertexNum * 3);

			lbfgs_parameter_t param;
			lbfgsfloatval_t fx = 0;
			lbfgs_parameter_init(&param);
			int opIterNum = 20000;
			param.max_iterations = opIterNum;
			param.epsilon = 1e-5f;

			int idx = i;
			int ret = lbfgs(Surface[i]->vertexNum * 3, tempU, &fx, evaluate, progress, (void*)(&idx), &param);

			char prefix[5];
			//itoa(iter,prefix,10);
            sprintf(prefix, "%d", iter);
			Surface[i]->if_name = fileName[i]+prefix+".off";
			writeOFF(Surface[i],fileName[i]+prefix+".off",i);
			cout<<"Finishing surface "<<i<<endl;

			if ((MESHLABOPTION > 0) && (iter % MESHLABOPTION == 0)) {
				cout<<"PostProcessing"<<endl;
				char cmdLine[1000];
				sprintf(cmdLine,"meshlabserver -i %s -o %s -s meshlabscript_demons.mlx\n",Surface[i]->if_name.c_str(),Surface[i]->if_name.c_str());
				cout<<cmdLine<<endl;
				system(cmdLine);
			}
		}

		REGWEIGHT *= 0.9;
		BENDWEIGHT *= 0.8;
		LANDMARKWEIGHT *= 0.9;
		EUCLWEIGHT *= 0.9;


		//if ((iter % 3 == 0)&&(iter != 0)) GEOMETRIC_RANGE *= 2;

		cout<<"Beginning next iteration"<<endl;

		for (int i = 0; i < surfaceNum; i++){
			Surface[i]->ComputeMeshProperty();
			Surface[i]->findSignatureAll();
			Surface[i]->constructEdgeList();
			Surface[i]->constructVertexList();
			Surface[i]->initialDeformation();
		}

		for (int i = 0; i < surfaceNum; i++){
			cout<<"Computing attraction force: "<<fileName[i]<<endl;

			int startSurface = max_zenyo(i - GEOMETRIC_RANGE, 0);
			int endSurface = min_zenyo(i + GEOMETRIC_RANGE, surfaceNum - 1);

			for (int j = startSurface; j <= endSurface; j++)
				if (i != j) {
					Surface[i]->findCorrespondenceBothWay(Surface[j],EUCLWEIGHT);
					Surface[i]->findMatch2(Surface[j],j);
				}
				Surface[i]->summerizeForce();
		}
	}

}

static lbfgsfloatval_t evaluate(
	void *instance,
	const lbfgsfloatval_t *u,
	lbfgsfloatval_t *g,
	const int n,
	const lbfgsfloatval_t step
	)
{
	int idx = *((int *)instance);//cout<<"begin optimization"<<endl;
	//cout<<idx<<endl;

	//double data = 0;
	double data = penalizeData(u,g,idx);

	/* stretching */
	double stretching = penalizeStretchQuadratic(u,g,idx);
	//double stretching = penalizeStretch(u,g,idx);
	//double stretching = 0;
	
	/* bending */
	double bending = penalizeBendQuadratic(u,g,idx);
	//double bending = penalizeBendAngleLP(u,g,idx);
	//double bending = 0;

	double landmark  = 0;
	if (Surface[idx]->landmarkNum != NULL) 
		landmark = penalizeLandmark(u,g,idx);

	return data + stretching +  bending + landmark;
}

lbfgsfloatval_t penalizeData(const lbfgsfloatval_t *u, lbfgsfloatval_t *g, int idx){
	PolyhedralSurf::Vertex_const_iterator vb = Surface[idx]->P.vertices_begin();
	PolyhedralSurf::Vertex_const_iterator ve = Surface[idx]->P.vertices_end();
	PolyhedralSurf::Vertex_const_handle vh;

	lbfgsfloatval_t fx = 0.0;

	int i = 0;
	//#pragma omp parallel for private(i)
	for (i = 0; vb != ve; vb++,i++){
		vh = vb;

		Point_3 deformed = vh->point() + Vector_3(u[i*3],u[i*3+1],u[i*3+2]);
		Point_3 cor = Point_3(Surface[idx]->targetPos[i].x(),Surface[idx]->targetPos[i].y(),Surface[idx]->targetPos[i].z());
		double weight = Surface[idx]->overallWeight[i];

		Vector_3 dis(deformed,cor);

		
		g[i*3] = -2 * dis.x() * weight;
		g[i*3+1] = -2 * dis.y() * weight;
		g[i*3+2] = -2 * dis.z() * weight;
		fx = fx + weight * dis.squared_length();
	}

	return fx;
}

lbfgsfloatval_t penalizeStretchQuadratic(const lbfgsfloatval_t *u, lbfgsfloatval_t *g, int surfaceIdx){
	    // edge-based  stretching energy
	double stretching = 0;
		for (int i=0; i < Surface[surfaceIdx]->edgeNum/2;i++){
			edge he = Surface[surfaceIdx]->edgeList[i];

			double weight =  REGWEIGHT * 2 * he.stretchWeight;
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

	facetStretching = stretching;
	return stretching;
}

lbfgsfloatval_t penalizeStretch(const lbfgsfloatval_t *u, lbfgsfloatval_t *g, int surfaceIdx){
	double stretching = 0;

	{ //triangle-based stretching energy
		facet* faceList = Surface[surfaceIdx]->faceList;

		int idx;
		for ( idx=0; idx < Surface[surfaceIdx]->faceNum;idx++){
			facet f = faceList[idx];
			double weight = f.area * REGWEIGHT;

			Point_3 newV1 = f.v1->point() + Vector_3(u[f.index[1]*3],u[f.index[1]*3+1],u[f.index[1]*3+2]);
			Point_3 newV2 = f.v2->point() + Vector_3(u[f.index[2]*3],u[f.index[2]*3+1],u[f.index[2]*3+2]);
			Point_3 newV3 = f.v3->point() + Vector_3(u[f.index[3]*3],u[f.index[3]*3+1],u[f.index[3]*3+2]);

			Vector_3 temp_v1 = Vector_3(newV1,newV2);
			Vector_3 temp_v2 = Vector_3(newV1,newV3);

			double v11 = 0;
			double v12 = sqrt(temp_v1.squared_length()+EPS);

			double k = temp_v1*temp_v2;
			double v22 = k / v12;

			double sqv21 = temp_v2.squared_length() - v22*v22;

			double v21 = - sqrt(sqv21+EPS);

			double J11 = 0*f.inverse[0][0] + v21*f.inverse[1][0];
			double J12 = 0*f.inverse[0][1] + v21*f.inverse[1][1];
			double J21 = v12*f.inverse[0][0] + v22*f.inverse[1][0];
			double J22 = v12*f.inverse[0][1] + v22*f.inverse[1][1];

			double S11 = J11*J11 + J21*J21; double S12 = J11*J12+J21*J22; double S21 = S12; double S22 = J12*J12 + J22*J22;
			double trace = (S11 + S22 - 2);
			
			double trace_2 = 0.5*(S11-1)*(S11-1) + 4*S12*S21 + 0.5*(S22-1)*(S22-1);

			stretching += weight *YOUNG/2/(1-POISSON*POISSON)*((1-POISSON)*trace_2+POISSON*trace*trace);

			/* COMPUTE GRADIENT */
			for (int i = 1; i < 4; i++)
				for (int j = 1; j < 4; j++){
					double Dtemp_v1_x = computeStretchGradient(1,1,i,j);double Dtemp_v1_y = computeStretchGradient(1,2,i,j);double Dtemp_v1_z = computeStretchGradient(1,3,i,j);
					double Dtemp_v2_x = computeStretchGradient(2,1,i,j);double Dtemp_v2_y = computeStretchGradient(2,2,i,j);double Dtemp_v2_z = computeStretchGradient(2,3,i,j);

					double dk = (Dtemp_v1_x*temp_v2.x()  +temp_v1.x()*Dtemp_v2_x) + (Dtemp_v1_y*temp_v2.y()  +temp_v1.y()*Dtemp_v2_y) + (Dtemp_v1_z*temp_v2.z()  +temp_v1.z()*Dtemp_v2_z);
					double Dv12 = 0.5 * 1 / v12 * (2*temp_v1.x()*Dtemp_v1_x + 2*temp_v1.y()*Dtemp_v1_y + 2*temp_v1.z()*Dtemp_v1_z);
					double Dv22 = (dk * v12 - Dv12 * k) / (v12*v12);
					double Dv21 = 0.5 * 1 / v21 * (2*temp_v2.x()*Dtemp_v2_x + 2*temp_v2.y()*Dtemp_v2_y + 2*temp_v2.z()*Dtemp_v2_z - 2*v22*Dv22);

					double DJ11 = Dv21*f.inverse[1][0];
					double DJ12 = Dv21*f.inverse[1][1];
					double DJ21 = Dv12*f.inverse[0][0] + Dv22*f.inverse[1][0];
					double DJ22 = Dv12*f.inverse[0][1] + Dv22*f.inverse[1][1];

					double DS11 = 2*J11*DJ11 + 2*J21*DJ21;
					double DS22 = 2*J12*DJ12 + 2*J22*DJ22;
					double DS12 = (DJ11*J12+J11*DJ12) + (DJ21*J22+J21*DJ22);
					double DS21 = DS12;

					g[f.index[i]*3+j-1] += weight *YOUNG/2/(1-POISSON*POISSON)*
						((1-POISSON)*(1*(S11-1)*DS11+4*S12*DS21+4*DS12*S21+1*(S22-1)*DS22)+
						 POISSON*2*trace*(DS11+DS22));

				}


		}
	}

	facetStretching = stretching;
	return stretching;
}

lbfgsfloatval_t penalizeBendQuadratic(const lbfgsfloatval_t *tu, lbfgsfloatval_t *g, int idx){
	PolyhedralSurf::Vertex_const_iterator vb = Surface[idx]->P.vertices_begin();
	PolyhedralSurf::Vertex_const_iterator ve = Surface[idx]->P.vertices_end();
	PolyhedralSurf::Vertex_const_handle vh,nb;
	vertex* vl = Surface[idx]->vertexList;

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

		double weight = BENDWEIGHT / (vl[i].area*vl[i].area);
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

lbfgsfloatval_t penalizeBendAngleLP(const lbfgsfloatval_t *u, lbfgsfloatval_t *g, int surfaceIdx){
	double bending = 0;
	facet* faceList = Surface[surfaceIdx]->faceList;

	Point_3 newV1,newV2,newV3,newNB1,newNB2,newNB3;
	double theta1,theta2,theta3;
	threeTuple db1,db2,db3,ds1,ds2,ds3,dnm,dn1,dn2,dn3;
	double dx1,dx2,dx3,Dtheta1,Dtheta2,Dtheta3,DSOD_00,DSOD_01,DSOD_10,DSOD_11;
	double weight,Dbend,trace,trace_2;

	for (int idx=0; idx < Surface[surfaceIdx]->faceNum;idx++){
		facet f = faceList[idx];
		if (f.isBorder) continue;
		weight = f.area * BENDWEIGHT;

		newV1 = f.v1->point() + Vector_3(u[f.index[1]*3],u[f.index[1]*3+1],u[f.index[1]*3+2]);
		newV2 = f.v2->point() + Vector_3(u[f.index[2]*3],u[f.index[2]*3+1],u[f.index[2]*3+2]);
		newV3 = f.v3->point() + Vector_3(u[f.index[3]*3],u[f.index[3]*3+1],u[f.index[3]*3+2]);
		newNB1 = f.nb1->point() + Vector_3(u[f.nbIdx1*3],u[f.nbIdx1*3+1],u[f.nbIdx1*3+2]);
		newNB2 = f.nb2->point() + Vector_3(u[f.nbIdx2*3],u[f.nbIdx2*3+1],u[f.nbIdx2*3+2]);
		newNB3 = f.nb3->point() + Vector_3(u[f.nbIdx3*3],u[f.nbIdx3*3+1],u[f.nbIdx3*3+2]);

		/* -------------- compute angle ------------------------- */
		/* ------------------------------------------------------ */
		Vector_3 b1(newV1,newV2);
		Vector_3 b2(newV2,newV3);
		Vector_3 b3(newV3,newV1);
		Vector_3 s1(newV1,newNB1);
		Vector_3 s2(newV2,newNB2);
		Vector_3 s3(newV3,newNB3);

		Vector_3 nm = cross_product(b1,b2);
		Vector_3 n1 = cross_product(s1,b1);
		Vector_3 n2 = cross_product(s2,b2);
		Vector_3 n3 = cross_product(s3,b3);

		double sign1 = s1*nm;
		double sign2 = s2*nm;
		double sign3 = s3*nm;

		double x1 = n1*nm / (sqrt(n1.squared_length()+EPS)*sqrt(nm.squared_length()+EPS));
		double x2 = n2*nm / (sqrt(n2.squared_length()+EPS)*sqrt(nm.squared_length()+EPS));
		double x3 = n3*nm / (sqrt(n3.squared_length()+EPS)*sqrt(nm.squared_length()+EPS));

		if (sign1 > 0) x1 = (2 - x1);
		if (sign2 > 0) x2 = (2 - x2);
		if (sign3 > 0) x3 = (2 - x3);

		if (x1 <= 1) theta1 = (-0.69813170079773212 * x1 * x1 - 0.87266462599716477) * x1 + 1.5707963267948966;
		else	     theta1 = (-0.69813170079773212 * (x1-2) * (x1-2) - 0.87266462599716477) * (x1-2) - 1.5707963267948966;
		if (x2 <= 1) theta2 = (-0.69813170079773212 * x2 * x2 - 0.87266462599716477) * x2 + 1.5707963267948966;
		else	     theta2 = (-0.69813170079773212 * (x2-2) * (x2-2) - 0.87266462599716477) * (x2-2) - 1.5707963267948966;
		if (x3 <= 1) theta3 = (-0.69813170079773212 * x3 * x3 - 0.87266462599716477) * x3 + 1.5707963267948966;
		else	     theta3 = (-0.69813170079773212 * (x3-2) * (x3-2) - 0.87266462599716477) * (x3-2) - 1.5707963267948966;

		/* */
		double SOD_00 = (theta1-f.theta1) * f.SO1[0][0] + (theta2-f.theta2) * f.SO2[0][0] + (theta3-f.theta3) * f.SO3[0][0];
		double SOD_01 = (theta1-f.theta1) * f.SO1[0][1] + (theta2-f.theta2) * f.SO2[0][1] + (theta3-f.theta3) * f.SO3[0][1];
		double SOD_10 = (theta1-f.theta1) * f.SO1[1][0] + (theta2-f.theta2) * f.SO2[1][0] + (theta3-f.theta3) * f.SO3[1][0];
		double SOD_11 = (theta1-f.theta1) * f.SO1[1][1] + (theta2-f.theta2) * f.SO2[1][1] + (theta3-f.theta3) * f.SO3[1][1];


			trace = (SOD_00 + SOD_11);
			trace_2 = SOD_00*SOD_00 + 2*SOD_01*SOD_10 + SOD_11*SOD_11;
			bending += weight *YOUNG/(1-POISSON*POISSON)*((1-POISSON)*trace_2+POISSON*trace*trace);
		/* -------------- compute gradient ------------------------- */
		/* --------------------------------------------------------- */

		for (int i = 1; i < 7; i++)
			for (int j = 1; j < 4; j++){
				/* */
				db1 = computeGradient(1,i,j);db2 = computeGradient(2,i,j);db3 = computeGradient(3,i,j);
				ds1 = computeGradient(4,i,j);ds2 = computeGradient(5,i,j);ds3 = computeGradient(6,i,j);

				dnm = computeCrossProductGradient(b1,b2,db1,db2);
				dn1 = computeCrossProductGradient(s1,b1,ds1,db1);
				dn2 = computeCrossProductGradient(s2,b2,ds2,db2);
				dn3 = computeCrossProductGradient(s3,b3,ds3,db3);

				dx1 = computeDotProductGradient(n1,nm,dn1,dnm);
				dx2 = computeDotProductGradient(n2,nm,dn2,dnm);
				dx3 = computeDotProductGradient(n3,nm,dn3,dnm);

				if (sign1 > 0) dx1 = -dx1;
				if (sign2 > 0) dx2 = -dx2;
				if (sign3 > 0) dx3 = -dx3;

				if (x1 <= 1) Dtheta1 = (-3 * 0.69813170079773212 * x1 * x1 - 0.87266462599716477)  * dx1;
				else		 Dtheta1 = (-3 * 0.69813170079773212 * (x1-2) * (x1-2) - 0.87266462599716477)  * dx1;
				if (x2 <= 1) Dtheta2 = (-3 * 0.69813170079773212 * x2 * x2 - 0.87266462599716477)  * dx2;
				else		 Dtheta2 = (-3 * 0.69813170079773212 * (x2-2) * (x2-2) - 0.87266462599716477)  * dx2;
				if (x3 <= 1) Dtheta3 = (-3 * 0.69813170079773212 * x3 * x3 - 0.87266462599716477)  * dx3;
				else		 Dtheta3 = (-3 * 0.69813170079773212 * (x3-2) * (x3-2) - 0.87266462599716477)  * dx3;

				DSOD_00 = Dtheta1 * f.SO1[0][0] + Dtheta2 * f.SO2[0][0] + Dtheta3 * f.SO3[0][0];
				DSOD_01 = Dtheta1 * f.SO1[0][1] + Dtheta2 * f.SO2[0][1] + Dtheta3 * f.SO3[0][1];
				DSOD_10 = Dtheta1 * f.SO1[1][0] + Dtheta2 * f.SO2[1][0] + Dtheta3 * f.SO3[1][0];
				DSOD_11 = Dtheta1 * f.SO1[1][1] + Dtheta2 * f.SO2[1][1] + Dtheta3 * f.SO3[1][1];

					Dbend = weight * YOUNG/(1-POISSON*POISSON)*
						((1-POISSON)*(2*SOD_00*DSOD_00+2*SOD_01*DSOD_10+2*DSOD_01*SOD_10+2*SOD_11*DSOD_11)+
						POISSON*2*trace*(DSOD_00+DSOD_11));
				//Dbend = weight * (2*SOD_00*DSOD_00 + 2*SOD_01*DSOD_01 + 2*SOD_10*DSOD_10 + 2*SOD_11*DSOD_11);

				if (i <= 3)
					g[f.index[i]*3+j-1] += Dbend;
				else if (i == 4)
					g[f.nbIdx1*3+j-1] += Dbend;
				else if (i == 5)
					g[f.nbIdx2*3+j-1] += Dbend;
				else if (i == 6)
					g[f.nbIdx3*3+j-1] += Dbend;
			}

	}

	facetBending = bending;
	return bending;
}

lbfgsfloatval_t penalizeBend(const lbfgsfloatval_t *u, lbfgsfloatval_t *g, int idx){
	double bending = 0;
	facet* faceList = Surface[idx]->faceList;

	Point_3 newV1,newV2,newV3,newNB1,newNB2,newNB3;
	double dv1x, dv1y, dv1z, dv2x, dv2y, dv2z, dv3x, dv3y, dv3z,
		   dnb1x,dnb1y,dnb1z,dnb2x,dnb2y,dnb2z,dnb3x,dnb3y,dnb3z;
	double Dtheta1,Dtheta2,Dtheta3,DSOD_00,DSOD_01,DSOD_10,DSOD_11;
	double weight,Dbend;

	for (int idx=0; idx < Surface[idx]->faceNum;idx++){
		facet f = faceList[idx];
		if (f.isBorder) continue;
		weight = BENDWEIGHT;

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

lbfgsfloatval_t penalizeLandmark(const lbfgsfloatval_t *u, lbfgsfloatval_t *g, int idx){
	PolyhedralSurf::Vertex_const_handle v1;
	double weight = LANDMARKWEIGHT;
	lbfgsfloatval_t fx = 0.0;

	int startIdx = max_zenyo(0,idx - min_zenyo(GEOMETRIC_RANGE,TEXTURE_RANGE));
	int endIdx = min_zenyo(surfaceNum-1,idx + min_zenyo(GEOMETRIC_RANGE,TEXTURE_RANGE));

	for (int i = startIdx; i <= endIdx; i++){
		if (i == idx) continue;

		for (int j = 0; j < Surface[idx]->landmarkNum[i]; j++){

			int i1 = Surface[idx]->landmarks1[i][j];
			int i2 = Surface[idx]->landmarks2[i][j];

			v1 = Surface[idx]->vertexIndex[i1];
			Point_3 deformed = v1->point() + Vector_3(u[i1*3],u[i1*3+1],u[i1*3+2]);
			Vector_3 dis(deformed,Surface[i]->vertexIndex[i2]->point());

			/*
			fx = fx + weight*dis.squared_length();
			g[i1*3] += -2 * dis.x() * weight;
			g[i1*3+1] += -2 * dis.y() * weight;
			g[i1*3+2] += -2 * dis.z() * weight;
			*/
			double len = sqrt(dis.squared_length());
			if (len<=3){
				g[i1*3]   += - dis.x() * weight;
				g[i1*3+1] += - dis.y() * weight;
				g[i1*3+2] += - dis.z() * weight;
				fx = fx + weight * 0.5 * dis.squared_length();
			}else{
				g[i1*3]   += 3 / (2 * len) * (2 * -dis.x()) * weight;
				g[i1*3+1] += 3 / (2 * len) * (2 * -dis.y()) * weight;
				g[i1*3+2] += 3 / (2 * len) * (2 * -dis.z()) * weight;
				fx = fx + weight * 3 * (len - 1.5);
			}
		}
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

inline threeTuple computeGradient(int v,int i,int j){
	if (v <= 3){ // computing derivatives w.r.t. b1,b2,b3
		if (i>3) return threeTuple(0,0,0);

		if (v == i){
			if (j == 1) return threeTuple(-1,0,0);
			if (j == 2) return threeTuple(0,-1,0);
			if (j == 3) return threeTuple(0,0,-1);
		}

		if ((v+1 == i) || (v-2 == i)){
			if (j == 1) return threeTuple(1,0,0);
			if (j == 2) return threeTuple(0,1,0);
			if (j == 3) return threeTuple(0,0,1);
		}

		return threeTuple(0,0,0);
	}

	if (v-3 == i){
		if (j == 1) return threeTuple(-1,0,0);
		if (j == 2) return threeTuple(0,-1,0);
		if (j == 3) return threeTuple(0,0,-1);
	}

	if (v == i){
		if (j == 1) return threeTuple(1,0,0);
		if (j == 2) return threeTuple(0,1,0);
		if (j == 3) return threeTuple(0,0,1);
	}

	return threeTuple(0,0,0);
}

inline threeTuple computeCrossProductGradient(Vector_3 b1,Vector_3 b2,threeTuple db1,threeTuple db2){
	double dx = (db1.y*b2.z()+b1.y()*db2.z)-(db1.z*b2.y()+b1.z()*db2.y);
	double dy = (db1.z*b2.x()+b1.z()*db2.x)-(db1.x*b2.z()+b1.x()*db2.z);
	double dz = (db1.x*b2.y()+b1.x()*db2.y)-(db1.y*b2.x()+b1.y()*db2.x);

	return threeTuple(dx,dy,dz);
}

inline double computeDotProductGradient(Vector_3 n1,Vector_3 n2,threeTuple dn1,threeTuple dn2){
	double g = n1*n2;
	double dg = dn1.x*n2.x() + n1.x()*dn2.x +
		dn1.y*n2.y() + n1.y()*dn2.y +
		dn1.z*n2.z() + n1.z()*dn2.z;

	double N1 = n1.squared_length() + EPS;
	double N2 = n2.squared_length() + EPS;
	double N1N2 = N1*N2;

	double dN1 = 2*n1.x()*dn1.x + 2*n1.y()*dn1.y + 2*n1.z()*dn1.z;
	double dN2 = 2*n2.x()*dn2.x + 2*n2.y()*dn2.y + 2*n2.z()*dn2.z;

	double dN1N2 = dN1*N2 + N1*dN2;

	double dSN1N2 = 0.5*dN1N2/sqrt(N1N2);

	return (sqrt(N1N2)*dg - dSN1N2*g)/(N1N2);
}

