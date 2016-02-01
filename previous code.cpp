#include "stdafx.h"
//select BESTMATCH correspondences
	/*
	for (int k = 0; k < BESTMATCH; k++){
		float maxA = -1;
		int maxI,maxJ;
		for (int i = 0; i < affinityM; i++){
			if (matchWeight[i] > 0.001) continue;

			for (int j = 0; j < affinityN; j++)
				if (VAL(affinity,i,j,affinityN) > maxA){
					maxA = VAL(affinity,i,j,affinityN);
					maxI = i;
					maxJ = j;
				}
		}

		bestMatch[maxI] = maxJ;
		matchWeight[maxI] = VAL(affinity,maxI,maxJ,affinityN);

		for (int i = 0; i < affinityM; i++)
			VAL(affinity,i,maxJ,affinityN) = 0;

	}
	

	meshColor.clear();
	for (int i = 0; i < affinityM; i++)
		if (matchWeight[i] > 0.00001){
			meshColor.insert(pair<Vertex_const_handle,float>(vertexIndex[i],1));
		}
	
	std::map<PolyhedralSurf::Vertex_const_handle, float>::iterator iter;

	for (int i = 0; i < faceNum; i++){
		iter = meshColor.find(faceList[i].v1);
		if (iter != meshColor.end())
			faceList[i].color1 = threeTuple(0,1,0);
		else 
			faceList[i].color1 = threeTuple(0.8,0.8,0.8);

		iter = meshColor.find(faceList[i].v2);
		if (iter != meshColor.end())
			faceList[i].color2 = threeTuple(0,1,0);
		else 
			faceList[i].color2 = threeTuple(0.8,0.8,0.8);

		iter = meshColor.find(faceList[i].v3);
		if (iter != meshColor.end())
			faceList[i].color3 = threeTuple(0,1,0);
		else 
			faceList[i].color3 = threeTuple(0.8,0.8,0.8);

	}
	*/
/*test convex */
// 		readOFF("temp1.off",tempU);
// 		PolyhedralSurf::Vertex_const_iterator vb = RCSurface->P.vertices_begin();
// 		PolyhedralSurf::Vertex_const_iterator ve = RCSurface->P.vertices_end();
// 		PolyhedralSurf::Vertex_const_handle vh;
// 
// 		for (int i = 0; vb != ve; vb++,i++){
// 			vh = vb;
// 			tempU[i*3] = tempU[i*3] - vh->point().x();
// 			tempU[i*3+1] = tempU[i*3+1] - vh->point().y();
// 			tempU[i*3+2] = tempU[i*3+2] - vh->point().z();
// 		}

// 		for (int i = 0; i < affinityM * 3; i++)
// 			tempU[i] = float(rand() % 100) / 200;

/* register to point cloud*/
// 		targetPos = new float[affinityM*3];
// 		memset(targetPos,0,sizeof(float)*affinityM*3);
// 		ifstream fin;
// 		fin.open("data/phantom/pt.txt");
// 		int targetNum;
// 		fin>>targetNum;
// 		cout<<"targetNum: "<<targetNum<<endl;
// 		for (int i = 0; i < targetNum; i++){
// 			int idx;
// 			float x,y,z;
// 			fin>>idx>>x>>y>>z;
// 			noStretch[idx] = true;
// 			targetPos[idx*3]   = x;
// 			targetPos[idx*3+1] = y;
// 			targetPos[idx*3+2] = z;
// 
// 		}
// 		fin.close();
// 		cout<<"targetNum: "<<targetNum<<endl;




//bending angle
/*
	for (int i=0; i < RCSurface->edgeNum/2;i++){
		edge he = RCSurface->edgeList[i];
		double weight = he.bendWeight * bendweight;
		
		Point_3 newV1 = he.v1->point() + Vector_3(u[he.index1*3],u[he.index1*3+1],u[he.index1*3+2]);
		Point_3 newV2 = he.v2->point() + Vector_3(u[he.index2*3],u[he.index2*3+1],u[he.index2*3+2]);
		Point_3 newVl = he.vl->point() + Vector_3(u[he.indexl*3],u[he.indexl*3+1],u[he.indexl*3+2]);
		Point_3 newVr = he.vr->point() + Vector_3(u[he.indexr*3],u[he.indexr*3+1],u[he.indexr*3+2]);

		Vector_3 v1(newV1,newVr);
		Vector_3 v2(newV1,newV2);
		Vector_3 v3(newV1,newVl);

		Vector_3 n1 = cross_product(v1,v2);
		Vector_3 n2 = cross_product(v2,v3);

		double sm1 = n1.squared_length();
		double sq1 = sqrt(sm1);
		double sm2 = n2.squared_length();
		double sq2 = sqrt(sm2);

		double p = n1*n2;
		double gg = sq1*sq2;

		double kos = p/gg;
		
		if (kos > 1) kos = 1;
		if (kos < -1) kos = -1;

		double theta = acos(kos);
		
		double w = (theta - he.angle)*(theta - he.angle);

		bending += w * weight;

		//compute bending gradient
		double dn1x,dn1y,dn1z,dn2x,dn2y,dn2z,dsm1,dsm2,dsq1,dsq2,dp,dg,dk,dtheta,dw;
		double dv1x,dv1y,dv1z,dv2x,dv2y,dv2z,dv3x,dv3y,dv3z;

		for (int j = 1; j <= 4; j++)
			for (int k = 1; k <= 3; k++){
				dv1x = computeGradient(1,1,j,k);dv1y = computeGradient(1,2,j,k); dv1z = computeGradient(1,3,j,k);
				dv2x = computeGradient(2,1,j,k);dv2y = computeGradient(2,2,j,k); dv2z = computeGradient(2,3,j,k);
				dv3x = computeGradient(3,1,j,k);dv3y = computeGradient(3,2,j,k); dv3z = computeGradient(3,3,j,k);

				dn1x = (dv1y*v2.z()+v1.y()*dv2z)-(dv1z*v2.y()+v1.z()*dv2y);
				dn1y = (dv1z*v2.x()+v1.z()*dv2x)-(dv1x*v2.z()+v1.x()*dv2z);
				dn1z = (dv1x*v2.y()+v1.x()*dv2y)-(dv1y*v2.x()+v1.y()*dv2x);
				dn2x = (dv2y*v3.z()+v2.y()*dv3z)-(dv2z*v3.y()+v2.z()*dv3y);
				dn2y = (dv2z*v3.x()+v2.z()*dv3x)-(dv2x*v3.z()+v2.x()*dv3z);
				dn2z = (dv2x*v3.y()+v2.x()*dv3y)-(dv2y*v3.x()+v2.y()*dv3x);

				dsm1 = 2*n1.x()*dn1x + 2*n1.y()*dn1y + 2*n1.z()*dn1z;
				dsm2 = 2*n2.x()*dn2x + 2*n2.y()*dn2y + 2*n2.z()*dn2z;

				dsq1 = 0.5 * 1 / sqrt(sm1) * dsm1;
				dsq2 = 0.5 * 1 / sqrt(sm2) * dsm2;

				dg = dsq1*sq2+sq1*dsq2;
				dp = (dn1x*n2.x()+n1.x()*dn2x) + (dn1y*n2.y()+n1.y()*dn2y) + (dn1z*n2.z()+n1.z()*dn2z);

				dk = (dp*gg-p*dg)/(gg*gg);

				if (abs(kos) == 1)
					dtheta = 0;
				else
					dtheta = -1 / sqrt(1-kos*kos)*dk;

				dw = 2*(theta - he.angle)*dtheta;

				g[he.index[j]*3+k-1] += dw * weight;
			}
		
	}
*/

//double computeGradient(int n,int m,int x,int y);

/*
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
*/

//bending determinant
/*
double bending = 0;

double yz1,xz1,xy1,xyz,det;
double dyz1,dxz1,dxy1,dxyz,ddet;
for (int i=0; i < RCSurface->edgeNum/2;i++){
	if (RCSurface->edgeList[i].isBoundary) continue;

	edge he = RCSurface->edgeList[i];
	double weight = he.bendWeight * bendweight;

	Point_3 v1 = he.v1->point() + Vector_3(u[he.index1*3],u[he.index1*3+1],u[he.index1*3+2]);
	Point_3 v2 = he.v2->point() + Vector_3(u[he.index2*3],u[he.index2*3+1],u[he.index2*3+2]);
	Point_3 vl = he.vl->point() + Vector_3(u[he.indexl*3],u[he.indexl*3+1],u[he.indexl*3+2]);
	Point_3 vr = he.vr->point() + Vector_3(u[he.indexr*3],u[he.indexr*3+1],u[he.indexr*3+2]);


	yz1 = v2.y()*vr.z()+v2.z()*vl.y()+vr.y()*vl.z()-vr.z()*vl.y()-v2.z()*vr.y()-v2.y()*vl.z();
	xz1 = v2.x()*vr.z()+v2.z()*vl.x()+vr.x()*vl.z()-vr.z()*vl.x()-v2.z()*vr.x()-v2.x()*vl.z();
	xy1 = v2.x()*vr.y()+v2.y()*vl.x()+vr.x()*vl.y()-vr.y()*vl.x()-v2.y()*vr.x()-v2.x()*vl.y();
	xyz = v2.x()*vr.y()*vl.z()+v2.y()*vr.z()*vl.x()+vr.x()*vl.y()*v2.z()-v2.z()*vr.y()*vl.x()-v2.y()*vr.x()*vl.z()-v2.x()*vr.z()*vl.y();

	double det = v1.x()*yz1-v1.y()*xz1+v1.z()*xy1-xyz;

	bending += (det-he.angle) * (det-he.angle) * weight;

	//compute gradient
	//dv1x
	dyz1 = 0; dxz1 = 0; dxy1 = 0; dxyz = 0; ddet = yz1;
	g[he.index1*3] += 2 * (det-he.angle) * ddet * weight;
	//dv1y
	dyz1 = 0; dxz1 = 0; dxy1 = 0; dxyz = 0; ddet = -xz1;
	g[he.index1*3+1] += 2 * (det-he.angle) * ddet * weight;
	//dv1z
	dyz1 = 0; dxz1 = 0; dxy1 = 0; dxyz = 0; ddet = xy1;
	g[he.index1*3+2] += 2 * (det-he.angle) * ddet * weight;

	//dv2x
	dyz1 = 0; dxz1 = vr.z()-vl.z(); dxy1 = vr.y()-vl.y(); dxyz = vr.y()*vl.z()-vr.z()*vl.y(); ddet = v1.x()*dyz1-v1.y()*dxz1+v1.z()*dxy1-dxyz;
	g[he.index2*3] += 2 * (det-he.angle) * ddet * weight;
	//dv2y
	dyz1 = vr.z()-vl.z(); dxz1 = 0; dxy1 = vl.x()-vr.x(); dxyz = vr.z()*vl.x()-vr.x()*vl.z(); ddet = v1.x()*dyz1-v1.y()*dxz1+v1.z()*dxy1-dxyz;
	g[he.index2*3+1] += 2 * (det-he.angle) * ddet * weight;
	//dv2z
	dyz1 = vl.y()-vr.y(); dxz1 = vl.x()-vr.x(); dxy1 = 0; dxyz = vr.x()*vl.y()-vr.y()*vl.x(); ddet = v1.x()*dyz1-v1.y()*dxz1+v1.z()*dxy1-dxyz;
	g[he.index2*3+2] += 2 * (det-he.angle) * ddet * weight;

	//dvlx
	dyz1 = 0; dxz1 = v2.z()-vr.z(); dxy1 = v2.y()-vr.y(); dxyz = v2.y()*vr.z()-v2.z()*vr.y(); ddet = v1.x()*dyz1-v1.y()*dxz1+v1.z()*dxy1-dxyz;		
	g[he.indexl*3] += 2 * (det-he.angle) * ddet * weight;
	//dvly
	dyz1 = v2.z()-vr.z(); dxz1 = 0; dxy1 = vr.x()-v2.x(); dxyz = v2.z()*vr.x()-v2.x()*vr.z(); ddet = v1.x()*dyz1-v1.y()*dxz1+v1.z()*dxy1-dxyz;
	g[he.indexl*3+1] += 2 * (det-he.angle) * ddet * weight;
	//dvlz
	dyz1 = vr.y()-v2.y(); dxz1 = vr.x()-v2.x(); dxy1 = 0; dxyz = v2.x()*vr.y()-v2.y()*vr.x(); ddet = v1.x()*dyz1-v1.y()*dxz1+v1.z()*dxy1-dxyz;
	g[he.indexl*3+2] += 2 * (det-he.angle) * ddet * weight;

	//dvrx
	dyz1 = 0; dxz1 = vl.z()-v2.z(); dxy1 = vl.y()-v2.y(); dxyz = vl.y()*v2.z()-vl.z()*v2.y(); ddet = v1.x()*dyz1-v1.y()*dxz1+v1.z()*dxy1-dxyz;
	g[he.indexr*3] += 2 * (det-he.angle) * ddet * weight;
	//dvry
	dyz1 = vl.z()-v2.z(); dxz1 = 0; dxy1 = v2.x()-vl.x(); dxyz = vl.z()*v2.x()-vl.x()*v2.z(); ddet = v1.x()*dyz1-v1.y()*dxz1+v1.z()*dxy1-dxyz;
	g[he.indexr*3+1] += 2 * (det-he.angle) * ddet * weight;
	//dvrz
	dyz1 = v2.y()-vl.y(); dxz1 = v2.x()-vl.x(); dxy1 = 0; dxyz = vl.x()*v2.y()-vl.y()*v2.x(); ddet = v1.x()*dyz1-v1.y()*dxz1+v1.z()*dxy1-dxyz;
	g[he.indexr*3+2] += 2 * (det-he.angle) * ddet * weight;
}

return bending;
*/


// void BasicMesh::ComputeAdjacency(){
// 	memset(adjacent,false,affinityM * affinityM * sizeof(bool));
// 
// 	for (int i = 0; i < vertexNum; i++){
// 		if (!attractor[i]) continue;
// 
// 		priority_queue<node> dijkstra;
// 
// 		PolyhedralSurf::Vertex_const_handle vh = vertexIndex[i];
// 		dijkstra.push(node(vh,0));
// 
// 		while (true){
// 			node curNode;
// 			while (true){
// 				curNode = dijkstra.top();
// 				dijkstra.pop();
// 
// 				int idx = indexMap.find(curNode.vh)->second;
// 				if (!VAL(adjacent,i,idx,affinityM)) break;
// 			}
// 
// 			if (curNode.dis > ADJACENTDIS) break;
// 
// 			int idx = indexMap.find(curNode.vh)->second;
// 			VAL(adjacent,i,idx,affinityM) = true;
// 
// 			PolyhedralSurf::Halfedge_around_vertex_const_circulator tempHfe = curNode.vh->vertex_begin();
// 
// 			int v1Degree = curNode.vh->degree();
// 			int degreeCount = 1;
// 			while (degreeCount<=v1Degree){
// 				PolyhedralSurf::Vertex_const_handle vnext = tempHfe->opposite()->vertex();
// 
// 				if (!VAL(adjacent,i,indexMap.find(vnext)->second,affinityM))
// 					dijkstra.push(node(vnext,curNode.dis + computeLength(tempHfe)));
// 
// 				tempHfe++;
// 				degreeCount++;
// 			}	
// 
// 		}
// 
// 		while (!dijkstra.empty()) dijkstra.pop();
// 
// 	}
// }

/* backwards reg */

// void BasicMesh::computeClosestPtDis(BasicMesh* secondMesh){
// 	PolyhedralSurf::Vertex_const_iterator vb = P.vertices_begin();
// 	PolyhedralSurf::Vertex_const_iterator ve = P.vertices_end();
// 
// 	PolyhedralSurf::Vertex_const_handle vh1,vh2;
// 
// 	for (int i = 0; vb != ve; vb++,i++){
// 		//cout<<"Working on vertex: "<<i<<endl;
// 		vh1 = vb;
// 
// 		for (int j = 0; j < secondMesh->vertexNum; j++){
// 			vh2 = secondMesh->vertexIndex[j];
// 
// 			VAL(closestDis,i,j,affinityN) = computeEuclideanDis(vh1->point(),vh2->point());
// 		}
// 	}
// }
// 
// void BasicMesh::findCorrespondenceBackwards(BasicMesh* secondMesh,int disWeight,BasicMesh* DFMesh){
// 	for (int i = 0; i < affinityM; i++)
// 		for (int j = 0; j < affinityN; j++)
// 			VAL(affinity,i,j,affinityN) = -1;
// 
// 
// 	PolyhedralSurf::Vertex_const_iterator vb = P.vertices_begin();
// 	PolyhedralSurf::Vertex_const_iterator ve = P.vertices_end();
// 
// 	PolyhedralSurf::Vertex_const_handle vh1,vh2,vh3;
// 	float dis,dif;
// 
// 	for (int i = 0; vb != ve; vb++,i++){
// 		//cout<<"Working on vertex: "<<i<<endl;
// 
// 		vh1 = vb;
// 
// 		for (int j = 0; j < secondMesh->vertexNum; j++){
// 			vh2 = secondMesh->vertexIndex[j];
// 
// 			dif = VAL(closestDis,i,j,affinityN);
// 			dis = computeEuclideanDis(vh1->point(),vh2->point());
// 
// 			dif = dif * (1-disWeight) + dis * disWeight;
// 
// 			VAL(affinity,i,j,affinityN) = dif;
// 		}
// 	}
// 
// 	for (int i = 0; i < affinityM; i++){
// 
// 		for (int j = 0; j < affinityN; j++){
// 			double scaledWeight;
// 
// 			if (VAL(affinity,i,j,affinityN) >= 0)
// 				scaledWeight = -VAL(affinity,i,j,affinityN)*0.75;
// 			else 
// 				scaledWeight = -100001;
// 
// 			if (scaledWeight < -100000)
// 				VAL(affinity,i,j,affinityN) = 0;
// 			else 
// 				VAL(affinity,i,j,affinityN) = exp(scaledWeight);
// 		}
// 
// 	}		
// }
// 
// void BasicMesh::findMatchBackwards(BasicMesh* secondMesh){
// 	//use all correspondences
// 	for (int i = 0; i < affinityM; i++){
// 
// 		float maxA = -1;
// 		for (int j = 0; j < affinityN; j++)
// 			if (VAL(affinity,i,j,affinityN) > maxA){
// 				maxA = VAL(affinity,i,j,affinityN);
// 				bestMatch[i] = j;
// 				matchWeight[i] = VAL(affinity,i,j,affinityN);
// 			}
// 
// 	}
// 
// }

/* optimization cmaes */
// void startCMAES();
// double evaluate_cmaes(double const *u, int N);
// double penalizeData_cmaes(double const *u, int N);
// double penalizeStretch_cmaes(double const *u, int N);
// double penalizeBend_cmaes(double const *u, int N);
// 
// /* optimization cs */
// void startOptimization_cg();
// void evaluate_cg(const alglib::real_1d_array &u, double &func, alglib::real_1d_array &g, void *ptr);
// double penalizeData(const alglib::real_1d_array &u, alglib::real_1d_array &g);
// double penalizeStretch(const alglib::real_1d_array &u, alglib::real_1d_array &g);
// double penalizeBend(const alglib::real_1d_array &u, alglib::real_1d_array &g);

// void startOptimization(){
// 	/* initialize optimization */
// 	char filename[50];
// 
// 	lbfgsfloatval_t* tempU = lbfgs_malloc(affinityM * 3);
// 	int ret;
// 
// 	double bestError = 100000;
// 	QueryPerformanceCounter(&t1);
// 	for (int i = 0; i < 9; i++){
// 		cout<<"Velocity Iteration: "<<i<<endl;
// 
// 		if (i > 0){
// 			cout<<"Reset surface..."<<endl;
// 
// 			RCSurface->ComputeMeshProperty(filename);
// 
// 			/* output errer */
// 			float finalErr,finalVar,finalBoundErr,finalBoundVar,maxErr;
// 			computeStatisticBool(RCSurface,CTSurface,finalErr,finalBoundErr,finalVar,finalBoundVar,maxErr);
// 
// 			if (finalBoundErr < bestError){
// 				bestError = finalBoundErr;
// 
// 				ofstream fout;
// 				fout.open("error.txt");
// 
// 				fout<<initialError<<' '<<initialVariance<<endl;
// 				fout<<initialBoundError<<' '<<initialBoundVar<<endl;
// 
// 
// 				fout<<finalErr<<' '<<finalVar<<endl;
// 				fout<<finalBoundErr<<' '<<finalBoundVar<<endl;
// 
// 				fout.close();
// 			}
// 
// 			/* output errer */
// 
// 			cout<<"Recompute geometric features..."<<endl;
// 			RCSurface->findSignatureAll();
// 
// 			memset(u,0,sizeof(lbfgsfloatval_t) * affinityM * 3);
// 			memset(affinity,0,affinityM * affinityN * sizeof(double));
// 			memset(bestMatch,0,sizeof(int)*affinityM);
// 			memset(matchWeight,0,sizeof(double)*affinityM);
// 
// 			cout<<"Recompute affinity map..."<<endl;
// 
// 			//RCSurface->constructEdgeList();
// 			/* prepare registration */
// 			RCSurface->constructLink();
// 			RCSurface->findCorrespondenceBothWay(CTSurface,0.9/*min_zenyo(i * 0.13,1)*/);
// 			RCSurface->findMatch(CTSurface);
// 			//RCSurface->findClosestPoint(CTSurface);
// 		}
// 
// 		cout<<"Begin optimization..."<<endl;
// 
// 		lbfgsfloatval_t lastFx = 100000001;
// 		lbfgsfloatval_t fx = 100000000;
// 
// 		int alterCount = 0;
// 		memset(tempU,0,sizeof(lbfgsfloatval_t)*affinityM*3);
// 
// 		/* register to point cloud*/
// 
// 		int opIterNum = 80000;
// 
// 		while (lastFx > fx + 0.01){
// 			memcpy(u,tempU,sizeof(lbfgsfloatval_t)*affinityM*3);
// 			lastFx = fx;
// 
// 			if (alterCount == 1) break;
// 			alterCount++;
// 			cout<<"--------------------------zig iteration: "<<alterCount<<endl;
// 
// 			// 			lbfgs_parameter_t param;
// 			// 			regweight = REGWEIGHT;
// 			// 			bendweight = 0;
// 			// 			ret = 0;
// 			// 
// 			// 			lbfgs_parameter_init(&param);
// 			// 			param.max_iterations = 10;
// 			// 			ret = lbfgs(affinityM * 3, tempU, &fx, evaluate, progress, NULL, &param);
// 			// 			
// 			// 			cout<<"function value: "<<fx<<endl;
// 
// 			////
// 
// 			cout<<"-------------------------zag iteration: "<<alterCount<<endl;
// 			lbfgs_parameter_t param1;
// 			regweight = REGWEIGHT * (1 - min_zenyo(i * 0.1,1));
// 			bendweight = BENDWEIGHT * (1 - min_zenyo(i * 0.1,1));
// 			ret = 0;
// 
// 			lbfgs_parameter_init(&param1);
// 			param1.max_iterations = opIterNum;
// 			ret = lbfgs(affinityM * 3, tempU, &fx, evaluate, progress, NULL, &param1);
// 
// 			cout<<"function value: "<<fx<<endl;
// 		}
// 
// 
// 		//output file
// 		filename[0] = '\0';
// 		strcat(filename,"temp");
// 		char iterNum[10];
// 		itoa(i,iterNum,10);
// 		strcat(filename,iterNum);
// 		strcat(filename,".off");
// 		writeOFF(RCSurface,filename);
// 
// 		cout<<"L-BFGS optimization terminated with status code: "<<ret<<"fx: "<<lastFx<<endl;
// 
// 	}
// 
// 	RCSurface->ComputeMeshProperty(filename);
// 	float finalErr,finalVar,finalBoundErr,finalBoundVar,maxErr;
// 	computeStatisticBool(RCSurface,CTSurface,finalErr,finalBoundErr,finalVar,finalBoundVar,maxErr);
// 
// 	if (finalBoundErr < bestError){
// 		bestError = finalBoundErr;
// 
// 		ofstream fout;
// 		fout.open("error.txt");
// 
// 		fout<<initialError<<' '<<initialVariance<<endl;
// 		fout<<initialBoundError<<' '<<initialBoundVar<<endl;
// 
// 
// 		fout<<finalErr<<' '<<finalVar<<endl;
// 		fout<<finalBoundErr<<' '<<finalBoundVar<<endl;
// 
// 		fout.close();
// 	}
// }