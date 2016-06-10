#include "stdafx.h"

#ifdef CGAL_USE_BOOST_PROGRAM_OPTIONS
#include <boost/program_options.hpp>
namespace po = boost::program_options;
#endif

#include <iostream>
#include <fstream>
using namespace std;

#include "MeshObject.h"


/* gather points around the vertex v using rings on the
   polyhedralsurf. the collection of points resorts to 3 alternatives:
   1. the exact number of points to be used
   2. the exact number of rings to be used
   3. nothing is specified
*/
void BasicMesh::gather_fitting_points(Vertex_const_handle v,std::vector<Point_3> &in_points, Poly_rings& poly_rings)
{
	//container to collect vertices of v on the PolyhedralSurf
	std::vector<Vertex_const_handle> gathered;
	//initialize
	in_points.clear();

	//OPTION -p nb_points_to_use, with nb_points_to_use != 0. Collect
	//enough rings and discard some points of the last collected ring to
	//get the exact "nb_points_to_use"
	if ( nb_points_to_use != 0 ) {
		poly_rings.collect_enough_rings(v, nb_points_to_use, gathered);//, vpm);
		if ( gathered.size() > nb_points_to_use ) gathered.resize(nb_points_to_use);
	}
	else {	// nb_points_to_use=0, this is the default and the option -p is not considered;
			// then option -a nb_rings is checked. If nb_rings=0, collect
			// enough rings to get the min_nb_points required for the fitting
			// else collect the nb_rings required
    if ( nb_rings == 0 )
		poly_rings.collect_enough_rings(v, min_nb_points, gathered);//, vpm);
    else poly_rings.collect_i_rings(v, nb_rings, gathered);//, vpm);
	}

	//store the gathered points
	std::vector<Vertex_const_handle>::const_iterator
    itb = gathered.begin(), ite = gathered.end();
	CGAL_For_all(itb,ite) in_points.push_back((*itb)->point());
}

/* Use the jet_fitting package and the class Poly_rings to compute
   diff quantities.
*/
void BasicMesh::compute_differential_quantities(PolyhedralSurf& P, Poly_rings& poly_rings)
{
	//container for approximation points
	std::vector<Point_3> in_points;

	//MAIN LOOP
	Vertex_const_iterator vitb = P.vertices_begin(), vite = P.vertices_end();

	int i = 0;
	for (; vitb != vite; vitb++,i++) {

		//initialize
		Vertex_const_handle v = vitb;
		
		in_points.clear();
		Monge_form monge_form;
		Monge_via_jet_fitting monge_fit;

		//gather points around the vertex using rings
		gather_fitting_points(v, in_points, poly_rings);

		//exit if the nb of points is too small
		if ( in_points.size() < min_nb_points )
		{
			cout<<v->point()<<endl;
			std::cerr << "Too few points to perform the fitting" << std::endl; 
			exit(1);
		}

		//For Ridges we need at least 3rd order info
		assert( d_monge >= 3);
		// run the main fct : perform the fitting
		monge_form = monge_fit(in_points.begin(), in_points.end(),d_fitting, d_monge);

		//switch min-max ppal curv/dir wrt the mesh orientation
		const Vector_3 normal_mesh = P.computeFacetsAverageUnitNormal(v);
		monge_form.comply_wrt_given_normal(normal_mesh);

		//Store monge data needed for ridge computations in property maps
		vertex2d1_map[v] = monge_form.maximal_principal_direction();
		vertex2d2_map[v] = monge_form.minimal_principal_direction();
		vertex2k1_map[v] = monge_form.coefficients()[0];
		vertex2k2_map[v] = monge_form.coefficients()[1];
		vertex2b0_map[v] = monge_form.coefficients()[2];
		vertex2b3_map[v] = monge_form.coefficients()[5];
		if ( d_monge >= 4) {
			//= 3*b1^2+(k1-k2)(c0-3k1^3)
		vertex2P1_map[v] =
			3*monge_form.coefficients()[3]*monge_form.coefficients()[3]
			+(monge_form.coefficients()[0]-monge_form.coefficients()[1])
			*(monge_form.coefficients()[6]
			-3*monge_form.coefficients()[0]*monge_form.coefficients()[0]
			*monge_form.coefficients()[0]);
		//= 3*b2^2+(k2-k1)(c4-3k2^3)
		vertex2P2_map[v] =
			3*monge_form.coefficients()[4]*monge_form.coefficients()[4]
			+(-monge_form.coefficients()[0]+monge_form.coefficients()[1])
			*(monge_form.coefficients()[10]
			-3*monge_form.coefficients()[1]*monge_form.coefficients()[1]
			*monge_form.coefficients()[1]);
		}
	} //END FOR LOOP
}

BasicMesh::BasicMesh(string filename,int idx){
	// additional parameters
	minx = 10000;
	miny = 10000;
	minz = 10000;

	maxx = -10000;
	maxy = -10000;
	maxz = -10000;

	//for CGAL computing
	d_fitting = 3;
	d_monge = 3;

	nb_rings = 0;
	nb_points_to_use = 0;

	int_tag = 3;
	tag_order = CGAL::Ridge_order_3;
	umb_size = 2;

	min_nb_points = (d_fitting + 1) * (d_fitting + 2) / 2;

	// file name
	if_name = filename;
	this->idx = idx;

	for (int i = 0; i < 8; i++){
		marchList[i] = NULL;
	}

	vertexIndex = NULL;
	faceList = NULL;
	edgeList = NULL;
	vertexList = NULL;

	u = NULL;
	bestMatch = NULL;
	matchWeight = NULL;
	targetPos = NULL;
	affinity = NULL;
	overallWeight = NULL;
	landmarkNum = NULL;
	landmarks1 = NULL;
	landmarks2 = NULL;

	dirNum = 8;
}

/* this function is borrowed for CGAL. however, translateMesh() is our own function for processing surfaces */
int BasicMesh::ComputeMeshProperty(){
	Vertex2FT_property_map vertex2k1_pm(vertex2k1_map), vertex2k2_pm(vertex2k2_map),
		vertex2b0_pm(vertex2b0_map), vertex2b3_pm(vertex2b3_map),
		vertex2P1_pm(vertex2P1_map), vertex2P2_pm(vertex2P2_map),vertex2P2_p();
	Vertex2Vector_property_map vertex2d1_pm(vertex2d1_map), vertex2d2_pm(vertex2d2_map);

	std::ifstream stream(if_name.c_str());
	P.clear();
	stream >> P;

	cout<<"vertex number: "<<P.size_of_vertices()<<endl;

	P.compute_facets_normals();

	//create a Poly_rings object
	Poly_rings poly_rings(P);
	compute_differential_quantities(P, poly_rings);

	//construct vertex<-->index mapping data
	indexMap.clear();
	signatureMap.clear();
	vertexBool.clear();
	meshColor.clear();
	if (vertexIndex != NULL) delete []vertexIndex;

	vertexIndex = new Vertex_const_handle[P.size_of_vertices()];

	PolyhedralSurf::Vertex_const_iterator vb = P.vertices_begin();
	PolyhedralSurf::Vertex_const_iterator ve = P.vertices_end();
	PolyhedralSurf::Vertex_const_handle vh,minh;

	for (int i = 0; vb != ve; vb++,i++){
		vh = vb;
		indexMap.insert(pair<Vertex_const_handle,int>(vh,i));
		vertexIndex[i] = vh;
	}
	
	/* translate surfaces (per faceet) */
	translateMesh(vertex2k1_pm,vertex2k2_pm);

	return 1;
}

void BasicMesh::ComputeStretchInverse(Vertex_const_handle v1,Vertex_const_handle v2,Vertex_const_handle v3,facet* face){
	Vector_3 temp_v1(v1->point(),v2->point());
	Vector_3 temp_v2(v1->point(),v3->point());

	double temp_v1_1 = 0;
	double temp_v1_2 = sqrt(temp_v1.squared_length()+EPS);

	double temp_v2_2 = temp_v1*temp_v2/temp_v1_2;

	double sqv21 = temp_v2.squared_length() - temp_v2_2*temp_v2_2;
	double temp_v2_1 = - sqrt(sqv21 + EPS);

	//local parametrization of vectors
	face->local_v1 = Vector_2(temp_v1_1,temp_v1_2);
	face->local_v2 = Vector_2(temp_v2_1,temp_v2_2-temp_v1_2);
	face->local_v3 = Vector_2(-temp_v2_1,-temp_v2_2);

	face->t1 = Vector_2(temp_v1_2,-temp_v1_1);
	face->t2 = Vector_2(temp_v2_2-temp_v1_2,-temp_v2_1);
	face->t3 = Vector_2(-temp_v2_2,temp_v2_1);

	//stretch inverse matrix
	double temp_scale = 1/(temp_v1_1*temp_v2_2-temp_v1_2*temp_v2_1);
	face->inverse[0][0] = temp_v2_2 * temp_scale;
	face->inverse[0][1] = -temp_v2_1 * temp_scale;
	face->inverse[1][0] = -temp_v1_2 * temp_scale;
	face->inverse[1][1] = temp_v1_1 * temp_scale;
}

void BasicMesh::ComputeShapeOperator(facet* face){
	/* determinant based */
// 	face->theta1 = -computeDeterminant(face->v1,face->v2,face->v3,face->nb1) / ((face->area*face->sideArea1)/face->l1);
// 	face->theta2 = -computeDeterminant(face->v1,face->v2,face->v3,face->nb2) / ((face->area*face->sideArea2)/face->l2);
// 	face->theta3 = -computeDeterminant(face->v1,face->v2,face->v3,face->nb3) / ((face->area*face->sideArea3)/face->l3);
	/* & */

	/* anlge based */
	face->theta1 = computeAngle(face->v1->point(),face->v2->point(),face->v3->point(),face->nb1->point());
	face->theta2 = computeAngle(face->v2->point(),face->v3->point(),face->v1->point(),face->nb2->point());
	face->theta3 = computeAngle(face->v3->point(),face->v1->point(),face->v2->point(),face->nb3->point());
	

	face->SO1[0][0] = face->t1.x() * face->t1.x() / (2*face->area*face->l1);
	face->SO1[0][1] = face->t1.x() * face->t1.y() / (2*face->area*face->l1);
	face->SO1[1][0] = face->t1.y() * face->t1.x() / (2*face->area*face->l1);
	face->SO1[1][1] = face->t1.y() * face->t1.y() / (2*face->area*face->l1);

	face->SO2[0][0] = face->t2.x() * face->t2.x() / (2*face->area*face->l2);
	face->SO2[0][1] = face->t2.x() * face->t2.y() / (2*face->area*face->l2);
	face->SO2[1][0] = face->t2.y() * face->t2.x() / (2*face->area*face->l2);
	face->SO2[1][1] = face->t2.y() * face->t2.y() / (2*face->area*face->l2);

	face->SO3[0][0] = face->t3.x() * face->t3.x() / (2*face->area*face->l3);
	face->SO3[0][1] = face->t3.x() * face->t3.y() / (2*face->area*face->l3);
	face->SO3[1][0] = face->t3.y() * face->t3.x() / (2*face->area*face->l3);
	face->SO3[1][1] = face->t3.y() * face->t3.y() / (2*face->area*face->l3);
}

void BasicMesh::translateMesh(Vertex2FT_property_map vertex2k1_pm,Vertex2FT_property_map vertex2k2_pm){
	faceNum = (int)P.size_of_facets();
	vertexNum = (int)P.size_of_vertices();
	edgeNum = (int)P.size_of_halfedges();

	if (faceList != NULL) delete []faceList;
	faceList = new facet[faceNum];

	std::map<Facet_const_handle, int> facetMap;

	PolyhedralSurf::Facet_const_iterator itfb = P.facets_begin();
	PolyhedralSurf::Facet_const_iterator itfe = P.facets_end();

	for (int i = 0;itfb != itfe; itfb++,i++){
		PolyhedralSurf::Facet_const_handle f = itfb;
		facetMap.insert(pair<Facet_const_handle,int>(f,i));

		Vector_3 fnorm = f->getUnitNormal();
		faceList[i].normal = threeTuple(fnorm);
		faceList[i].area = computeArea(f);

		PolyhedralSurf::Halfedge_const_handle h1 = f->halfedge();
		PolyhedralSurf::Vertex_const_handle v1 = h1->vertex();
		faceList[i].nb3 = h1->opposite()->next()->vertex();
		faceList[i].sideArea3 = computeArea(h1->opposite()->facet());

		PolyhedralSurf::Halfedge_const_handle h2 = h1->next();
		PolyhedralSurf::Vertex_const_handle v2 = h2->vertex();
		faceList[i].nb1 = h2->opposite()->next()->vertex();
		faceList[i].sideArea1 = computeArea(h2->opposite()->facet());

		PolyhedralSurf::Halfedge_const_handle h3 = h2->next();
		PolyhedralSurf::Vertex_const_handle v3 = h3->vertex();
		faceList[i].nb2 = h3->opposite()->next()->vertex();
		faceList[i].sideArea2 = computeArea(h3->opposite()->facet());

		faceList[i].v1 = v1;
		faceList[i].v2 = v2;
		faceList[i].v3 = v3;
		faceList[i].index[1] = indexMap.find(v1)->second;
		faceList[i].index[2] = indexMap.find(v2)->second;
		faceList[i].index[3] = indexMap.find(v3)->second;
		faceList[i].nbIdx1	 = indexMap.find(faceList[i].nb1)->second;
		faceList[i].nbIdx2	 = indexMap.find(faceList[i].nb2)->second;
		faceList[i].nbIdx3	 = indexMap.find(faceList[i].nb3)->second;

		faceList[i].p1 = threeTuple(v1->point());
		faceList[i].p2 = threeTuple(v2->point());
		faceList[i].p3 = threeTuple(v3->point());

		//compute stretch inverse matrix
		ComputeStretchInverse(v1,v2,v3,&faceList[i]);

		faceList[i].l1 = sqrt(faceList[i].local_v1.squared_length());
		faceList[i].l2 = sqrt(faceList[i].local_v2.squared_length());
		faceList[i].l3 = sqrt(faceList[i].local_v3.squared_length());

		//compute shape operator
		faceList[i].isBorder = false;
		if ((h1->is_border_edge()) || (h2->is_border_edge()) || (h3->is_border_edge()))
			faceList[i].isBorder = true;
		else
			ComputeShapeOperator(&faceList[i]);
	}

	itfb = P.facets_begin();
	itfe = P.facets_end();
	std::map<Facet_const_handle, int>::iterator iter;

	for (int i = 0;itfb != itfe; itfb++,i++){
		PolyhedralSurf::Facet_const_handle fh = itfb;

		PolyhedralSurf::Halfedge_const_handle h = fh->halfedge();
		if (h->is_border_edge())
			faceList[i].faceIdx3 = -1;
		else{
			Facet_const_handle f = h->opposite()->facet();
			iter = facetMap.find(f);
			faceList[i].faceIdx3 = iter->second;
		}


		h = h->next();
		if (h->is_border_edge())
			faceList[i].faceIdx1 = -1;
		else{
			Facet_const_handle f = h->opposite()->facet();
			iter = facetMap.find(f);
			faceList[i].faceIdx1 = iter->second;
		}

		h = h->next();
		if (h->is_border_edge())
			faceList[i].faceIdx2 = -1;
		else{	
			Facet_const_handle f = h->opposite()->facet();
			iter = facetMap.find(f);
			faceList[i].faceIdx2 = iter->second;
		}
	}

	for (int i = 0; i < faceNum; i++){
		minx = min(minx, min(faceList[i].p1.x,min(faceList[i].p2.x,faceList[i].p3.x)));
		miny = min(miny, min(faceList[i].p1.y,min(faceList[i].p2.y,faceList[i].p3.y)));
		minz = min(minz, min(faceList[i].p1.z,min(faceList[i].p2.z,faceList[i].p3.z)));

		maxx = max(maxx, max(faceList[i].p1.x,max(faceList[i].p2.x,faceList[i].p3.x)));
		maxy = max(maxy, max(faceList[i].p1.y,max(faceList[i].p2.y,faceList[i].p3.y)));
		maxz = max(maxz, max(faceList[i].p1.z,max(faceList[i].p2.z,faceList[i].p3.z)));
	}
	max_all = max(max(maxx-minx,maxy-miny),maxz-minz);
}

void BasicMesh::constructVertexList(){
	if (vertexList != NULL) delete []vertexList;
	vertexList = new vertex[P.size_of_vertices()];

	PolyhedralSurf::Vertex_const_iterator vb = P.vertices_begin();
	PolyhedralSurf::Vertex_const_iterator ve = P.vertices_end();
	PolyhedralSurf::Vertex_const_handle vh,nb;

	for (int i = 0; vb != ve; vb++,i++){
		vh = vb;

		PolyhedralSurf::Halfedge_around_vertex_const_circulator tempHfe = vh->vertex_begin();

		vertexList[i].idx = indexMap.find(vh)->second;
		vertexList[i].neighbourNum = vh->degree();
		vertexList[i].nbIdx = new int[vh->degree()];
		vertexList[i].weight = new double[vh->degree()];
		vertexList[i].area = 0;

		for (int j = 0; j<vh->degree(); j++){
			nb = tempHfe->opposite()->vertex();
			vertexList[i].nbIdx[j] = indexMap.find(nb)->second;
			vertexList[i].area += computeArea(tempHfe->facet());
			tempHfe++;
		}

		for (int j = 0; j<vh->degree(); j++){
			int previous = (j-1+vh->degree())%vh->degree();
			int next = (j+1)%vh->degree();

			Vector_3 v1(vertexIndex[vertexList[i].nbIdx[previous]]->point(),vh->point());
			Vector_3 v2(vertexIndex[vertexList[i].nbIdx[previous]]->point(),
				vertexIndex[vertexList[i].nbIdx[j]]->point());

			double angle1 = computeAngle(v1,v2);

			Vector_3 v3(vertexIndex[vertexList[i].nbIdx[next]]->point(),vh->point());
			Vector_3 v4(vertexIndex[vertexList[i].nbIdx[next]]->point(),
				vertexIndex[vertexList[i].nbIdx[j]]->point());

			double angle2 = computeAngle(v3,v4);

			vertexList[i].weight[j] = cos(angle1)/sin(angle1)+cos(angle2)/sin(angle2);
			vertexList[i].weight[j] = max_zenyo(0.1,vertexList[i].weight[j]);
		}
	}
}

void BasicMesh::constructEdgeList(){
	if (edgeList != NULL) delete []edgeList;
	edgeList = new edge[edgeNum/2];
	std::map<PolyhedralSurf::Halfedge_const_handle, int> edgeMap;
	edgeMap.clear();

	PolyhedralSurf::Halfedge_const_iterator iteb = P.halfedges_begin();
	PolyhedralSurf::Halfedge_const_iterator itee = P.halfedges_end();
	std::map<PolyhedralSurf::Halfedge_const_handle, int>::iterator iterEdge;
	std::map<PolyhedralSurf::Vertex_const_handle, int>::iterator iterV;

	int idx = 0;

	for (int i = 0;iteb != itee; iteb++,i++){
		PolyhedralSurf::Halfedge_const_handle eh = iteb;

		iterEdge = edgeMap.find(eh->opposite());
		if (iterEdge == edgeMap.end()){
			edgeList[idx].v1 = eh->opposite()->vertex();
			edgeList[idx].v2 = eh->vertex();
			
			iterV = indexMap.find(edgeList[idx].v1);
			edgeList[idx].index1 = iterV->second;

			iterV = indexMap.find(edgeList[idx].v2);
			edgeList[idx].index2 = iterV->second;

			double mainArea,sideArea;

			if (eh->is_border()) mainArea = 0; else mainArea = computeArea(eh->facet());
			if (eh->opposite()->is_border()) sideArea = 0; else sideArea = computeArea(eh->opposite()->facet());
			
			double area = mainArea + sideArea;
			double l = computeLength(eh);

			edgeList[idx].stretchWeight = area / l / l;

			edgeMap.insert(pair<PolyhedralSurf::Halfedge_const_handle,int>(eh,idx));
			edgeList[idx].eh = eh;

			// boundary doesn't yield bending energy
			edgeList[idx].isBoundary = false;
			if (eh->is_border_edge()) {
				edgeList[idx].isBoundary = true;
				idx++;
				continue;
			}

			edgeList[idx].vl = eh->opposite()->next()->vertex();
			edgeList[idx].vr = eh->next()->vertex();

			iterV = indexMap.find(edgeList[idx].vl);
			edgeList[idx].indexl = iterV->second;

			iterV = indexMap.find(edgeList[idx].vr);
			edgeList[idx].indexr = iterV->second;

			edgeList[idx].index[1] = edgeList[idx].index1;
			edgeList[idx].index[2] = edgeList[idx].index2;
			edgeList[idx].index[3] = edgeList[idx].indexl;
			edgeList[idx].index[4] = edgeList[idx].indexr;

			//compute angle
			/*
			Vector_3 v1(edgeList[idx].v1->point(),edgeList[idx].vr->point());
			Vector_3 v2(edgeList[idx].v1->point(),edgeList[idx].v2->point());
			Vector_3 v3(edgeList[idx].v1->point(),edgeList[idx].vl->point());

			Vector_3 n1 = cross_product(v1,v2);
			Vector_3 n2 = cross_product(v2,v3);

			double sm1 = n1.squared_length();
			double sq1 = sqrt(sm1);
			double sm2 = n2.squared_length();
			double sq2 = sqrt(sm2);

			double p = n1*n2;
			double g = sq1*sq2;

			double k = p/g;

			if (k > 1) k = 1;
			if (k < -1) k = -1;

			edgeList[idx].angle = acos(k);

			if (_isnan(edgeList[idx].angle)) cout<<"~~~prepration"<<endl;
			*/

			Point_3 v1 = edgeList[idx].v1->point();
			Point_3 v2 = edgeList[idx].v2->point();
			Point_3 vr = edgeList[idx].vr->point();
			Point_3 vl = edgeList[idx].vl->point();
			
			double yz1,xz1,xy1,xyz;
			yz1 = v2.y()*vr.z()+v2.z()*vl.y()+vr.y()*vl.z()-vr.z()*vl.y()-v2.z()*vr.y()-v2.y()*vl.z();
			xz1 = v2.x()*vr.z()+v2.z()*vl.x()+vr.x()*vl.z()-vr.z()*vl.x()-v2.z()*vr.x()-v2.x()*vl.z();
			xy1 = v2.x()*vr.y()+v2.y()*vl.x()+vr.x()*vl.y()-vr.y()*vl.x()-v2.y()*vr.x()-v2.x()*vl.y();
			xyz = v2.x()*vr.y()*vl.z()+v2.y()*vr.z()*vl.x()+vr.x()*vl.y()*v2.z()-v2.z()*vr.y()*vl.x()-v2.y()*vr.x()*vl.z()-v2.x()*vr.z()*vl.y();

			double det = v1.x()*yz1-v1.y()*xz1+v1.z()*xy1-xyz;
			edgeList[idx].angle = det;
			
			//edgeList[idx].bendWeight = l * l / area;
			double angleDominator = (mainArea*sideArea)/l;
			edgeList[idx].bendWeight =  l * l / area / angleDominator / angleDominator;

			idx++;
		}
	}
}

void BasicMesh::initialDeformation(){
	if (u != NULL) delete[] u;
	u = lbfgs_malloc(vertexNum * 3);
	memset(u,0,sizeof(lbfgsfloatval_t)*vertexNum*3);

	if (targetPos != NULL) delete[] targetPos;
	targetPos = new Vector_3[vertexNum];

	if (overallWeight != NULL) delete[] overallWeight;
	overallWeight = new double[vertexNum];
	memset(overallWeight,0,sizeof(double)*vertexNum);


	if (bestMatch != NULL) {
		for (int i = 0; i < surfaceNum; i++)
			if (bestMatch[i] != NULL) delete[] bestMatch[i];

		delete[] bestMatch;
	}
	if (matchWeight != NULL) delete[] matchWeight;
	bestMatch = new Vector_3*[surfaceNum];
	matchWeight = new float*[surfaceNum];

	for (int i = 0; i < surfaceNum; i++){
		bestMatch[i] = new Vector_3[vertexNum];
		matchWeight[i] = new float[vertexNum];

		for (int j = 0; j < vertexNum; j++)
			bestMatch[i][j] = Vector_3(0,0,0);

		memset(matchWeight[i],0,sizeof(float)*vertexNum);
	}
}

void BasicMesh::summerizeForce(){
	int startSurface = max_zenyo(idx - GEOMETRIC_RANGE, 0);
	int endSurface = min_zenyo(idx + GEOMETRIC_RANGE, surfaceNum - 1);

	for (int j = 0; j < vertexNum; j++){
		Vector_3 force = Vector_3(0,0,0);
		double weightSum = 0;

		for (int i = startSurface; i <= endSurface; i++)
			if (i != idx){
				force = force + bestMatch[i][j] * matchWeight[i][j];

				weightSum += matchWeight[i][j];
			}
		
			if (weightSum > 0.001){
				targetPos[j] = force/weightSum;
				overallWeight[j] = weightSum;
			}
			else
				overallWeight[j] = 0;
		
	}
}


void BasicMesh::constructLandmark(string filename){
	ifstream fin;
	fin.open(filename.c_str());
	double temp;

	landmarkNum = new int[surfaceNum];
	landmarks1 = new int*[surfaceNum];
	landmarks2 = new int*[surfaceNum];

	int startIdx = max_zenyo(0,idx - TEXTURE_RANGE);
	int endIdx = min_zenyo(surfaceNum-1,idx + TEXTURE_RANGE);
	for (int i = startIdx; i <= endIdx; i++){
		if (i == idx) continue;
		fin>>temp;
		fin>>landmarkNum[i];

		landmarks1[i] = new int[landmarkNum[i]];
		landmarks2[i] = new int[landmarkNum[i]];

		for (int j = 0; j < landmarkNum[i]; j++){
			fin >> temp;
			landmarks1[i][j] = int(temp);
		}
		for (int j = 0; j < landmarkNum[i]; j++){
			fin >> temp;
			landmarks2[i][j] = int(temp);
		}

	}

	fin.close();
}
