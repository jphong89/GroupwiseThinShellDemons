#include "stdafx.h"
#include "GlutManager.h"
#include <iostream>
using namespace std;

// 
// int beginGlut(int argc, char* argv[]){
// 	glutInit(&argc,argv);
// 	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH ); 
// 	glutInitWindowSize(400, 400); 
// 
// 	glutCreateWindow("OpenGL"); 
// 	glutReshapeFunc(reshape);
// 	glutDisplayFunc(display);
// 	init();
// 	glutKeyboardFunc(ProcessKeyboard);
// 	glutSpecialFunc(ProcessSpecialKeyboead);
// 	glutMouseFunc(Mouse);
// 
// 	if (pthread_mutex_init(&lock, NULL) != 0)
// 	{
// 		cout<<"mutex init failed\n"<<endl;
// 		return 1;
// 	}
// 
// 	cout<<"Begin Rendering..."<<endl;
// 	glutMainLoop();
// 	return 0;
// }
// 
// void ProcessKeyboard(unsigned char key,int x,int y)
// {
// 	if (key == '=')
// 		camera_Scale = camera_Scale + 0.1;
// 	else if (key == '-')
// 		camera_Scale = camera_Scale - 0.1;
// 	else if (key == 'w')
// 		rotate_X = rotate_X - 2;
// 	else if (key == 's')
// 		rotate_X = rotate_X + 2;
// 	else if (key == 'a')
// 		rotate_Y = rotate_Y - 2;
// 	else if (key == 'd')
// 		rotate_Y = rotate_Y + 2;
// 	else if (key == 'q')
// 		rotate_Z = rotate_Z - 2;
// 	else if (key == 'e')
// 		rotate_Z = rotate_Z + 2;
// 	else if (key == 'b'){
// 		if (GPU_ON)
// 			pthread_t tid;
// 		else{
// 			pthread_t tid;
// 			int err = pthread_create(&tid, NULL, &startOptimization, NULL);
// 			if (err != 0){
// 				cout<<"cannot creat thread"<<endl;
// 				return;
// 			}
// 
// 		}
// 	}
// 	else if (key == '1')
// 		switch1 = !switch1;
// 	else if (key == '2')
// 		switch2 = !switch2;
// 	else if (key == '3')
// 		switch3 = !switch3;
// 	else if (key == '4')
// 		switch4 = !switch4;
// 
// 	display();
// }
// 
// void ProcessSpecialKeyboead(int key, int x, int y) 
// {
// 	if (key == GLUT_KEY_LEFT)
// 		camera_Right = camera_Right - 0.1;
// 	else if (key == GLUT_KEY_RIGHT)
// 		camera_Right = camera_Right + 0.1;
// 	else if (key == GLUT_KEY_UP)
// 		camera_Up = camera_Up + 0.1;
// 	else if (key == GLUT_KEY_DOWN)
// 		camera_Up = camera_Up - 0.1;
// 
// 	display();
// }
// 
// void Mouse(int button, int state, int x, int y) //处理鼠标点击
// {
// 	GLint viewport[4];   
// 	GLdouble modelview[16];   
// 	GLdouble projection[16];   
// 	GLfloat winX, winY, winZ;   
// 	GLdouble posX, posY, posZ;   
// 
// 	glPushMatrix(); 
// 	glLoadIdentity();    
//   
// 	glTranslatef(camera_Right,camera_Up,-6.0f);
// 
// 	glRotatef(rotate_X,1.0f,0.0f,0.0f);	
// 	glRotatef(rotate_Y,0.0f,1.0f,0.0f);	
// 	glRotatef(rotate_Z,0.0f,0.0f,1.0f);	
// 
// 	glGetIntegerv(GL_VIEWPORT, viewport);
// 	glGetDoublev(GL_MODELVIEW_MATRIX, modelview);   
// 	glGetDoublev(GL_PROJECTION_MATRIX, projection);   
// 	glPopMatrix();   
// 
// 	winX = x;   
// 	winY = viewport[3] - (float)y;   
// 	glReadPixels((int)winX, (int)winY, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ);   
// 	gluUnProject(winX, winY, winZ, modelview, projection, viewport, &posX, &posY, &posZ); 
// 
// 	threeTuple a;
// 	a.x = posX / camera_Scale * max_all + (maxx - minx) / 2 + minx;
// 	a.y = posY / camera_Scale * max_all + (maxy - miny) / 2 + miny;
// 	a.z = posZ / camera_Scale * max_all + (maxz - minz) / 2 + minz;
// 	
// 	cout<<"seleted coordinates:    "<<a.x<<' '<<a.y<<' '<<a.z<<endl;
// 
// 	PolyhedralSurf::Vertex_const_handle vh;
// 
// 	if (switch1) vh = mouseEvent(RCSurface,a);
// 	if (switch2) vh = mouseEvent(CTSurface,a);
// 	
// 	//draw CT affinity
// 	if (switch1){
// 		CTSurface->meshColor.clear();
// 
// 		int center = RCSurface->indexMap.find(vh)->second;
// 		cout<<matchWeight[center]<<endl;
// 		
// 		PolyhedralSurf::Vertex_const_iterator vb = CTSurface->P.vertices_begin();
// 		PolyhedralSurf::Vertex_const_iterator ve = CTSurface->P.vertices_end();
// 
// 		double maxAffinity = 0;
// 		int maxI = 0;
// 		for (int i = 0; vb != ve; vb++,i++){
// 			vh = vb;
// 			CTSurface->meshColor.insert(pair<Vertex_const_handle,float>(vh,VAL(affinity,center,i,affinityN)));
// 			if (VAL(affinity,center,i,affinityN) > maxAffinity){
// 				maxAffinity = VAL(affinity,center,i,affinityN);
// 				maxI = i;
// 			}
// 		}
// 
// 		CTSurface->centre = CTSurface->vertexIndex[maxI];
// 
// 		//cout<<center<<endl;
// 		CTSurface->deleteMarchList();
// 		CTSurface->vertexBool.clear();
// 		for (int j = 0; j < CTSurface->dirNum; j++)
// 			CTSurface->findSignature(j,CTSurface->centre,0);
// 
// 		for (int i = 0; i < CTSurface->faceNum; i++){
// 			float color = CTSurface->meshColor.find(CTSurface->faceList[i].v1)->second;
// 			CTSurface->faceList[i].color1 = getJetColor(color);
// 
// 			color = CTSurface->meshColor.find(CTSurface->faceList[i].v2)->second;
// 			CTSurface->faceList[i].color2 = getJetColor(color);
// 
// 			color = CTSurface->meshColor.find(CTSurface->faceList[i].v3)->second;
// 			CTSurface->faceList[i].color3 = getJetColor(color);
// 		}
// 	}
//  	
// 	display();
// }
// 
// PolyhedralSurf::Vertex_const_handle mouseEvent(BasicMesh* mesh,threeTuple a){
// 	PolyhedralSurf::Vertex_const_handle vh = mesh->findVertexHandle(a);
// 
// 	int center = mesh->indexMap.find(vh)->second;
// 
// 	mesh->deleteMarchList();
// 	mesh->vertexBool.clear();
// 	for (int j = 0; j < mesh->dirNum; j++)
// 		mesh->findSignature(j,mesh->centre,0);
// 
// 	return vh;
// 
// }
// 
// void init(){
// 	glEnable(GL_DEPTH_TEST);			//开启深度测试
// 	
// 	glEnable(GL_LIGHT0);								// Enable Light One
// 	glEnable(GL_LIGHTING); 
// 
// 	GLfloat light_ambient [] = { 0.3,0.3,0.3, 1.0 };
// 	GLfloat light_diffuse [] = { 0.3, 0.3, 0.3, 1.0 };
// 	GLfloat light_specular [] = { 0.2, 0.2, 0.2, 1.0 };
// 	GLfloat light_position[] = {0.0f, 0.0f, 0.0f, 1.0 };
// 
// 	glLightfv(GL_LIGHT0, GL_AMBIENT,   light_ambient);
// 	glLightfv(GL_LIGHT0, GL_DIFFUSE,   light_diffuse);
// 	glLightfv(GL_LIGHT0, GL_SPECULAR,  light_specular);
// 	glLightfv(GL_LIGHT0, GL_POSITION,  light_position);
// 	glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, 1.0f);
// 
// 	glEnable(GL_COLOR_MATERIAL);
// 
// 	GLfloat ambient [] = { 0.029412, 0.223529, 0.327451, 1.0 };
// 	GLfloat diffuse [] = { 0.180392, 0.268627, 0.813725, 1.0 };
// 	GLfloat specular[] = { 0.492157, 0.441176, 0.807843, 1.0 };
// 	GLfloat shininess[] = { 50.8974f };
// 	GLfloat emission[]  = {0.3f,0.3f,0.3f,1.0f};
// 
// 	glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT,   ambient);
// 	glMaterialfv( GL_FRONT_AND_BACK, GL_DIFFUSE,   diffuse);
// 	glMaterialfv( GL_FRONT_AND_BACK, GL_SPECULAR,  specular);
// 	glMaterialfv( GL_FRONT_AND_BACK, GL_SHININESS, shininess);
// 	glMaterialfv( GL_FRONT_AND_BACK, GL_EMISSION,  emission);
// 
// 	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
// 
// 	cout<<"Initializing..."<<endl;
// 	
// 	paintColor(RCSurface,threeTuple(0.8,0.8,0.8),threeTuple(0.8,0.8,0.8),threeTuple(0.8,0.8,0.8));
// }
// 
// void display(){
// 	pthread_mutex_lock(&lock);
// 	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// Clear The Screen And The Depth Buffer
// 	glLoadIdentity();
// 	//m_Camera.DisplayScene();
// 	glTranslatef(camera_Right,camera_Up,-6.0f);
// 
// 
// 	glRotatef(rotate_X,1.0f,0.0f,0.0f);	
// 	glRotatef(rotate_Y,0.0f,1.0f,0.0f);	
// 	glRotatef(rotate_Z,0.0f,0.0f,1.0f);	
// 
// 	if (switch1) {
// 		drawSurface(RCSurface,switch3);
// 		drawExtra(RCSurface);
// 	}
// 	if (switch2) {
// 		drawSurface(CTSurface,0);
// 		drawExtra(CTSurface);
// 	}
// 
// 	glFlush();
// 	glutSwapBuffers();
// 	pthread_mutex_unlock(&lock);
// }
// 
// void reshape(int w,int h)
// {
// 	glViewport( 0, 0, w, h );
// 	glMatrixMode( GL_PROJECTION );
// 	glLoadIdentity();
// 	gluPerspective(75.0f, (float)w/h, 1.0f, 1000.0f);
// 	glMatrixMode( GL_MODELVIEW );
// }
// 
// void drawSurface(BasicMesh* mesh,int deform){
// 	threeTuple color;
// 
// 	Vector_3 pd1,pd2;
// 
// 	float colorRatio = 1;
// 
// 	glLineWidth(5.0);
// 	glColor3f(1.0f,0.0f,0.0f);
// 
// 	for (int i = 0; i < mesh->linkNum; i++){
// 		Point_3 p1 = mesh->vertexIndex[linkList[i].index1]->point();
// 		Point_3 p2 = mesh->vertexIndex[linkList[i].index2]->point();
// 		int index1 = linkList[i].index1;
// 		int index2 = linkList[i].index2;
// 
// 		threeTuple a;
// 		a.x = (p1.x() + deform*u[index1*3] -  (maxx - minx) / 2 - minx) / max_all * camera_Scale;
// 		a.y = (p1.y() + deform*u[index1*3+1]  - (maxy - miny) / 2 - miny) / max_all * camera_Scale;
// 		a.z = (p1.z() + deform*u[index1*3+2]  - (maxz - minz) / 2 - minz) / max_all * camera_Scale;
// 
// 		threeTuple b;
// 		b.x = (p2.x() + deform*u[index2*3] - (maxx - minx) / 2 - minx) / max_all * camera_Scale;
// 		b.y = (p2.y() + deform*u[index2*3+1] - (maxy - miny) / 2 - miny) / max_all * camera_Scale;
// 		b.z = (p2.z() + deform*u[index2*3+2] - (maxz - minz) / 2 - minz) / max_all * camera_Scale;
// 
// 		glBegin(GL_LINES);
// 			glVertex3f(a.x,a.y,a.z);
// 			glVertex3f(b.x,b.y,b.z);
// 		glEnd();
// 	}
// 
// 	for (int i = 0; i < mesh->faceNum; i++){
// 		facet f = mesh->faceList[i];
// 
// 		f.p1.x += deform*u[f.index[1]*3];
// 		f.p1.y += deform*u[f.index[1]*3+1];
// 		f.p1.z += deform*u[f.index[1]*3+2];
// 
// 		f.p2.x += deform*u[f.index[2]*3];
// 		f.p2.y += deform*u[f.index[2]*3+1];
// 		f.p2.z += deform*u[f.index[2]*3+2];
// 
// 		f.p3.x += deform*u[f.index[3]*3];
// 		f.p3.y += deform*u[f.index[3]*3+1];
// 		f.p3.z += deform*u[f.index[3]*3+2];
// 		
// 		threeTuple a;
// 		a.x = (f.p1.x-  (maxx - minx) / 2 - minx) / max_all * camera_Scale;
// 		a.y = (f.p1.y - (maxy - miny) / 2 - miny) / max_all * camera_Scale;
// 		a.z = (f.p1.z - (maxz - minz) / 2 - minz) / max_all * camera_Scale;
// 
// 		threeTuple b;
// 		b.x = (f.p2.x - (maxx - minx) / 2 - minx) / max_all * camera_Scale;
// 		b.y = (f.p2.y - (maxy - miny) / 2 - miny) / max_all * camera_Scale;
// 		b.z = (f.p2.z - (maxz - minz) / 2 - minz) / max_all * camera_Scale;
// 
// 		threeTuple c;
// 		c.x = (f.p3.x - (maxx - minx) / 2 - minx) / max_all * camera_Scale;
// 		c.y = (f.p3.y - (maxy - miny) / 2 - miny) / max_all * camera_Scale;
// 		c.z = (f.p3.z - (maxz - minz) / 2 - minz) / max_all * camera_Scale;
// 
// 		glLineWidth(1.0);
// 		glBegin(GL_TRIANGLES);
// 		glNormal3f(f.normal.x,f.normal.y,f.normal.z);
// 		
// 		if (occluded[f.index[1]]) glColor3f(1.0f,0.0f,0.0f);
// 		else glColor3f(f.color1.x,f.color1.y,f.color1.z);
// 
// 		glVertex3f(a.x,a.y,a.z);
// 		
// 		
// 		if (occluded[f.index[2]]) glColor3f(1.0f,0.0f,0.0f);
// 		else glColor3f(f.color2.x,f.color2.y,f.color2.z);
// 
// 		glVertex3f(b.x,b.y,b.z);
// 
// 
// 		if (occluded[f.index[3]]) glColor3f(1.0f,0.0f,0.0f);
// 		else glColor3f(f.color3.x,f.color3.y,f.color3.z);
// 
// 		glVertex3f(c.x,c.y,c.z);
// 		glEnd();
// 
// 		glLineWidth(3.0);
// 
// 	}
// 
// }
// 
// void drawExtra(BasicMesh* mesh){
// 	glLineWidth(5.0);
// 
// 	threeTuple color;
// 	for (int i = 0; i < mesh->dirNum; i++){
// 
// 		color = getJetColor(0.125*(i+1)); 
// 		glColor3f(color.x,color.y,color.z);
// 
// 		threeTuple* marchList = mesh->marchList[i];
// 
// 		glBegin(GL_LINE_STRIP);
// 		while (marchList != NULL){
// 
// 			threeTuple a;
// 			a.x = (marchList->x - (maxx - minx) / 2 - minx) / max_all * camera_Scale;
// 			a.y = (marchList->y - (maxy - miny) / 2 - miny) / max_all * camera_Scale;
// 			a.z = (marchList->z - (maxz - minz) / 2 - minz) / max_all * camera_Scale;
// 			glVertex3f(a.x,a.y,a.z);
// 
// 			marchList = marchList->next;
// 		}
// 		glEnd();
// 	}
// 
// 	glColor3f(1,0,0);
// 	
// 	for (int i = 0; i < mesh->landmarkNum; i++){
// 		Point_3 p = mesh->vertexIndex[mesh->landmark[i]]->point();
// 
// 		threeTuple a;
// 		a.x = (p.x() - (maxx - minx) / 2 - minx) / max_all * camera_Scale;
// 		a.y = (p.y() - (maxy - miny) / 2 - miny) / max_all * camera_Scale;
// 		a.z = (p.z() - (maxz - minz) / 2 - minz) / max_all * camera_Scale;
// 		drawSphere(a.x,a.y,a.z,0.5/max_all*camera_Scale,10,10);
// 	}
// }
// 
// void drawSphere(GLfloat xx, GLfloat yy, GLfloat zz, GLfloat radius, GLfloat M, GLfloat N)
// {
// 	float step_z = PI/M;
// 	float step_xy = 2*PI/N;
// 	float x[4],y[4],z[4];
// 	float angle_z = 0.0;
// 	float angle_xy = 0.0;
// 	int i=0, j=0;
// 	glBegin(GL_QUADS);
// 	for(i=0; i<M; i++)
// 	{
// 		angle_z = i * step_z;
// 
// 		for(j=0; j<N; j++)
// 		{
// 			angle_xy = j * step_xy;
// 			x[0] = radius * sin(angle_z) * cos(angle_xy);
// 			y[0] = radius * sin(angle_z) * sin(angle_xy);
// 			z[0] = radius * cos(angle_z);
// 			x[1] = radius * sin(angle_z + step_z) * cos(angle_xy);
// 			y[1] = radius * sin(angle_z + step_z) * sin(angle_xy);
// 			z[1] = radius * cos(angle_z + step_z);
// 			x[2] = radius * sin(angle_z + step_z) * cos(angle_xy + step_xy);
// 			y[2] = radius * sin(angle_z + step_z) * sin(angle_xy + step_xy);
// 			z[2] = radius * cos(angle_z + step_z);
// 			x[3] = radius * sin(angle_z) * cos(angle_xy + step_xy);
// 			y[3] = radius * sin(angle_z) * sin(angle_xy + step_xy);
// 			z[3] = radius * cos(angle_z);
// 			for(int k=0; k<4; k++)
// 			{
// 				glVertex3f(xx+x[k], yy+y[k],zz+z[k]);
// 			}
// 		}
// 	}
// 	glEnd();
// }
// 
// void paintColor(BasicMesh* mesh, threeTuple color1,threeTuple color2,threeTuple color3){
// 	for (int i = 0; i < mesh->faceNum; i++){
// 		mesh->faceList[i].color1 = color1;
// 		mesh->faceList[i].color2 = color2;
// 		mesh->faceList[i].color3 = color3;
// 	}
// }
// 
// 
