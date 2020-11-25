

#include "Render.h"
#include <math.h>
#include <Windows.h>
#include <GL\GL.h>
#include <GL\GLU.h>
#include "MyVEC3D.h"
#include <vector>
#include <array>

std::vector<std::vector<Vector3>> points(
	{
		{
			{3, 0, 3},
			{3, -1, 2},
			{3, -2, 2},
			{3, -3, 1}
		},
		{
			{4, 0, 1},
			{4, -1, 2},
			{4, -2, -2},
			{4, -3, 1}
		},
		{
			{5, 0, 1},
			{5, -1, 2},
			{5, -2, 2},
			{5, -3, 1}
		},
		{
			{6, 0, 1},
			{6, -1, 2},
			{6, -2, 2},
			{6, -3, 1}
		}
	});

unsigned long int factorial(int i)
{
	if (i == 0) return 1;
	else return i * factorial(i - 1);
}

double Bernstein(double u, double n, int index) {
	return (factorial(n) / (factorial(index) * factorial(n - index))) * pow(u, index) * pow(1 - u, n - index);
}

void BezeSurfacePoint(double u, double v, Vector3& vec) {
	Vector3 new_v;
	int n = 3, m = 3;
	for (size_t i = 0; i < points.size(); ++i) {
		for (size_t j = 0; j < points[i].size(); ++j)
		{
			new_v += points[i][j] * Bernstein(u, n, i) * Bernstein(v, m, j);
		}
	}
	vec = new_v;
}

std::array<Vector3, 4> Triangle(double u, double v, double h) {
	std::array<Vector3, 4> tmp;
	BezeSurfacePoint(u, v, tmp[0]);
	BezeSurfacePoint(u, v + h, tmp[1]);
	BezeSurfacePoint(u + h, v, tmp[2]);
	BezeSurfacePoint(u + h, v + h, tmp[3]);

	return tmp;
}
inline double Bez(double p0, double p1, double p2, double p3, double t)
{
	return p0 * (1 - t) * (1 - t) * (1 - t) + 3 * p1 * t * (1 - t) * (1 - t) + 3 * t * t * p2 * (1 - t) + p3 * t * t * t; //посчитаная формула
}
inline double Ermit(double p1, double p4, double r1, double r4, double t)
{
	return p1 * (2 * t * t * t - 3 * t * t + 1) + p4 * (-2 * t * t * t + 3 * t * t) + (r1 - p1) * (t * t * t - 2 * t * t + t) + (r4 - p4) * (t * t * t -  t * t); //посчитаная формула
}
double f(double p1, double p2, double p3, double t)
{
	return p1 * (1 - t) * (1 - t) + 2 * p2 * t * (1 - t) + p3 * t * t; //посчитанная формула
}

double t_max = 0;
bool flagReverse = false;

Vector3 Bize(double p0[], double p1[], double p2[], double p3[], double t)
{
	Vector3 Vec;
	Vec.setCoords(Bez(p0[0], p1[0], p2[0], p3[0], t), Bez(p0[1], p1[1], p2[1], p3[1], t), Bez(p0[2], p1[2], p2[2], p3[2], t));
	return Vec;
}



void BizePlane()
{
	double h = 0.1;
	for (double u = 0; u < 1.01 - h; u += h) {
		glBegin(GL_TRIANGLES);
		for (double v = 0; v < 1.01 - h; v += h) {
			std::array<Vector3, 4> arr = Triangle(u, v, h);

			/*glNormal3dv(FindNormal(arr[0].toArray(), arr[1].toArray(), arr[2].toArray(), 1));
			glTexCoord3dv(TextureCoords(arr[0]));*/
			glVertex3dv(arr[0].toArray());
			//glTexCoord3dv(TextureCoords(arr[1]));
			glVertex3dv(arr[1].toArray());
			//glTexCoord3dv(TextureCoords(arr[3]));
			glVertex3dv(arr[3].toArray());

			//glNormal3dv(FindNormal(arr[0].toArray(), arr[1].toArray(), arr[2].toArray(), 1));
			//glTexCoord3dv(TextureCoords(arr[0]));
			glVertex3dv(arr[0].toArray());
			//glTexCoord3dv(TextureCoords(arr[2]));
			glVertex3dv(arr[2].toArray());
			//glTexCoord3dv(TextureCoords(arr[3]));
			glVertex3dv(arr[3].toArray());
		}
		glEnd();
	}
}

void Render(double delta_time)
{
	glColor3d(0.6, 0.5, 0);
	BizePlane();
	glColor3d(0, 0, 0);
	double P0[] = { 0,0,0 };
	double P1[] = { 1,-2,3 }; 
	double P2[] = { -4,6,7 };
	double P3[] = { 10,10,0 };

	double P00[] = { -4,0,3 };
	double P11[] = { 1,-3,4 }; 
	double P22[] = { -10,3,5 };
	double P33[] = { -5,5,0 };

	double Erm0[] = { 5,-12, 0 };
	double Erm1[] = { 2,-8,-10 };
	double Erm2[] = { -8,-10,-12 };
	double Erm3[] = { 6,-17,-5 };

	double Erm00[] = { -14,-10,3 };
	double Erm11[] = { -8,-12,0 };
	double Erm22[] = { -6,-8,5 };
	double Erm33[] = { -5,-14,0 };
	glLineWidth(3); //ширина линии
	glBegin(GL_LINE_STRIP);
	for (double t = 0; t <= 1.0001; t += 0.01)
	{
		double P[3];
		P[0] = Bez(P0[0], P1[0], P2[0], P3[0], t);
		P[1] = Bez(P0[1], P1[1], P2[1], P3[1], t);
		P[2] = Bez(P0[2], P1[2], P2[2], P3[2], t);
		glVertex3dv(P); //Рисуем точку P
	}
	glEnd();
	glLineWidth(1);
	glBegin(GL_LINE_STRIP);
	glVertex3d(P0[0], P0[1], P0[2]);
	glVertex3d(P1[0], P1[1], P1[2]);
	glVertex3d(P2[0], P2[1], P2[2]);
	glVertex3d(P3[0], P3[1], P3[2]);
	glEnd();
	
	if (!flagReverse)
	{
		t_max += delta_time / 5; //t_max становится = 1 за 5 секунд
		if (t_max > 1)
		{
			t_max = 1; //после обнуляется
			flagReverse = !flagReverse;
		}
	}
	else
	{
		t_max -= delta_time / 5; //t_max становится = 1 за 5 секунд
		if (t_max < 0)
		{
			t_max = 0; //после обнуляется
			flagReverse = !flagReverse;
		}
	}
	double coord[] = { 0,0,0 };
	
	coord[0] = Bez(P0[0], P1[0], P2[0], P3[0], t_max);
	coord[1] = Bez(P0[1], P1[1], P2[1], P3[1], t_max);
	coord[2] = Bez(P0[2], P1[2], P2[2], P3[2], t_max);
	Vector3 P_old = Bize(P0, P1, P2, P3, !flagReverse ? t_max - delta_time : t_max + delta_time);
	Vector3 P = Bize(P0, P1, P2, P3, t_max);
	Vector3 VecP_P_old = (P - P_old).normolize();

	Vector3 rotateX(VecP_P_old.X(), VecP_P_old.Y(), 0);
	rotateX = rotateX.normolize();

	Vector3 VecPrX = Vector3(1, 0, 0).vectProisvedenie(rotateX);
	double CosX = Vector3(1, 0, 0).ScalarProizv(rotateX);
	double SinAngleZ = VecPrX.Z() / abs(VecPrX.Z());
	double AngleOZ = acos(CosX) * 180 / 3.14 * SinAngleZ;

	double AngleOY = acos(VecP_P_old.Z()) * 180 / 3.14 - 90;


	glColor3d(0, 0, 0);
	glPushMatrix();
	glTranslated(coord[0], coord[1], coord[2]);
	glRotated(AngleOZ, 0, 0, 1);
	glRotated(AngleOY, 0, 1, 0);
	//Cube
	// White side - BACK
	glBegin(GL_POLYGON);
	glColor3f(0, 0.33, 0.555);
	glVertex3f(0.5, -0.5, 0.5);
	glVertex3f(0.5, 0.5, 0.5);
	glVertex3f(-0.5, 0.5, 0.5);
	glVertex3f(-0.5, -0.5, 0.5);
	glEnd();

	// Purple side - RIGHT
	glBegin(GL_POLYGON);
	glVertex3f(0.5, -0.5, -0.5);
	glVertex3f(0.5, 0.5, -0.5);
	glVertex3f(0.5, 0.5, 0.5);
	glVertex3f(0.5, -0.5, 0.5);
	glEnd();

	// Green side - LEFT
	glBegin(GL_POLYGON);
	glVertex3f(-0.5, -0.5, 0.5);
	glVertex3f(-0.5, 0.5, 0.5);
	glVertex3f(-0.5, 0.5, -0.5);
	glVertex3f(-0.5, -0.5, -0.5);
	glEnd();

	// Blue side - TOP
	glBegin(GL_POLYGON);
	glVertex3f(0.5, 0.5, 0.5);
	glVertex3f(0.5, 0.5, -0.5);
	glVertex3f(-0.5, 0.5, -0.5);
	glVertex3f(-0.5, 0.5, 0.5);
	glEnd();

	// Red side - BOTTOM
	glBegin(GL_POLYGON);
	glVertex3f(0.5, -0.5, -0.5);
	glVertex3f(0.5, -0.5, 0.5);
	glVertex3f(-0.5, -0.5, 0.5);
	glVertex3f(-0.5, -0.5, -0.5);
	glEnd();
	glBegin(GL_POLYGON);

	glVertex3f(0.5, -0.5, -0.5);     
	glVertex3f(0.5, 0.5, -0.5);     
	glVertex3f(-0.5, 0.5, -0.5);      
	glVertex3f(-0.5, -0.5, -0.5);      
	glEnd();
	glPopMatrix();
	//glLineWidth(3);
	//glBegin(GL_LINE_STRIP);
	//for (double t = 0; t <= 1.0001; t += 0.01)
	//{
	//	double PP[3];
	//	PP[0] = Bez(P00[0], P11[0], P22[0], P33[0], t);
	//	PP[1] = Bez(P00[1], P11[1], P22[1], P33[1], t);
	//	PP[2] = Bez(P00[2], P11[2], P22[2], P33[2], t);
	//	glVertex3dv(PP); //Рисуем точку P
	//}
	//glEnd();
	//glLineWidth(1); //возвращаем ширину линии = 1
	//glBegin(GL_LINE_STRIP);
	//glVertex3d(P00[0], P00[1], P00[2]);
	//glVertex3d(P11[0], P11[1], P11[2]);
	//glVertex3d(P22[0], P22[1], P22[2]);
	//glVertex3d(P33[0], P33[1], P33[2]);
	//glEnd(); 
	//glLineWidth(3);
	//glBegin(GL_LINE_STRIP);
	//for (double t = 0; t <= 1.0001; t += 0.01)
	//{
	//	double Erm[3];
	//	Erm[0] = Ermit(Erm0[0], Erm1[0], Erm2[0], Erm3[0], t);
	//	Erm[1] = Ermit(Erm0[1], Erm1[1], Erm2[1], Erm3[1], t);
	//	Erm[2] = Ermit(Erm0[2], Erm1[2], Erm2[2], Erm3[2], t);
	//	glVertex3dv(Erm); 
	//}
	//glEnd();
	//glLineWidth(1); //возвращаем ширину линии = 1
	//glBegin(GL_LINE_STRIP);
	//glVertex3d(Erm0[0], Erm0[1], Erm0[2]);
	//glVertex3d(Erm2[0], Erm2[1], Erm2[2]);
	//glEnd();
	//glBegin(GL_LINE_STRIP);
	//glVertex3d(Erm1[0], Erm1[1], Erm1[2]);
	//glVertex3d(Erm3[0], Erm3[1], Erm3[2]);
	//glEnd();

	//glLineWidth(3);
	//glBegin(GL_LINE_STRIP);
	//for (double t = 0; t <= 1.0001; t += 0.01)
	//{
	//	double Ermm[3];
	//	Ermm[0] = Ermit(Erm00[0], Erm11[0], Erm22[0], Erm33[0], t);
	//	Ermm[1] = Ermit(Erm00[1], Erm11[1], Erm22[1], Erm33[1], t);
	//	Ermm[2] = Ermit(Erm00[2], Erm11[2], Erm22[2], Erm33[2], t);
	//	glVertex3dv(Ermm);
	//}
	//glEnd();
	//glLineWidth(1); //возвращаем ширину линии = 1
	//glBegin(GL_LINE_STRIP);
	//glVertex3d(Erm00[0], Erm00[1], Erm00[2]);
	//glVertex3d(Erm22[0], Erm22[1], Erm22[2]);
	//glEnd();
	//glBegin(GL_LINE_STRIP);
	//glVertex3d(Erm11[0], Erm11[1], Erm11[2]);
	//glVertex3d(Erm33[0], Erm33[1], Erm33[2]);
	//glEnd();
	//t_max += delta_time / 5; //t_max становится = 1 за 5 секунд
	//if (t_max > 1) t_max = 0; //после обнуляется
	//double P1[] = { 0,0,0 }; //Наши точки
	//double P2[] = { -4,6,7 };
	//double P3[] = { 10,10,0 };
	//double P[3];
	//glBegin(GL_LINE_STRIP); //построим отрезки P1P2 и P2P3
	//glVertex3dv(P1);
	//glVertex3dv(P2);
	//glVertex3dv(P3);
	//glEnd();
	//glLineWidth(3); //ширина линии
	//glColor3d(0, 1, 0);
	//glBegin(GL_LINE_STRIP);
	//for (double t = 0; t <= t_max; t += 0.01)
	//{
	//	P[0] = f(P1[0], P2[0], P3[0], t);
	//	P[1] = f(P1[1], P2[1], P3[1], t);
	//	P[2] = f(P1[2], P2[2], P3[2], t);
	//	glVertex3dv(P); //Рисуем точку P
	//}
	//glEnd();
	//glColor3d(1, 0, 1);
	//glLineWidth(1); //возвращаем ширину линии = 1
	////нарисуем все точки
	//glPointSize(10);
	//glBegin(GL_POINTS);
	//glVertex3dv(P);
	//glEnd();
	//glColor3d(1, 0, 0);
	//glBegin(GL_POINTS);
	//glVertex3dv(P1);
	//glVertex3dv(P2);
	//glVertex3dv(P3);
	//glEnd();
}

