/*
Author: @Johanan Luo
Time: begin at 2018/2/10
Version:
This is a sysytem to solve linear problem in 3-D.
This is a new way about prune&&search
It takes linear time in theory.
Input:some lines like 5.5a1 + 6a2 +7 a3 >= b
Output:the optimal x's posotion
*/

//some header
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <algorithm>
#include <vector>
#include <time.h>
#include <string.h>
#include <limits.h>
#include <float.h>

typedef int BOOL;
#define TRUE  1
#define FALSE 0

#define PI 3.1415926

using namespace std;
unsigned int randomSeed = 0;
//This is line
struct flat {
	//a1x + a2y + a3z >= b
	double a1, a2, a3, b;
	double op;
	double VX, VY, VZ;
	bool symbol;
	int groupNumber;
};

// Object Function
struct Obfunc {
	//xo = r1x +r2y +r3z
	double r1, r2, r3;
};

// the Line of two flats
struct Line {
	//d1x +d2y =e
	double d1, d2, e;
	double DX, DY, DZ;
	struct flat *f1, *f2;
};
struct Pair
{
	int index;
	struct flat f1, f2;
	struct Line line;
	struct Line_2D linexy;
	bool symbol;
};

//two lines' vertex
struct Vertex
{
	double x, y,z;
};
struct Line_2D {
	// a1x + a2y >= b
	double a1, a2, b;
	double slope;
	bool symbol;
};

// Object Function
struct Objfunc_2D {
	// xd = c1x + c2y
	double c1, c2;
};

// Vertex_2D structure
struct Vertex_2D {
	double x, y;
};

// PAIR structure
struct Pair_2D {
	int index;
	struct Line_2D line1, line2;
	flat f1, f2;
	Vertex_2D v1;
	int indexPair1, indexPair2;
	bool symbol;
};

struct PairXY {
	int index;
	int indexOfP1;
	int indexOfP2;
	struct Line_2D line1, line2;
	Vertex_2D point;
	bool  symbol;
};

typedef struct Line_2D Line_2D;
typedef struct Objfunc_2D Objfunc_2D;
typedef struct Vertex_2D Vertex_2D;
typedef struct PairXY PairXY;

typedef struct flat flat;
typedef struct Line line;
typedef struct Obfunc obfunc;
typedef struct Pair pair;
typedef struct Vertex vertex;
typedef struct PairXY PairXY;
vector<struct flat> originalC;

struct Vertex_2D Answer;
struct Vertex ver_answer;
struct Line line_answer;
struct flat flat_answer;

vector<struct flat> OriginalC;
vector<struct Line_2D> lineC;
//make objective function conincides with x-y axis
double transformationAngle(struct Obfunc *obfunc, struct Line *line) {
	double angle;
	double obNx, obNy, obNz;
	if (obfunc->r3 > 0) {
		obNx = -obfunc->r1;
		obNy = -obfunc->r2;
		obNz = -obfunc->r3;
	}
	angle = atan(sqrt(obNx*obNx + obNy*obNy*obNy) / obNz);

	struct flat *flatXY = (struct flat*)malloc(sizeof(struct flat));
	flatXY->a1 = 0;
	flatXY->a2 = 0;
	flatXY->a3 = 1;
	flatXY->b = 0;
	Vvector(flatXY);

	struct flat *flat1 = (struct flat*)malloc(sizeof(struct flat));
	flat1->a1 = obfunc->r1;
	flat1->a2 = obfunc->r2;
	flat1->a3 = obfunc->r3;
	flat1->b = 0;
	Vvector(flat1);
	findLines(flat1, flatXY, line);

	free(flat1);
	free(flatXY);

	return angle;


}
bool transformation(struct flat *flat, struct flat *nflat, struct Obfunc *obfucn) {
	//double transformationSlope = -(object.r1 / object.r2);
	struct Line *line = (struct Line*)malloc(sizeof(struct Line));
	double angle = transformationAngle(obfucn, line);

	double symbol = sqrt(line->DX*line->DX + line->DY*line->DY + line->DZ*line->DZ);
	double a = line->DX / symbol;
	double b = line->DY / symbol;
	double c = line->DZ / symbol;
	double A = flat->a1;
	double B = flat->a2;
	double C = flat->a3;

	double Matrix[3][3];

	Matrix[0][0] = a*a + (1 - a*a)*cos(angle);
	Matrix[0][1] = a*b + (1 - cos(angle) + c*sin(angle));
	Matrix[0][2] = a*c + (1 - cos(angle) - b*sin(angle));
	Matrix[1][0] = a*b + (1 - cos(angle) - c*sin(angle));
	Matrix[1][1] = b*b + (1 - b*b)*cos(angle);
	Matrix[1][2] = b*c + (1 - cos(angle) + a*sin(angle));
	Matrix[2][0] = a*c + (1 - cos(angle) + b*sin(angle));
	Matrix[2][1] = b*c + (1 - cos(angle) - a*sin(angle));
	Matrix[2][2] = c*c + (1 - c*c)*cos(angle);

	nflat->b = flat->b;

	nflat->a1 = Matrix[0][0] * A + Matrix[0][1] * B + Matrix[0][2] * C;
	nflat->a2 = Matrix[1][0] * A + Matrix[1][1] * B + Matrix[1][2] * C;
	nflat->a3 = Matrix[2][0] * A + Matrix[2][1] * B + Matrix[2][2] * C;

	Vvector(nflat);
	nflat->symbol = true;
	free(line);
	return true;
}

void Vvector(struct flat *flat) {
	if (flat->a3 > 0) {
		flat->VX = flat->a1;
		flat->VY = flat->a2;
		flat->VZ = flat->a3;
	}
	else {
		flat->VX = -flat->a1;
		flat->VY = -flat->a2;
		flat->VZ = -flat->a3;
	}
}

//find the two flats' line
bool findLines(struct flat *flat1, struct flat *flat2, struct Line *l1) {
	if (abs(flat1->a1*flat2->a2 - flat1->a2*flat2->a1) < 1e-6 && abs(flat1->a3*flat2->a2 - flat1->a2*flat2->a3)) {
		l1 = NULL;
		return false;
	}
	l1->f1 = flat1;
	l1->f2 = flat2;
	l1->DX = flat1->VY*flat2->VZ - flat1->VY*flat1->VZ;
	l1->DY = flat1->VX*flat2->VZ - flat1->VX*flat1->VZ;
	l1->DZ = flat1->VX*flat2->VY - flat1->VX*flat1->VY;

	if (abs(l1->DZ>1e-6)) {
		l1->e = 0;
		l1->d1 = -((flat1->b*flat2->a2 - flat2->b*flat1->a2) / (flat1->a1*flat2->a2 - flat2->a1*flat1->a2));
		l1->d2 = -((flat1->b*flat2->a1 - flat2->b*flat1->a1) / (flat1->a2*flat2->a1 - flat1->a1*flat2->a2));
	}
	else {
		l1->e = ((flat1->b*flat2->a2 - flat1->a2*flat2->b) / (flat1->a3*flat2->a1 - flat1->a2*flat2->a3));
		double bT1 = flat1->b - l1->e*flat1->a3;
		double bT2 = flat2->b - l1->e*flat2->a3;
		if (abs(l1->DX - 0) < 1e-6) {
			l1->d1 = ((flat1->b*flat2->a3 - flat2->b*flat1->a3) / (flat1->a1*flat2->a3 - flat2->a1*flat1->a3));
			l1->d2 = 1;
		}
		else {
			l1->d1 = 0;
			l1->d2 = (flat1->b - bT1) / (flat1->a2);
		}
	}
	return true;
}

bool allocate(struct flat flats[],struct Obfunc *objfunc,int index, int *numG, int *numH, int *numE ) {
	(*numG) = (*numH) = (*numE) = 0;
	int i = 0;
		for (; i < index; i++) {
			transformation(&originalC[i], &flats[i], objfunc);
			if (abs(flats[i].a3) < 1e-6)
				(*numG)++;
			else if (flats[i].a3 < 0)
				(*numH)++;
			else if (flats[i].a3 > 0)
				(*numE)++;
	}
		return true;
}

//do some segemention
void JudgeGHE(struct flat F1[], struct flat F2[], struct flat F3[], struct flat flats[], int *numG, int *numH, int *numE ) {
	int g=0, h = 0,e = 0;
	int numC = (*numG) + (*numH) + (*numE);
	int i = 0;
	for (; i < numC; i++) {
		if (abs(flats[i].a3) < 1e-6) {
			//#g++;
			//a3z>a1x+a2y+b
			F1[g].a1 = -flats[i].a1 / flats[i].a3;
			F1[g].a2 = -flats[i].a2 / flats[i].a3;
			F1[g].a3 = 1;
			F1[g].b = flats[i].b / flats[i].a3;
			F1[g].VX = flats[i].VX;
			F1[g].VY = flats[i].VY;
			F1[g].VZ = flats[i].VZ;
			F1[g].symbol = true;
			g++;
		}
		else if (abs(flats[i].a3) < 0) {
			//#h++;
			//a3z<a1x+a2y+b
			F2[h].a1 = -flats[i].a1 / flats[i].a3;
			F2[h].a2 = -flats[i].a2 / flats[i].a3;
			F2[h].a3 = 1;
			F2[h].b = flats[i].b / flats[i].a3;
			F2[h].VX = flats[i].VX;
			F2[h].VY = flats[i].VY;
			F2[h].VZ = flats[i].VZ;	
			F2[g].symbol = true;
			h++;
		}
		else if (abs(flats[i].a3) > 0) {
			//#g++;
			//0<a1x+a2y+b
			F3[e].a1 = -flats[i].a1 / flats[i].a3;
			F3[e].a2 = -flats[i].a2 / flats[i].a3;
			F3[e].a3 = 1;
			F3[e].b = flats[i].b / flats[i].a3;
			F3[e].VX = flats[i].VX;
			F3[e].VY = flats[i].VY;
			F3[e].VZ = flats[i].VZ;
			e++;
		}
	}

	return;
}

bool Judgeparallel(struct flat *f1, struct flat *f2) {
	double judge1 = abs(f1->VX*f2->VY - f2->VX*f1->VY);
	double judge2 = abs(f1->VX*f2->VZ - f2->VX*f1->VZ);
	double judge3= abs(f1->VY*f2->VZ - f2->VY*f1->VZ);
	if (judge1 < 1e-6&&judge2 < 1e-6&&judge3 < 1e-6) 
		return true;
	
	else
		return false;
}
//find lines in xy-planes
bool intersection(struct Line_2D *l1, struct Line_2D *l2, struct Vertex_2D *v1) 
{
	if (abs(l1->a1*l2->a2 - l2->a1*l1->a2) < 1e-6)
	{
		v1 = NULL;
		return false;
	}
	v1->x = -(l1->b*l2->a2 - l2->b*l1 -> a2) / (l1->a1*l2->a2 - l2->a1*l1->a2);
	v1->y = (l1->b*l2->a1 - l2->b*l1->a1) / (l1->a2*l2->a1 - l1->a1*l2->a2);
	return true;

}
bool projectionXY(struct Line *line, struct Line_2D *linxy) {
	struct Vertex_2D *vxy = (struct Vertex_2D*)malloc(sizeof(struct Vertex_2D));
	vxy->x = line->d1;
	vxy->y = line->d2;
	double xprime = vxy->x + line->DX;
	double yprime = vxy->y + line->DY;
	if (abs(xprime - vxy->x) < 1e-6 && abs(yprime - vxy->y) < 1e-6)
		return false;
	else if (abs(xprime - vxy->x) < 1e-6) {
		linxy->a2 = 0;
		linxy ->a1 = 1;
		linxy->b = vxy->x;
		if ((xprime > vxy->x&&yprime < vxy->y) || (xprime<vxy->x&&yprime>vxy->y))
			linxy->slope = -FLT_MAX;
		else
			linxy->slope = FLT_MAX;
	}
	else if (abs(yprime - vxy->y) < 1e-6) {
		linxy->a1 = 0;
		linxy->a2 = 1;
		linxy->b = vxy->y;
		linxy->slope = 0;
	}
	else {
		linxy->slope = line->DX / line->DY;
		linxy->a1 = linxy->slope;
		linxy->a2 = 1;
		linxy->b = vxy->y-linxy->a1*vxy->x;

	}
	linxy->symbol = true;
	free(vxy);
	return true;
}
// Intersection Vertex_2D
/*bool Intersection(struct Line_2D *l1, struct Line_2D *l2, struct Vertex_2D *v1)
{
	if (abs(l1->a1 * l2->a2 - l2->a1 * l1->a2) < 1e-6)
	{
		v1 = NULL;
		return false;
	}
	v1->x = -(l1->b * l2->a2 - l2->b * l1->a2) / (l1->a1 * l2->a2 - l2->a1 * l1->a2);
	v1->y = (l1->b * l2->a1 - l2->b * l1->a1) / (l1->a2 * l2->a1 - l1->a1 * l2->a2);
	return true;
}*/

// Slope Line_2D
bool Slope(struct Line_2D *l)
{
	if (fabs(l->a2 - 0.0) < 1e-6)
	{
		if ((l->a1 > 0 && l->a2 < 0) || (l->a1 < 0 && l->a2 > 0))
		{
			l->slope = FLT_MAX;
		}
		else if ((l->a1 < 0 && l->a2 < 0) || (l->a1 > 0 && l->a2 > 0))
		{
			l->slope = -FLT_MAX;
		}
		else
		{
			l->slope = -l->a1 / l->a2;
		}
		return false;
	}
	l->slope = -l->a1 / l->a2;
	return true;
}

// Compare
int Compare(const void *a, const void *b)
{
	struct Line_2D *aa = (struct Line_2D *)a;
	struct Line_2D *bb = (struct Line_2D *)b;
	return ((aa->slope > bb->slope) ? 1 : -1);
}

// transformation - O(n)
bool transformation_2D(struct Line_2D Line_2Ds[], struct Objfunc_2D object, int index, int *numG, int *numH)
{
	double thetaArc = atan(-object.c1 / object.c2);
	double thetaDec = atan(-object.c1 / object.c2) * 180 / PI;

	int i;
	double a1Temp, a2Temp, bTemp;

	(*numG) = 0;
	(*numH) = 0;

	for (i = 0; i < index; i += 1) {
		a1Temp = originalC[i].a1;
		a2Temp = originalC[i].a2;
		bTemp = originalC[i].b;

		Line_2Ds[i].a1 = cos(thetaArc) * a1Temp + sin(thetaArc) * a2Temp;
		Line_2Ds[i].a2 = cos(thetaArc) * a2Temp - sin(thetaArc) * a1Temp;
		Line_2Ds[i].b = bTemp;

		if (Line_2Ds[i].a2 > 0) {
			(*numG)++;
		}
		else if (Line_2Ds[i].a2 < 0) {
			(*numH)++;
		}
		else {
			return false;
		}

		Slope(&Line_2Ds[i]);
		Line_2Ds[i].symbol = true;

	}

	if ((*numG) + (*numH) != index) {
		printf("Fatal Error at transformation()!\n");
		exit(-1);
	}

	return true;
}

// Separation - O(n)
bool Judge(struct Line_2D I1[], struct Line_2D I2[], struct Line_2D Line_2Ds[], int numG, int numH)
{
	int index = numG + numH;
	int i, g = 0, h = 0;
	for (i = 0; i < index; i++) {
		if (Line_2Ds[i].a2 > 0) {
			I1[g].a1 = -Line_2Ds[i].a1 / Line_2Ds[i].a2;
			I1[g].a2 = 1;
			I1[g].b = Line_2Ds[i].b / Line_2Ds[i].a2;
			Slope(&I1[g]);
			I1[g].slope = -I1[g].slope;
			I1[g].symbol = true;
			g++;
		}
		else if (Line_2Ds[i].a2 < 0) {
			I2[h].a1 = -Line_2Ds[i].a1 / Line_2Ds[i].a2;
			I2[h].a2 = 1;
			I2[h].b = Line_2Ds[i].b / Line_2Ds[i].a2;
			Slope(&I2[h]);
			I2[h].slope = -I2[h].slope;
			I2[h].symbol = true;
			h++;
		}
		else {
			return false;
		}
	}
	return true;
}



// Make pairs
bool MakePairs(struct flat f1[], struct flat f2[],struct flat f3[], struct Pair pairsG[], struct Pair pairsH[], struct Pair pairsE[],
				int numG, int numH, int numE, int *indexG, int *indexE)
{
	int g, h, e, gT, eT;
	(*indexG) = (*indexE) = 0;
	for (g = 0; g < numG; g++) {
		if (f1[g].symbol == false) {
			continue;
		}
		for (gT= g+1; gT < numG; gT++) {
			if (f1[gT].symbol == true) {
				break;
			}
		}
		if (gT == numG) break;
		if (Judgeparallel(&f1[g], &f1[gT])) {
			if (f1[g].b > f1[gT].b) {
				f1[gT].symbol = false;
			}
			else{
				f1[g].symbol = false;
			}
			g = gT - 1;
			continue;
		}
		pairsG[(*indexG)].index = (*indexG);
		pairsG[(*indexG)].f1 = f1[g];
		pairsG[(*indexG)].f2 = f1[gT];
		pairsG[(*indexG)].symbol = true;
		findLines(&(pairsG[(*indexG)].f1), &(pairsG[(*indexG)].f2), &(pairsG[(*indexG)].line));
		projectionXY(&(pairsG[(*indexG)].line), &(pairsG[(*indexG)].linexy));
		(*indexG)++;
		g++;
	}
	for (e = 0; e < numG; e++) {
		if (f3[e].symbol == false) {
			continue;
		}
		for (eT = e + 1; eT < numE; eT++)
		{
			if (f3[eT].symbol == true) {
				break;
			}

		}
		if (eT == numE) break;

		if (Judgeparallel(&f3[e], &f3[eT])) {
			if (f3[e].b > f3[eT].b) {
				f3[eT].symbol = false;
			}
			else {
				f3[e].symbol = false;
			}
			e = eT - 1;
			continue;
		}
		pairsE[(*indexE)].index = (*indexE);
		pairsE[(*indexE)].f1 = f3[e];
		pairsE[(*indexE)].f2 = f3[eT];
		pairsE[(*indexE)].symbol = true;
		findLines(&(pairsE[(*indexE)].f1), &(pairsE[(*indexE)].f2), &(pairsE[(*indexE)].line));
		projectionXY(&(pairsE[(*indexE)].line), &(pairsE[(*indexE)].linexy));
		(*indexE)++;
		e++;
		return true;
	}
}



// sg, Sg, sh, Sh
struct Vertex_2D *TestingLine_2D(struct Pair_2D pairsG[], struct Pair_2D pairsH[],
	struct Line_2D I1[], struct Line_2D I2[],
	int numG, int numH, int numDot,
	double *leftBound, double *rightBound)
{
	// Randomly choose a v1

	randomSeed += numG;
	srand((unsigned int)time(NULL));
	int index = (numDot == 0) ? 0 : (rand() % numDot);
	//int index = round ? 1 : 0;
	double xPrimeG = pairsG[index].v1.x;   // x' - xPrime
	double yPrimeG = pairsG[index].v1.y;
	double yPrimeH;


	struct Line_2D *sg = NULL;
	struct Line_2D *Sg = NULL;
	struct Line_2D *sh = NULL;
	struct Line_2D *Sh = NULL;

	vector<int> Line_2DsG;
	vector<int> Line_2DsH;

	// Finding g(x') and H(x')
	for (int i = 0; i < numG; i++) {
		if (I1[i].symbol == true) {
			if ((abs(yPrimeG - (I1[i].a1 * xPrimeG + I1[i].b)) > 1e-6 && yPrimeG < (I1[i].a1 * xPrimeG + I1[i].b)) || (sg == NULL || Sg == NULL)) {




				yPrimeG = I1[i].a1 * xPrimeG + I1[i].b;
				sg = &I1[i];
				Sg = &I1[i];
			}
		}
	}
	for (int i = 0; i < numH; i++) {
		if (I2[i].symbol == true) {
			if (sh == NULL || Sh == NULL) {
				sh = &I2[i];
				Sh = &I2[i];
				yPrimeH = I2[i].a1 * xPrimeG + I2[i].b;
			}
			else if (abs(yPrimeH - (I2[i].a1 * xPrimeG + I2[i].b)) > 1e-6 && yPrimeH > (I2[i].a1 * xPrimeG + I2[i].b)) {
				yPrimeH = I2[i].a1 * xPrimeG + I2[i].b;
				sh = &I2[i];
				Sh = &I2[i];
			}
		}
	}
	if (numH == 0) {
		yPrimeH = yPrimeG + 1000.0;
	}

	for (int i = 0; i < numG; i++) {
		double currentLine_2DValueG = I1[i].a1 * xPrimeG + I1[i].b;
		if (I1[i].symbol == false || abs(currentLine_2DValueG - yPrimeG) >= 1e-6) {
			continue;
		}

		if (I1[i].a1 < sg->a1) {
			sg = &I1[i];
		}
		if (I1[i].a1 > Sg->a1) {
			Sg = &I1[i];
		}
	}
	// Finding sh - min h(x') && Finding Sh - max h(x')
	for (int i = 0; i < numH; i++) {
		double currentLine_2DValueH = I2[i].a1 * xPrimeG + I2[i].b;
		if (I2[i].symbol == false || abs(currentLine_2DValueH - yPrimeH) >= 1e-6) {
			continue;
		}

		if (I2[i].a1 < sh->a1) {
			sh = &I2[i];
		}
		if (I2[i].a1 > Sh->a1) {
			Sh = &I2[i];
		}
	}

	// Is feasible
	if (abs(yPrimeG - yPrimeH) < 1e-6) {
		if (sg->a1 > 0 && sg->a1 >= Sh->a1) {
			// x* < x'
			if (sh != Sh) {
				sh->symbol = false;
			}
			Sg->symbol = false;
			(*rightBound) = xPrimeG;
			return NULL;
		}
		else if (Sg->a1 < 0 && Sg->a1 <= sh->a1) {
			// x* > x'
			if (sh != Sh) {
				Sh->symbol = false;
			}
			sg->symbol = false;
			(*leftBound) = xPrimeG;
			return NULL;
		}
		else {
			// x* = x'
			Answer.x = xPrimeG;
			Answer.y = yPrimeG;
			return &(Answer);
		}
	}
	else if (yPrimeG > yPrimeH) {   // infeasible
		if (sg->a1 > Sh->a1) {
			// x* < x'
			if (sh != Sh) {
				sh->symbol = false;
			}
			Sg->symbol = false;
			(*rightBound) = xPrimeG;
			return NULL;
		}
		else if (Sg->a1 < sh->a1) {
			// x* > x'
			if (sh != Sh) {
				Sh->symbol = false;
			}
			sg->symbol = false;
			(*leftBound) = xPrimeG;
			return NULL;
		}
		else if ((sg->a1 - Sh->a1) <= 0 && 0 <= (Sg->a1 - sh->a1)) {
			// no feasible
			printf("No feasible Answer!\n");
			exit(0);
			return NULL;
		}
	}
	else if (yPrimeG < yPrimeH) {   // feasible
		if (sg->a1 > 0) {
			// x* < x'
			if (sg != Sg) {
				Sg->symbol = false;
			}
			else {
				if (pairsG[index].line1.a1 < pairsG[index].line2.a1) {
					pairsG[index].line1.symbol = false;
				}
				else if (pairsG[index].line1.a1 > pairsG[index].line2.a1) {
					pairsG[index].line2.symbol = false;
				}
			}
			(*rightBound) = xPrimeG;
			return NULL;
		} 
		else if (Sg->a1 < 0) {
			// x* > x'
			if (sg != Sg) {
				sg->symbol = false;
			}
			else {
				if (pairsG[index].line1.a1 < pairsG[index].line2.a1) {
					pairsG[index].line1.symbol = false;
				}
				else if (pairsG[index].line1.a1 > pairsG[index].line2.a1) {
					pairsG[index].line2.symbol = false;
				}
			}
			(*leftBound) = xPrimeG;
			return NULL;
		}
		else if (sg->a1 <= 0 && 0 <= Sg->a1) {
			// x* = x'
			Answer.x = xPrimeG;
			Answer.y = yPrimeG;
			return &(Answer);
		}
	}
}

/*void Line_2DarProgramming(void)
{
	int indexRecord = 0;
	int numGRecord;
	int numHRecord;
	int indexPair;
	double leftBound, rightBound;
	double aTemp, bTemp, cTemp;
	bool judge = false;
	struct Objfunc_2D object;

	//int round = 0;

	while (1) {
		scanf("%lf%lf%lf", &aTemp, &bTemp, &cTemp);
		if (aTemp == 0.0 && bTemp == 0.0 && cTemp == 0.0) {
			break;
		}
		struct Line_2D Line_2DTemp;
		Line_2DTemp.a1 = aTemp;
		Line_2DTemp.a2 = bTemp;
		Line_2DTemp.b = cTemp;
		//originalC.push_back(Line_2DTemp);
		indexRecord++;
	}
	scanf("%lf%lf", &object.c1, &object.c2);
	scanf("%lf%lf", &leftBound, &rightBound);

	struct Line_2D *Line_2Ds = (struct Line_2D *)malloc(indexRecord * sizeof(struct Line_2D));
	struct Line_2D *I1 = (struct Line_2D *)malloc(indexRecord * sizeof(struct Line_2D));
	struct Line_2D *I2 = (struct Line_2D *)malloc(indexRecord * sizeof(struct Line_2D));
	struct Pair *pairG = (struct Pair *)malloc(indexRecord * sizeof(struct Pair) / 2);
	struct Pair *pairH = (struct Pair *)malloc(indexRecord * sizeof(struct Pair) / 2);
	struct Vertex_2D *sln = NULL;

	judge = transformation_2D(Line_2Ds, object, indexRecord, &numGRecord, &numHRecord);
	if (judge == false) {
		printf("Fatal Error at Line_2DarProgramming() - transformation()!\n");
		exit(-1);
	}

	judge = Judge(I1, I2, Line_2Ds, numGRecord, numHRecord);
	if (judge == false) {
		printf("Fatal Error at Line_2DarProgramming() - Judge()!\n");
		exit(-1);
	}

	while (1) {
		judge = MakePairs(I1, I2, pairG, pairH, numGRecord, numHRecord, &indexPair, leftBound, rightBound);
		if (judge == false) {
			printf("Fatal Error at Line_2DarProgramming() - MakePairs()!\n");
			exit(-1);
		}

		sln = TestingLine_2D(pairG, pairH, I1, I2, numGRecord, numHRecord, indexPair, &leftBound, &rightBound);
		if (sln != NULL) {
			break;
		}
	}

	printf("The optimal answer under these constraints is: %lf %lf", sln->x, sln->y);

	return;

}

double findMaxMin(struct flat *flat) {
	int j = 0;
	//g(x,y)=MAX
	//h(x,y)=Min
	//e(x,y)=MAX
	//fx=gx-fx,ex;
	double gx = 0, hx = 0, ex = 0, fx;
	double tempG, tempH, tempE;
	while (1) {
		if (flat[j].groupNumber = 1) {
			tempG = flat[j].b;
			if (gx < flat[j].b) {
				gx = flat[j].b;
			}
		}
		else if (flat[j].groupNumber = 2) {
			//find min
			tempH = flat[j].b;
			if (hx > flat[j].b) {
				hx = flat[j].b;
			}
		}
		else {
			//find max
			tempE = flat[j].b;
			if (ex > flat[j].b) {
				ex = flat[j].b;
			}
		}
		j++;
		if (flat[j].groupNumber != 0 || 1 || 2)
			break;
	}
	fx = gx - hx;
}
// to acheive the line rotaotion to judge
bool rotation_3D(struct Obfunc *objfun, struct Line *line, struct flat *flat, struct flat *nflat) {
	// you should achieve 3D first
	
}*/

// judge the line and drop some flats



// the main function
double LP_3d() {

}

// test a line and judge the correctness

bool Test(struct flat f1[], struct flat f2[], struct flat f3[], struct Pair pairsG[], struct Pair pairsE[], struct Pair pairsH[]
	, int numG, int numH, int numE, int PairG, int PairH, int PairE)
{
	vector<struct PairXY> pairLineG;
	vector<struct PairXY> pairLineH;
	vector<struct PairXY> pairLineE;

	int indexG = 0;
	int indexE = 0;
	int indexH = 0;
	int indexI = 0;

	pairLineG.clear();
	pairLineH.clear();
	pairLineE.clear();

	int g, gT;
	int e, eT;

	for (g = 0; g < (PairG); g++) {
		if (pairsG[g].symbol == false) {
			continue;
		}
		for (gT = g + 1; gT < PairG; gT++) {
			if (pairsG[gT].symbol == true)
				break;
		}
		if (gT == PairG) break;
		if (abs(pairsG[g].linexy.slope - pairsG[gT].linexy.slope) < 1e-6) {
			if (pairsG[g].linexy.b > pairsG[gT].linexy.b) {
				pairsG[gT].linexy.symbol = false;
			}
			else {
				pairsG[gT].linexy.symbol = true;
			}
			g = gT - 1;
			continue;
		}
		struct Vertex_2D *pG = (struct Vertex_2D *)malloc(sizeof(struct Vertex_2D));
		intersection(&(pairsG[g].linexy), &(pairsG[gT].linexy), pG);
		struct PairXY pairxyg;
		pairxyg.symbol = false;
		pairxyg.index = indexG;
		pairxyg.indexOfP1 = g;
		pairxyg.indexOfP2 = gT;
		pairxyg.line1 = pairsG[g].linexy;
		pairxyg.line2 = pairsG[gT].linexy;
		pairxyg.point.x = pG->x;
		pairxyg.point.y = pG->y;
		pairLineG.push_back(pairxyg);
		indexG++;
	}
	for (e = 0; e < PairE; e++) {
		if (pairsE[e].symbol == false) {
			continue;
		}
		for (eT = e + 1; eT < PairE; eT++) {
			if (pairsE[eT].symbol == true) {
				break;
			}
		}
		if (eT == PairE) break;

		if (abs(pairsE[e].linexy.slope - pairsE[eT].linexy.slope) < 1e-6) {
			if (pairsE[e].linexy.b > pairsE[eT].linexy.b) {
				pairsE[eT].linexy.symbol = false;
			}
			else {
				pairsE[e].linexy.symbol = false;
			}
			e = eT - 1;
			continue;
		}
		struct Vertex_2D *pE = (struct Vertex_2D *)malloc(sizeof(struct Vertex_2D));
		intersection(&(pairsE[e].linexy), &(pairsE[eT].linexy), pE);
		struct PairXY pairxye;
		pairxye.symbol = false;
		pairxye.index = indexE;
		pairxye.indexOfP1 = e;
		pairxye.indexOfP2 = eT;
		pairxye.line1 = pairsE[e].linexy;
		pairxye.line2 = pairsE[eT].linexy;
		pairxye.point.x = pE->x;
		pairxye.point.y = pE->y;
		pairLineE.push_back(pairxye);
		indexE++;
	}
	int i, iT;
	for (i = 0; i < lineC.size(); i++) {
		if (lineC[i].symbol == false || lineC[i].a2 < 0) {
			continue;
		}
		for (iT = i + 1; iT < lineC.size(); iT++) {
			if (lineC[i].symbol == true && lineC[i].a2 > 0) {
				break;
			}
		}
		if (iT == lineC.size()) break;
		if (abs(lineC[i].slope - lineC[iT].slope) < 1e-6) {
			if (lineC[i].b > lineC[iT].b) {
				lineC[iT].symbol = false;
			}
			else {
				lineC[i].symbol = false;
			}
			i = iT - 1;
			continue;
		}
		struct Vertex_2D *vb = (struct Vertex_2D *)malloc(sizeof(struct Vertex_2D));
		intersection(&(lineC[i]), &(lineC[iT]), vb);
		struct Line_2D line1, line2;

		line1.a1 = -lineC[i].a1 / lineC[i].a2;
		line1.a2 = 1;
		line1.b = lineC[i].b / lineC[i].a2;
		line1.symbol = true;
		line1.symbol = lineC[i].slope;

		line2.a1 = -lineC[iT].a1 / lineC[iT].a2;
		line2.a2 = 1;
		line2.b = lineC[iT].b / lineC[iT].a2;
		line2.symbol = true;
		line2.symbol = lineC[iT].slope;

		struct PairXY pb;
		pb.index = -indexI;
		pb.indexOfP1 = i;
		pb.indexOfP2 = iT;
		pb.line1 = line1;
		pb.line2 = line2;
		pb.point.x = vb->x;
		pb.point.y = vb->y;
		pb.symbol = false;
		pairLineG.push_back(pb);
	}

	srand((unsigned int)time(NULL));
	int index = ((pairLineG.size() + pairLineE.size()) == 0) ? 0 : (rand() % (pairLineG.size() + pairLineE.size()));

	double Originalxm;
	double Originalym;
	bool gTeF = true;
	if (index >= pairLineG.size()) {
		Originalxm = pairLineE[index - pairLineG.size()].point.x;
		Originalym = pairLineE[index - pairLineG.size()].point.y;
		gTeF = false;
	}
	else {
		Originalxm = pairLineE[index].point.x;
		Originalym = pairLineE[index].point.y;
		gTeF = true;
	}

	double gmax = -FLT_MAX, hmin = FLT_MAX, emax = -FLT_MAX;
	vector<struct flat *>gAfterTrans;
	vector<struct flat *>hAfterTrans;
	vector<struct flat *>eAfterTrans;
	vector<struct flat *>gAfterTransoptimal;
	vector<struct flat *>hAfterTransoptimal;
	vector<struct flat *>eAfterTransoptimal;
	for (int i = 0; i < numG; i++) {
		if (f1[i].symbol == false) {
			continue;
		}
		struct flat *temp = (struct flat *)malloc(sizeof(struct flat));
		temp->a1 = f1[i].a1;
		temp->a2 = f1[i].a2;
		temp->a3 = f1[i].a3;
		temp->b = f1[i].b-(f1[i].a1*(-Originalxm))-(f1[i].a2*(-Originalym)) ;
		Vvector(temp);
		temp->symbol = false;
		gAfterTrans.push_back(temp);

		if (gmax < temp->b) gmax = temp->b;
	}
	for (int i = 0; i < numH; i++) {
		if (f2[i].symbol == false) {
			continue;
		}
		struct flat *temp = (struct flat *)malloc(sizeof(struct flat));
		temp->a1 = f2[i].a1;
		temp->a2 = f2[i].a2;
		temp->a3 = f2[i].a3;
		temp->b = f2[i].b - (f2[i].a1*(-Originalxm)) - (f2[i].a2*(-Originalym));
		Vvector(temp);
		temp->symbol = false;
		gAfterTrans.push_back(temp);

		if (gmax < temp->b) hmin = temp->b;
	}
	for (int i = 0; i < numE; i++) {
		if (f3[i].symbol == false) {
			continue;
		}
		struct flat *temp = (struct flat *)malloc(sizeof(struct flat));
		temp->a1 = f3[i].a1;
		temp->a2 = f3[i].a2;
		temp->a3 = f3[i].a3;
		temp->b = f3[i].b - (f3[i].a1*(-Originalxm)) - (f3[i].a2*(-Originalym));
		Vvector(temp);
		temp->symbol = false;
		gAfterTrans.push_back(temp);

		if (gmax < temp->b) emax = temp->b;
	}

	double fx;
	fx = gmax - hmin;
	double xoptimal,yoptimal;
	//use the equations in the paper to judge the optimal's position.
	//it should be minde that The lamda should be minded.
	for (int j = 0; j < lineC.size(); j++) {
		if (f2[j].symbol = false) {
			if (f1[j].b - f2[j].b < 0) {
				if (gmax < 0 && gmax <= hmin&&emax <= 0) {
					//minimize g can get the optimal value
					f2[j].symbol = true;
				}
				else if (gmax > 0 && gmax - hmin >= 0 && emax >= 0) {
					f2[j].symbol =false;
				}
			}
			
		}
		for (int m = 0; m < lineC.sizeof(); m++) {
			if (f2[m].symbol = false) {
				if (f1[m].b - f2[m].b > 0) {
					if (gmax > 0 && gmax >= hmin&&emax >= 0) {
						//minimize g can get the optimal value
						f2[m].symbol = true;
					}
					else if (  gmax - hmin > 0 && emax < 0) {
						f2[m].symbol = false;
					}
				}
			
	}
		
	}
		for (int h = 0; h < lineC.sizeof();h++) {
			if (f2[h].symbol = false) {
				if (f1[h].b - f2[h].b > 0) {
					if ( gmax >= hmin&&emax > 0) {
						//minimize g can get the optimal value
						f2[h].symbol = true;
					}
				}
				
			}


}
//main function
int main(void) {

	LP_3d();
	struct flat *test = (struct flat *)malloc(sizeof(struct flat));
	test->a1 = test->a2 = test->b = 1;
	Vvector(test);

	struct flat *flat= (struct flat *)malloc(sizeof(struct flat));
	struct Line *line= (struct Line *)malloc(sizeof(struct Line));
	struct Obfunc *obfunc = (struct Obfunc *)malloc(sizeof(struct Obfunc));

	obfunc->r1 = obfunc->r2 = obfunc->r3 = 1;

	


	return 0;
}
