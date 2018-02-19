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
	double d1, d2 ,e;
	double DX, DY, DZ;
	struct flat *f1, *f2;
};
struct Pair
{
	int index;
	struct flat f1, f2;
	struct Line line;
	bool symbol;
};

//two lines' vertex
struct Vertex
{
	double x, y, z;
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
struct Pair {
	int index;
	Line_2D Line_2D1, Line_2D2;
	Vertex_2D v1;
	bool symbol;
};

typedef struct Line_2D Line_2D;
typedef struct Objfunc_2D Objfunc_2D;
typedef struct Vertex_2D Vertex_2D;

vector<struct Line_2D> originalC;


typedef struct flat flat;
typedef struct Line line;
typedef struct Obfunc obfunc;
typedef struct Pair pair;
typedef struct Vertex vertex;

vector<struct flat> originalC;

struct Vertex_2D Answer;
struct Vertex ver_answer;
struct Line line_answer;
struct flat flat_answer;
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
bool transformation(struct flat *flat, struct flat *nflat, struct Obfunc *obfucn, struct Line *line) {
	//double transformationSlope = -(object.r1 / object.r2);
	double angle = transformationAngle(obfucn, line);

	double symbol = sqrt(line->DX*line->DX + line->DY*line->DY + line->DZ*line->DZ);
	double a = line->DX/symbol;
	double b = line->DY/symbol;
	double c = line->DZ/symbol;
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
bool findLines(struct flat *flat1,struct flat *flat2, struct Line *l1) {
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

bool allocate(struct flat *flat) {
	int i = 0;
		while (1) {
			if (flat[i].op == '>='&&flat[i].b != 0) {
				flat[i].groupNumber = 1;
			}
			else if (flat[i].op == '<='&&flat[i].b != 0) {
				flat[i].groupNumber = 1;
			}
			else {
				flat[i].groupNumber = 3;
			}
			flat[i].symbol = false;
			i++;
	}
	
}

// Intersection Vertex_2D
bool Intersection(struct Line_2D *l1, struct Line_2D *l2, struct Vertex_2D *v1)
{
	if (abs(l1->a1 * l2->a2 - l2->a1 * l1->a2) < 1e-6)
	{
		v1 = NULL;
		return false;
	}
	v1->x = -(l1->b * l2->a2 - l2->b * l1->a2) / (l1->a1 * l2->a2 - l2->a1 * l1->a2);
	v1->y = (l1->b * l2->a1 - l2->b * l1->a1) / (l1->a2 * l2->a1 - l1->a1 * l2->a2);
	return true;
}

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
bool transformation(struct Line_2D Line_2Ds[], struct Objfunc_2D object, int index, int *numG, int *numH)
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
bool MakePairs(struct Line_2D I1[], struct Line_2D I2[],struct Pair pairsG[], struct Pair pairsH[],	int numG, int numH, int *index,	double leftBound, double rightBound)
{
	int g, h, gtemp;
	(*index) = 0;
	for (g = 0; g < numG; g += 1) {
		// drop
		if (I1[g].symbol == false) {
			continue;
		}
		for (gtemp = g + 1; gtemp < numG; gtemp++) {
			if (I1[gtemp].symbol == true) {
				break;
			}
		}
		if (gtemp == numG) break;

		if (abs(I1[g].slope - I1[gtemp].slope) < 1e-6) {
			if (I1[g].b > I1[gtemp].b) {
				I1[gtemp].symbol = false;
			}
			else {
				I1[g].symbol = false;
			}
			g = gtemp - 1;
			continue;
		}
		struct Vertex_2D *p = (struct Vertex_2D *)malloc(sizeof(struct Vertex_2D));
		Intersection(&I1[g], &I1[gtemp], p);
		if (p->x < leftBound || p->x > rightBound) {
			if (abs(I1[g].slope) > abs(I1[gtemp].slope)) {
				I1[gtemp].symbol = false;
			}
			else if (abs(I1[gtemp].slope) < abs(I1[gtemp].slope)) {
				I1[g].symbol = false;
			}
			g = gtemp - 1;
			continue;
		}
		pairsG[(*index)].index = (*index);
		pairsG[(*index)].Line_2D1 = I1[g];
		pairsG[(*index)].Line_2D2 = I1[gtemp];
		pairsG[(*index)].v1.x = p->x; pairsG[(*index)].v1.y = p->y;

		(*index)++;
		g++;
	}

	return true;
}

// sg, Sg, sh, Sh
struct Vertex_2D *TestingLine_2D(struct Pair pairsG[], struct Pair pairsH[],
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
				if (pairsG[index].Line_2D1.a1 < pairsG[index].Line_2D2.a1) {
					pairsG[index].Line_2D2.symbol = false;
				}
				else if (pairsG[index].Line_2D1.a1 > pairsG[index].Line_2D2.a1) {
					pairsG[index].Line_2D1.symbol = false;
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
				if (pairsG[index].Line_2D1.a1 < pairsG[index].Line_2D2.a1) {
					pairsG[index].Line_2D1.symbol = false;
				}
				else if (pairsG[index].Line_2D1.a1 > pairsG[index].Line_2D2.a1) {
					pairsG[index].Line_2D2.symbol = false;
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

void Line_2DarProgramming(void)
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
		originalC.push_back(Line_2DTemp);
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

	judge = transformation(Line_2Ds, object, indexRecord, &numGRecord, &numHRecord);
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
bool rotation_3D() {
	// you should achieve 3D first

}

// judge the line and drop some flats
bool judge(struct flat *flat  ) {
	
}


// the main function
double LP_3d() {

}



//main function
int main(void) {

	LP_3d();

	return 0;
}
