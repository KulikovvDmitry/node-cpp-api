  // contents of addon_source.cc

  // This is a basic addon - a binding.gyp file 
  // would need to include this file as it's source.

  #include <nan.h>
  using namespace std;
  using namespace Nan;
  using namespace v8;

#include <iostream>
#include <stdio.h>
#include <vector>

class Point
{
	public:
		double x, y;
		bool color;
		Point(){}

		void setCoords(double x, double y)
		{
			this->x = x;
			this->y = y;
		}

		void setColor(bool color)
		{
			this->color = color;
		}
};

class Compute
{
	public: Compute(int points)
	{
		this->points = points;
		masX = new double[points];
		masY = new double[points];

		masA = new double*[points];
		masB = new double*[points];
		masC = new double*[points];
		masDD = new double*[points];
		masomega = new double*[points];
		mask = new double*[points];
		masV = new double*[points];
		masS = new double*[points];
		masL = new double*[points];
		masc = new double*[points];
		masN = new double*[points];
		mash = new double*[points];
		massigma = new double*[points];
		masa = new double*[points];
		masb = new double*[points];
		masp = new double*[points];
		masalpha = new double*[points];
		masbetta = new double*[points];

		q0 = new double*[points];
		q1 = new double*[points];
		q2 = new double*[points];
		q3 = new double*[points];
		q4 = new double*[points];
		q5 = new double*[points];
		q6 = new double*[points];
		det0 = new double*[points];
		det1 = new double*[points];
		det2 = new double*[points];
		det3 = new double*[points];
		det4 = new double*[points];
		det5 = new double*[points];
		det6 = new double*[points];
		outPoints = new Point*[points];

		for (int i = 0; i < points; i++)
		{
			masA[i] = new double[points];
			masB[i] = new double[points];
			masC[i] = new double[points];
			masDD[i] = new double[points];
			masomega[i] = new double[points];
			mask[i] = new double[points];
			masV[i] = new double[points];
			masS[i] = new double[points];
			masL[i] = new double[points];
			masc[i] = new double[points];
			masN[i] = new double[points];
			mash[i] = new double[points];
			massigma[i] = new double[points];
			masa[i] = new double[points];
			masb[i] = new double[points];
			masp[i] = new double[points];
			masalpha[i] = new double[points];
			masbetta[i] = new double[points];

			q0[i] = new double[points];
			q1[i] = new double[points];
			q2[i] = new double[points];
			q3[i] = new double[points];
			q4[i] = new double[points];
			q5[i] = new double[points];
			q6[i] = new double[points];
			det0[i] = new double[points];
			det1[i] = new double[points];
			det2[i] = new double[points];
			det3[i] = new double[points];
			det4[i] = new double[points];
			det5[i] = new double[points];
			det6[i] = new double[points];
			outPoints[i] = new Point[points];
		}
	}

	public:~Compute()
	{
		delete[] masX;
		delete[] masY;

		for (int i = 0; i < points; i++)
		{
			delete[] masA[i];
			delete[] masB[i];
			delete[] masC[i];
			delete[] masDD[i];
			delete[] masomega[i];
			delete[] mask[i];
			delete[] masV[i];
			delete[] masS[i];
			delete[] masL[i];
			delete[] masc[i];
			delete[] masN[i];
			delete[] mash[i];
			delete[] massigma[i];
			delete[] masa[i];
			delete[] masb[i];
			delete[] masp[i];
			delete[] masalpha[i];
			delete[] masbetta[i];

			delete[] q0[i];
			delete[] q1[i];
			delete[] q2[i];
			delete[] q3[i];
			delete[] q4[i];
			delete[] q5[i];
			delete[] q6[i];
			delete[] det0[i];
			delete[] det1[i];
			delete[] det2[i];
			delete[] det3[i];
			delete[] det4[i];
			delete[] det5[i];
			delete[] det6[i];
		}

		delete[] masA;
		delete[] masB;
		delete[] masC;
		delete[] masDD;
		delete[] masomega;
		delete[] mask;
		delete[] masV;
		delete[] masS;
		delete[] masL;
		delete[] masc;
		delete[] masN;
		delete[] mash;
		delete[] massigma;
		delete[] masa;
		delete[] masb;
		delete[] masp;
		delete[] masalpha;
		delete[] masbetta;

		delete[] q0;
		delete[] q1;
		delete[] q2;
		delete[] q3;
		delete[] q4;
		delete[] q5;
		delete[] q6;
		delete[] det0;
		delete[] det1;
		delete[] det2;
		delete[] det3;
		delete[] det4;
		delete[] det5;
		delete[] det6;
	}

	public:
		double* masX;
		double* masY;
		Point** outPoints;

		double** masA;
		double** masB;
		double** masC;
		double** masDD;
		double** masomega;
		double** mask;
		double** masV;
		double** masS;
		double** masL;
		double** masc;
		double** masN;
		double** mash;
		double** massigma;
		double** masa;
		double** masb;
		double** masp;
		double** masalpha;
		double** masbetta;

		double** q0;
		double** q1;
		double** q2;
		double** q3;
		double** q4;
		double** q5;
		double** q6;
		double** det0;
		double** det1;
		double** det2;
		double** det3;
		double** det4;
		double** det5;
		double** det6;

		int points;
		double A;
		double B;
		double C;
		double DD;
		double omega;
		double k;
		double V;
		double S;
		double L;
		double c;
		double N;
		double h;
		double sigma;
		double a;
		double b;
		double p;
		double alpha;
		double betta;
		int x, y;
		double koeffXMin, koeffXMax, koeffYMin, koeffYMax;
};

using namespace std;

void formulas(Compute* compute, double*** masI)
{
	int tr = 0;
	int fls = 0;
	for (int i = 0; i < compute->points; i++)
	{
		for (int j = 0; j < compute->points; j++)
		{
			compute->q6[i][j] = (-masI[0][i][j] * masI[1][i][j]) + (masI[3][i][j] * masI[3][i][j]);
			compute->q5[i][j] = ((-masI[0][i][j] * masI[11][i][j]) - (((masI[0x11][i][j] * masI[6][i][j]) * masI[0][i][j]) * masI[1][i][j])) + (((masI[0x11][i][j] * masI[6][i][j]) * masI[3][i][j]) * masI[3][i][j]);
			compute->q4[i][j] = ((((((((((((((((masI[3][i][j] * masI[3][i][j]) * masI[0x10][i][j]) * masI[6][i][j]) * masI[6][i][j]) - (((masI[0][i][j] * masI[9][i][j]) * masI[9][i][j]) * masI[13][i][j])) - ((masI[15][i][j] * masI[10][i][j]) * masI[1][i][j])) + ((masI[8][i][j] * masI[10][i][j]) * masI[1][i][j])) - (((2.0 * masI[3][i][j]) * masI[9][i][j]) * masI[10][i][j])) - (masI[14][i][j] * masI[0][i][j])) - (masI[5][i][j] * masI[1][i][j])) - (((masI[2][i][j] * masI[2][i][j]) * masI[4][i][j]) * masI[4][i][j])) - ((((masI[0][i][j] * masI[0x10][i][j]) * masI[6][i][j]) * masI[6][i][j]) * masI[1][i][j])) - ((((2.0 * masI[8][i][j]) * masI[10][i][j]) * masI[12][i][j]) * masI[1][i][j])) + ((((2.0 * masI[3][i][j]) * masI[8][i][j]) * masI[13][i][j]) * masI[9][i][j])) + ((((2.0 * masI[3][i][j]) * masI[12][i][j]) * masI[10][i][j]) * masI[9][i][j])) - (((masI[0x11][i][j] * masI[6][i][j]) * masI[0][i][j]) * masI[11][i][j])) - (((masI[8][i][j] * masI[8][i][j]) * masI[13][i][j]) * masI[1][i][j]);
			compute->q3[i][j] = (((((((((((((((((((((-masI[14][i][j] * masI[0][i][j]) * masI[0x10][i][j]) * masI[6][i][j]) * masI[9][i][j]) + ((((masI[14][i][j] * masI[8][i][j]) * masI[0x10][i][j]) * masI[6][i][j]) * masI[3][i][j])) - (((((masI[0x11][i][j] * masI[6][i][j]) * masI[0][i][j]) * masI[9][i][j]) * masI[9][i][j]) * masI[13][i][j])) - ((((masI[0x11][i][j] * masI[6][i][j]) * masI[15][i][j]) * masI[10][i][j]) * masI[1][i][j])) + ((((masI[0x11][i][j] * masI[6][i][j]) * masI[8][i][j]) * masI[10][i][j]) * masI[1][i][j])) - (((((2.0 * masI[0x11][i][j]) * masI[6][i][j]) * masI[3][i][j]) * masI[9][i][j]) * masI[10][i][j])) - (((((masI[0x11][i][j] * masI[6][i][j]) * masI[8][i][j]) * masI[8][i][j]) * masI[13][i][j]) * masI[1][i][j])) - ((masI[15][i][j] * masI[10][i][j]) * masI[11][i][j])) + ((masI[8][i][j] * masI[10][i][j]) * masI[11][i][j])) - ((((((2.0 * masI[0x11][i][j]) * masI[6][i][j]) * masI[8][i][j]) * masI[10][i][j]) * masI[12][i][j]) * masI[1][i][j])) + ((((((2.0 * masI[0x11][i][j]) * masI[6][i][j]) * masI[3][i][j]) * masI[8][i][j]) * masI[13][i][j]) * masI[9][i][j])) + ((((((2.0 * masI[0x11][i][j]) * masI[6][i][j]) * masI[3][i][j]) * masI[12][i][j]) * masI[10][i][j]) * masI[9][i][j])) - (masI[5][i][j] * masI[11][i][j])) - (((2.0 * masI[14][i][j]) * masI[6][i][j]) * masI[3][i][j])) - ((((masI[0][i][j] * masI[0x10][i][j]) * masI[6][i][j]) * masI[6][i][j]) * masI[11][i][j])) - ((((2.0 * masI[8][i][j]) * masI[10][i][j]) * masI[12][i][j]) * masI[11][i][j])) - (((masI[0x11][i][j] * masI[6][i][j]) * masI[5][i][j]) * masI[1][i][j])) - (((((masI[0x11][i][j] * masI[6][i][j]) * masI[2][i][j]) * masI[2][i][j]) * masI[4][i][j]) * masI[4][i][j])) - (((masI[8][i][j] * masI[8][i][j]) * masI[13][i][j]) * masI[11][i][j]);
			compute->q2[i][j] = ((((((((((((((((((((((((((((((((((((masI[9][i][j] * masI[10][i][j]) * masI[12][i][j]) * masI[0x10][i][j]) * masI[6][i][j]) * masI[6][i][j]) * masI[3][i][j]) - ((((((masI[8][i][j] * masI[10][i][j]) * masI[12][i][j]) * masI[0x10][i][j]) * masI[6][i][j]) * masI[6][i][j]) * masI[1][i][j])) + (((((masI[14][i][j] * masI[8][i][j]) * masI[0x10][i][j]) * masI[6][i][j]) * masI[2][i][j]) * masI[4][i][j])) - ((((((2.0 * masI[0x11][i][j]) * masI[6][i][j]) * masI[8][i][j]) * masI[10][i][j]) * masI[12][i][j]) * masI[11][i][j])) - ((masI[14][i][j] * masI[15][i][j]) * masI[10][i][j])) + (((((masI[9][i][j] * masI[9][i][j]) * masI[10][i][j]) * masI[10][i][j]) * masI[12][i][j]) * masI[12][i][j])) - (((masI[5][i][j] * masI[9][i][j]) * masI[9][i][j]) * masI[13][i][j])) - (((((2.0 * masI[9][i][j]) * masI[9][i][j]) * masI[10][i][j]) * masI[10][i][j]) * masI[12][i][j])) - (((masI[14][i][j] * masI[8][i][j]) * masI[8][i][j]) * masI[13][i][j])) + ((masI[14][i][j] * masI[8][i][j]) * masI[10][i][j])) - (((((masI[15][i][j] * masI[10][i][j]) * masI[0x10][i][j]) * masI[6][i][j]) * masI[6][i][j]) * masI[1][i][j])) + (((((masI[8][i][j] * masI[10][i][j]) * masI[0x10][i][j]) * masI[6][i][j]) * masI[6][i][j]) * masI[1][i][j])) - ((((((2.0 * masI[3][i][j]) * masI[0x10][i][j]) * masI[6][i][j]) * masI[6][i][j]) * masI[9][i][j]) * masI[10][i][j])) - ((((masI[0x11][i][j] * masI[6][i][j]) * masI[15][i][j]) * masI[10][i][j]) * masI[11][i][j])) + ((((masI[0x11][i][j] * masI[6][i][j]) * masI[8][i][j]) * masI[10][i][j]) * masI[11][i][j])) - (((((masI[0x11][i][j] * masI[6][i][j]) * masI[8][i][j]) * masI[8][i][j]) * masI[13][i][j]) * masI[11][i][j])) - (((((masI[0x11][i][j] * masI[6][i][j]) * masI[6][i][j]) * masI[0][i][j]) * masI[9][i][j]) * masI[13][i][j])) + (((((masI[0x11][i][j] * masI[6][i][j]) * masI[6][i][j]) * masI[3][i][j]) * masI[8][i][j]) * masI[13][i][j])) + (((((masI[0x11][i][j] * masI[6][i][j]) * masI[6][i][j]) * masI[3][i][j]) * masI[12][i][j]) * masI[10][i][j])) - (masI[14][i][j] * masI[5][i][j])) - ((((masI[15][i][j] * masI[10][i][j]) * masI[9][i][j]) * masI[9][i][j]) * masI[13][i][j])) - ((((masI[8][i][j] * masI[10][i][j]) * masI[9][i][j]) * masI[9][i][j]) * masI[13][i][j])) - (((masI[0x11][i][j] * masI[6][i][j]) * masI[5][i][j]) * masI[11][i][j])) - ((((2.0 * masI[14][i][j]) * masI[6][i][j]) * masI[2][i][j]) * masI[4][i][j])) + (((((2.0 * masI[6][i][j]) * masI[6][i][j]) * masI[9][i][j]) * masI[13][i][j]) * masI[3][i][j])) - (((((2.0 * masI[6][i][j]) * masI[6][i][j]) * masI[8][i][j]) * masI[13][i][j]) * masI[1][i][j])) - (((((2.0 * masI[6][i][j]) * masI[6][i][j]) * masI[12][i][j]) * masI[10][i][j]) * masI[1][i][j])) - ((((masI[14][i][j] * masI[0][i][j]) * masI[0x10][i][j]) * masI[6][i][j]) * masI[6][i][j])) - ((((2.0 * masI[14][i][j]) * masI[8][i][j]) * masI[10][i][j]) * masI[12][i][j])) - ((((masI[5][i][j] * masI[0x10][i][j]) * masI[6][i][j]) * masI[6][i][j]) * masI[1][i][j])) - ((((((masI[2][i][j] * masI[2][i][j]) * masI[4][i][j]) * masI[4][i][j]) * masI[0x10][i][j]) * masI[6][i][j]) * masI[6][i][j])) + (((masI[9][i][j] * masI[9][i][j]) * masI[10][i][j]) * masI[10][i][j]);
			compute->q1[i][j] = ((((((((((((((((((((((((-masI[14][i][j] * masI[5][i][j]) * masI[0x10][i][j]) * masI[6][i][j]) * masI[9][i][j]) + ((((2.0 * masI[14][i][j]) * masI[6][i][j]) * masI[9][i][j]) * masI[10][i][j])) - (((((masI[15][i][j] * masI[10][i][j]) * masI[0x10][i][j]) * masI[6][i][j]) * masI[6][i][j]) * masI[11][i][j])) + (((((masI[8][i][j] * masI[10][i][j]) * masI[0x10][i][j]) * masI[6][i][j]) * masI[6][i][j]) * masI[11][i][j])) - (((((2.0 * masI[6][i][j]) * masI[6][i][j]) * masI[8][i][j]) * masI[13][i][j]) * masI[11][i][j])) - (((((2.0 * masI[6][i][j]) * masI[6][i][j]) * masI[12][i][j]) * masI[10][i][j]) * masI[11][i][j])) - (((((masI[0x11][i][j] * masI[6][i][j]) * masI[5][i][j]) * masI[9][i][j]) * masI[9][i][j]) * masI[13][i][j])) - (((((((2.0 * masI[0x11][i][j]) * masI[6][i][j]) * masI[9][i][j]) * masI[9][i][j]) * masI[10][i][j]) * masI[10][i][j]) * masI[12][i][j])) + (((((((masI[0x11][i][j] * masI[6][i][j]) * masI[9][i][j]) * masI[9][i][j]) * masI[10][i][j]) * masI[10][i][j]) * masI[12][i][j]) * masI[12][i][j])) - (((((2.0 * masI[14][i][j]) * masI[6][i][j]) * masI[8][i][j]) * masI[13][i][j]) * masI[9][i][j])) - (((((2.0 * masI[14][i][j]) * masI[6][i][j]) * masI[12][i][j]) * masI[10][i][j]) * masI[9][i][j])) - (((((masI[14][i][j] * masI[15][i][j]) * masI[10][i][j]) * masI[0x10][i][j]) * masI[6][i][j]) * masI[9][i][j])) - ((((((masI[8][i][j] * masI[10][i][j]) * masI[12][i][j]) * masI[0x10][i][j]) * masI[6][i][j]) * masI[6][i][j]) * masI[11][i][j])) + ((((((2.0 * masI[6][i][j]) * masI[6][i][j]) * masI[9][i][j]) * masI[13][i][j]) * masI[2][i][j]) * masI[4][i][j])) - ((((((masI[0x11][i][j] * masI[6][i][j]) * masI[6][i][j]) * masI[2][i][j]) * masI[4][i][j]) * masI[8][i][j]) * masI[13][i][j])) - ((((((masI[0x11][i][j] * masI[6][i][j]) * masI[6][i][j]) * masI[2][i][j]) * masI[4][i][j]) * masI[12][i][j]) * masI[10][i][j])) - ((((((masI[0x11][i][j] * masI[6][i][j]) * masI[15][i][j]) * masI[10][i][j]) * masI[9][i][j]) * masI[9][i][j]) * masI[13][i][j])) - ((((((masI[0x11][i][j] * masI[6][i][j]) * masI[8][i][j]) * masI[10][i][j]) * masI[9][i][j]) * masI[9][i][j]) * masI[13][i][j])) - ((((((masI[14][i][j] * masI[8][i][j]) * masI[10][i][j]) * masI[12][i][j]) * masI[0x10][i][j]) * masI[6][i][j]) * masI[9][i][j])) - ((((masI[5][i][j] * masI[0x10][i][j]) * masI[6][i][j]) * masI[6][i][j]) * masI[11][i][j])) + (((((masI[0x11][i][j] * masI[6][i][j]) * masI[9][i][j]) * masI[9][i][j]) * masI[10][i][j]) * masI[10][i][j])) + (((((((masI[9][i][j] * masI[10][i][j]) * masI[12][i][j]) * masI[0x10][i][j]) * masI[6][i][j]) * masI[6][i][j]) * masI[2][i][j]) * masI[4][i][j]);
			compute->q0[i][j] = (((((((((((((((((-masI[14][i][j] * masI[8][i][j]) * masI[10][i][j]) * masI[12][i][j]) * masI[0x10][i][j]) * masI[6][i][j]) * masI[6][i][j]) - ((((((masI[0x11][i][j] * masI[6][i][j]) * masI[6][i][j]) * masI[15][i][j]) * masI[10][i][j]) * masI[9][i][j]) * masI[13][i][j])) - ((((masI[14][i][j] * masI[5][i][j]) * masI[0x10][i][j]) * masI[6][i][j]) * masI[6][i][j])) + ((((((masI[9][i][j] * masI[9][i][j]) * masI[10][i][j]) * masI[10][i][j]) * masI[0x10][i][j]) * masI[6][i][j]) * masI[6][i][j])) - (((((masI[14][i][j] * masI[15][i][j]) * masI[10][i][j]) * masI[0x10][i][j]) * masI[6][i][j]) * masI[6][i][j])) + (((((masI[14][i][j] * masI[8][i][j]) * masI[10][i][j]) * masI[0x10][i][j]) * masI[6][i][j]) * masI[6][i][j])) - (((((2.0 * masI[14][i][j]) * masI[6][i][j]) * masI[6][i][j]) * masI[8][i][j]) * masI[13][i][j])) - (((((2.0 * masI[14][i][j]) * masI[6][i][j]) * masI[6][i][j]) * masI[12][i][j]) * masI[10][i][j])) - (((((((masI[9][i][j] * masI[9][i][j]) * masI[10][i][j]) * masI[10][i][j]) * masI[12][i][j]) * masI[0x10][i][j]) * masI[6][i][j]) * masI[6][i][j])) - ((((((2.0 * masI[6][i][j]) * masI[6][i][j]) * masI[9][i][j]) * masI[9][i][j]) * masI[13][i][j]) * masI[10][i][j])) - (((((masI[0x11][i][j] * masI[6][i][j]) * masI[6][i][j]) * masI[5][i][j]) * masI[9][i][j]) * masI[13][i][j])) + (((((((masI[0x11][i][j] * masI[6][i][j]) * masI[6][i][j]) * masI[9][i][j]) * masI[10][i][j]) * masI[10][i][j]) * masI[12][i][j]) * masI[12][i][j])) - ((((((masI[0x11][i][j] * masI[6][i][j]) * masI[6][i][j]) * masI[9][i][j]) * masI[10][i][j]) * masI[10][i][j]) * masI[12][i][j]);

			compute->det0[i][j] = compute->q0[i][j];
			compute->det1[i][j] = compute->q1[i][j];
			compute->det2[i][j] = compute->q1[i][j] * compute->q2[i][j] - compute->q3[i][j] * compute->q0[i][j];
			compute->det3[i][j] = compute->q1[i][j] * compute->q2[i][j] * compute->q3[i][j] + compute->q0[i][j] * compute->q1[i][j] * compute->q5[i][j] - compute->q1[i][j] * compute->q1[i][j] * compute->q4[i][j] - compute->q0[i][j] * compute->q3[i][j] * compute->q3[i][j];
			compute->det4[i][j] = -compute->q0[i][j] * compute->q0[i][j] * compute->q5[i][j] * compute->q5[i][j] - compute->q0[i][j] * compute->q1[i][j] * compute->q3[i][j] * compute->q6[i][j] + 2 * compute->q0[i][j] * compute->q1[i][j] * compute->q4[i][j] * compute->q5[i][j] + compute->q0[i][j] * compute->q2[i][j] * compute->q3[i][j] * compute->q5[i][j] - compute->q0[i][j] * compute->q3[i][j] * compute->q3[i][j] * compute->q4[i][j] + compute->q1[i][j] * compute->q1[i][j] * compute->q2[i][j] * compute->q6[i][j] - compute->q1[i][j] * compute->q1[i][j] * compute->q4[i][j] * compute->q4[i][j] - compute->q1[i][j] * compute->q2[i][j] * compute->q2[i][j] * compute->q5[i][j] + compute->q1[i][j] * compute->q2[i][j] * compute->q3[i][j] * compute->q4[i][j];
			compute->det5[i][j] = -compute->q0[i][j] * compute->q0[i][j] * compute->q5[i][j] * compute->q5[i][j] * compute->q5[i][j] - 3 * compute->q0[i][j] * compute->q1[i][j] * compute->q3[i][j] * compute->q5[i][j] * compute->q6[i][j] + 2 * compute->q0[i][j] * compute->q1[i][j] * compute->q4[i][j] * compute->q5[i][j] * compute->q5[i][j] + compute->q0[i][j] * compute->q2[i][j] * compute->q3[i][j] * compute->q5[i][j] * compute->q5[i][j] + compute->q0[i][j] * compute->q3[i][j] * compute->q3[i][j] * compute->q3[i][j] * compute->q6[i][j] - compute->q0[i][j] * compute->q3[i][j] * compute->q3[i][j] * compute->q4[i][j] * compute->q5[i][j] - compute->q1[i][j] * compute->q1[i][j] * compute->q1[i][j] * compute->q6[i][j] * compute->q6[i][j] + 2 * compute->q1[i][j] * compute->q1[i][j] * compute->q2[i][j] * compute->q5[i][j] * compute->q6[i][j] + compute->q1[i][j] * compute->q1[i][j] * compute->q3[i][j] * compute->q4[i][j] * compute->q6[i][j] - compute->q1[i][j] * compute->q1[i][j] * compute->q4[i][j] * compute->q4[i][j] * compute->q5[i][j] - compute->q1[i][j] * compute->q2[i][j] * compute->q2[i][j] * compute->q5[i][j] * compute->q5[i][j] - compute->q1[i][j] * compute->q2[i][j] * compute->q3[i][j] * compute->q3[i][j] * compute->q6[i][j] + compute->q1[i][j] * compute->q2[i][j] * compute->q3[i][j] * compute->q4[i][j] * compute->q5[i][j];
			compute->det6[i][j] = compute->q6[i][j] * compute->det5[i][j];

			if (compute->det0[i][j] < 0 || compute->det1[i][j] < 0 || compute->det2[i][j] < 0 || compute->det3[i][j] < 0 || compute->det4[i][j] < 0 || compute->det5[i][j] < 0 || compute->det6[i][j] < 0)
			{
				compute->outPoints[i][j].setColor(false);
			}
			else
			{
				compute->outPoints[i][j].setColor(true);
			}
		}
	}
}

void changeAxes(Compute* compute, int x, int y, double*** masI, vector<double> items)
{
	double minX = items[x] - items[x] / (100 / compute->koeffXMin);
	double maxX = items[x] + items[x] / (100 / compute->koeffXMax);
	double minY = items[y] - items[y] / (100 / compute->koeffYMin);
	double maxY = items[y] + items[y] / (100 / compute->koeffYMax);

	for (int i = 0; i < compute->points; i++)
	{
		for (int j = 0; j < compute->points; j++)
		{
			masI[x][i][j] = minX + j * items[x] / compute->points;
			masI[y][i][j] = minY + i * items[y] / compute->points;
			compute->outPoints[i][j].setCoords(minX + j * items[x] / compute->points, minY + i * items[y] / compute->points);
		}
	}

	formulas(compute, masI);
}

void computing(Compute* compute)
{
	vector<double> items;
	items.push_back(compute->A);
	items.push_back(compute->B);
	items.push_back(compute->C);
	items.push_back(compute->DD);
	items.push_back(compute->omega);
	items.push_back(compute->k);
	items.push_back(compute->V);
	items.push_back(compute->S);
	items.push_back(compute->L);
	items.push_back(compute->c);
	items.push_back(compute->N);
	items.push_back(compute->h);
	items.push_back(compute->sigma);
	items.push_back(compute->a);
	items.push_back(compute->b);
	items.push_back(compute->p);
	items.push_back(compute->alpha);
	items.push_back(compute->betta);

	double*** masI = new double**[18];
	for (int i = 0; i < 18; i++)
	{
		masI[i] = new double*[compute->points];
		for (int j = 0; j < compute->points; j++)
		{
			masI[i][j] = new double[compute->points];
			for (int k = 0; k < compute->points; k++)
			{
				masI[i][j][k] = items.at(i);
			}
		}
	}

	changeAxes(compute, compute->x, compute->y, masI, items);
}

Point** calculate(int points, double A, double B, double C, double DD, double omega, double k, double V, double S, double L, double c, double N,
	double h, double sigma, double a, double b, double p, double alpha, double betta, int x, int y, double koeffXMin, double koeffXMax, double koeffYMin, double koeffYMax)
{
	Compute* compute = new Compute(points);
	compute->A = A;
	compute->B = B;
	compute->C = C;
	compute->DD = DD;
	compute->omega = omega;
	compute->k = k;
	compute->V = V;
	compute->S = S;
	compute->L = L;
	compute->c = c;
	compute->N = N;
	compute->h = h;
	compute->sigma = sigma;
	compute->a = a;
	compute->b = b;
	compute->p = p;
	compute->alpha = alpha;
	compute->betta = betta;
	compute->x = x;
	compute->y = y;
	compute->koeffXMin = koeffXMin;
	compute->koeffXMax = koeffXMax;
	compute->koeffYMin = koeffYMin;
	compute->koeffYMax = koeffYMax;

	computing(compute);
	compute->~Compute();
	
	return compute->outPoints;
}


  // Accepts 1 number from JavaScript, adds 42 and returns the result.
  NAN_METHOD(PassNumber) { 
    // bool value = false;
    // Local<Boolean> retval = Nan::New(!value);

    // Local<Array> squares = New<v8::Array>(3);
    // for (unsigned int i = 0; i < 3; i++ ) {
    //   Nan::Set(squares, i, Nan::New<Boolean>(value));
    // }
    Nan::Maybe<int> points = Nan::To<int>(info[0]); 
    Nan::Maybe<double> A = Nan::To<double>(info[1]); 
    Nan::Maybe<double> B = Nan::To<double>(info[2]); 
    Nan::Maybe<double> C = Nan::To<double>(info[3]); 
    Nan::Maybe<double> DD = Nan::To<double>(info[4]); 
    Nan::Maybe<double> omega = Nan::To<double>(info[5]); 
    Nan::Maybe<double> k = Nan::To<double>(info[6]); 
    Nan::Maybe<double> V = Nan::To<double>(info[7]); 
    Nan::Maybe<double> S = Nan::To<double>(info[8]); 
    Nan::Maybe<double> L = Nan::To<double>(info[9]); 
    Nan::Maybe<double> c = Nan::To<double>(info[10]); 
    Nan::Maybe<double> N = Nan::To<double>(info[11]); 
    Nan::Maybe<double> h = Nan::To<double>(info[12]); 
    Nan::Maybe<double> sigma = Nan::To<double>(info[13]); 
    Nan::Maybe<double> a = Nan::To<double>(info[14]); 
    Nan::Maybe<double> b = Nan::To<double>(info[15]); 
    Nan::Maybe<double> p = Nan::To<double>(info[16]); 
    Nan::Maybe<double> alpha = Nan::To<double>(info[17]); 
    Nan::Maybe<double> betta = Nan::To<double>(info[18]); 
    Nan::Maybe<int> x = Nan::To<int>(info[19]); 
    Nan::Maybe<int> y = Nan::To<int>(info[20]); 
    Nan::Maybe<double> koeffXMin = Nan::To<double>(info[21]); 
    Nan::Maybe<double> koeffXMax = Nan::To<double>(info[22]); 
    Nan::Maybe<double> koeffYMin = Nan::To<double>(info[23]); 
    Nan::Maybe<double> koeffYMax = Nan::To<double>(info[24]); 

    Point** retval = calculate(
        points.FromJust(),
        A.FromJust(), 
        B.FromJust(), 
        C.FromJust(), 
        DD.FromJust(), 
        omega.FromJust(), 
        k.FromJust(), 
        V.FromJust(), 
        S.FromJust(), 
        L.FromJust(), 
        c.FromJust(), 
        N.FromJust(), 
        h.FromJust(), 
        sigma.FromJust(), 
        a.FromJust(), 
        b.FromJust(), 
        p.FromJust(), 
        alpha.FromJust(), 
        betta.FromJust(), 
        x.FromJust(), 
        y.FromJust(), 
        koeffXMin.FromJust(), 
        koeffXMax.FromJust(), 
        koeffYMin.FromJust(), 
        koeffYMax.FromJust()
    );

    Local<Array> result = New<v8::Array>(points.FromJust());
    for (unsigned int i = 0; i < points.FromJust(); i++) {
      Local<Array> tempResult = New<v8::Array>(points.FromJust());
      for (unsigned int j = 0; j < points.FromJust(); j++) {   
		Local<Object> value = Nan::New<Object>();

		Local<String> x = Nan::New<String>("x").ToLocalChecked();
        Local<String> y = Nan::New<String>("y").ToLocalChecked();
        Local<String> color = Nan::New<String>("color").ToLocalChecked();
		
		Nan::Set(value, x, Nan::New<Number>(retval[i][j].x));
		Nan::Set(value, y, Nan::New<Number>(retval[i][j].y));
		Nan::Set(value, color, Nan::New<Boolean>(retval[i][j].color));

        Nan::Set(tempResult, j, value);
      }

      Nan::Set(result, i, tempResult);
    }

    info.GetReturnValue().Set(result); 
  }

  // Called by the NODE_MODULE macro below, 
  // exposes a pass_number method to JavaScript, which maps to PassNumber 
  // above.
  NAN_MODULE_INIT(Init) {
     Nan::Set(target, New<String>("calculate").ToLocalChecked(),
        GetFunction(New<FunctionTemplate>(PassNumber)).ToLocalChecked());
  }

  // macro to load the module when require'd
  NODE_MODULE(my_addon, Init)