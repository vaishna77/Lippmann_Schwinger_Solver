#include <omp.h>
#include <boost/math/special_functions/bessel.hpp>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <bits/stdc++.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <stdlib.h>
#include <Eigen/Dense>
#include <vector>
#include <complex>
#include <set>
#include <algorithm>
#include <filesystem>

const double epsilonRoundOff = pow(10,-8);
// #define EIGEN_DONT_PARALLELIZE
using namespace std;
using namespace Eigen;
const double PI	=	3.1415926535897932384;
int Qchoice;

typedef MatrixXcd Mat;
typedef VectorXcd Vec;
double kappa;
struct pts2D {
	double x,y;
};
//const std::complex<double> I(0.0, 1.0);
typedef std::complex<double> kernel_dtype;

double besselJ(int n, double x) {
	//cout << "n: " << n << "	x: " << x << endl;
	if (n >= 0) {
		double temp = boost::math::cyl_bessel_j(double(n), x);
		return temp;
	}
	else {
		double temp = boost::math::cyl_bessel_j(double(-n), x);
		if (-n%2 == 0)
			return temp;
		else
			return -temp;
	}
}

double besselY(int n, double x) {
	if (n >= 0) {
		double temp = boost::math::cyl_neumann(double(n), x);
		return temp;
	}
	else {
		double temp = boost::math::cyl_neumann(double(-n), x);
		if (-n%2 == 0)
			return temp;
		else
			return -temp;
	}
}

#include "domain2D.hpp"
#include "FarFieldInteraction.hpp"
#include "FarFieldInteraction2.hpp"
#include "NearFieldInteraction1.hpp"
#include "NearFieldInteraction2.hpp"
#include "ACA.hpp"

pts2D SFN_centers[20] = {{-2.5, -2.5},//separatedFineNeighbors
												{-1.5, -2.5},
												{-0.5, -2.5},
												{0.5, -2.5},
												{1.5, -2.5},
												{2.5, -2.5},
												{2.5, -1.5},
												{2.5, -0.5},
												{2.5, 0.5},
												{2.5, 1.5},
												{2.5, 2.5},
												{1.5, 2.5},
												{0.5, 2.5},
												{-0.5, 2.5},
												{-1.5, 2.5},
												{-2.5, 2.5},
												{-2.5, 1.5},
												{-2.5, 0.5},
												{-2.5, -0.5},
												{-2.5, -1.5}
												};
pts2D N_centers[9] = {{-2, -2},//colleagueNeighbors
											 {0, -2},
											 {2, -2},
											 {2, 0},
											 {2, 2},
											 {0, 2},
											 {-2, 2},
											 {-2, 0},
											 {0, 0}};
pts2D FN_centers[12] = {{-1.5, -1.5},//fineNeighbors
												 {-0.5, -1.5},
												 {0.5, -1.5},
												 {1.5, -1.5},
												 {1.5, -0.5},
												 {1.5, 0.5},
												 {1.5, 1.5},
												 {0.5, 1.5},
												 {-0.5, 1.5},
												 {-1.5, 1.5},
												 {-1.5, 0.5},
												 {-1.5, -0.5}};
pts2D CN_centers[12] = {{-3, -3},//coarseNeighbors
												{-1, -3},
												{1, -3},
												{3, -3},
												{3, -1},
												{3, 1},
												{3, 3},
												{1, 3},
												{-1, 3},
												{-3, 3},
												{-3, 1},
												{-3, -1}};

double gaussian(const pts2D r) {
	double R2 = r.x*r.x + r.y*r.y;
	return 1.5*exp(-160.0*R2);
}

double cavity(const pts2D r) {
	double x1 = r.x;
	double x2 = r.y;
	double R = sqrt(x1*x1 + x2*x2);
	double theta = atan2(x2, x1);
	return (1.0-pow(sin(0.5*theta),500))*exp(-2000*(0.1-R*R)*(0.1-R*R));
}

double lens(const pts2D r) {
	double x1 = r.x;
	double x2 = r.y;
	double R = sqrt(x1*x1 + x2*x2);
	return 4*(x2-0.1)*(1.0-erf(25.0*(R-0.3)));
}

double multipleGaussians(const pts2D r) {
	int noOfGaussians = 20;
	double centersX[20] = {-1.17182014706864, //randomly but well-separated gaussians
-0.861066625441127,
-1.14533939124057,
-0.856639124568525,
-0.620066350133170,
-0.457064031294315,
-0.482782006538245,
-0.635070123186999,
0.108594058866164,
-0.0315722590552472,
-0.177274727610547,
0.0342987980936911,
0.369662097773613,
0.591444264226524,
0.401225715100565,
0.666822699314234,
1.10327781188617,
0.827519268441487,
0.873411278211136,
1.09482906183464};
double centersY[20] = {-0.361066625441127,
0.0719520938231629,
0.354660608759432,
1.14336087543147,
-1.12006635013317,
-0.457064031294315,
0.364929876813001,
0.802261227250738,
-0.394084977497864,
-0.0315722590552472,
0.322725272389453,
1.03429879809369,
-1.13033790222639,
0.0137164199764117,
0.401225715100565,
1.16682269931423,
-0.345187707108795,
-0.172480731558513,
0.373411278211136,
1.09482906183464};
	double q = 0.0;
	double a = 0.0013;
	for (size_t i = 0; i < noOfGaussians; i++) {
		double R2 = (r.x-centersX[i])*(r.x-centersX[i]) + (r.y-centersY[i])*(r.y-centersY[i]);
		q += 1.5*exp(-R2/a);
	}
	return q;
}

double flatBump(const pts2D r) {
	double R2 = r.x*r.x + r.y*r.y;
	return 0.5*erfc(5.0*(R2-1.0));
}

double plasma(const pts2D r) {
	double C = 0.4987;
	double psi = 1.0-(r.x-0.15*(1-r.x*r.x))*(r.x-0.15*(1-r.x*r.x)) - C*(1.0+0.3*r.x)*(1.0+0.3*r.x)*r.y*r.y;
	double a[5] = {0.45, 0.195, 0.51, 0.195, 0.63};
	double x[5] = {0.8, 0.54, -0.14, -0.5, 0.18};
	double y[5] = {0.0, -0.28, 0.7, -0.01, 0.8};
	if (psi <= 0.05) {
		return 0.0;
	}
	else {
		double g = 0.0;
		for (size_t i = 0; i < 5; i++) {
			double temp = (r.x-x[i])*(r.x-x[i])+(r.y-y[i])*(r.y-y[i]);
			g += a[i]*exp(-temp/0.01);
		}
		return -1.5*(psi-0.05)-g*cos(0.9*r.y);
	}
}

class userkernel{
public:
	userkernel() {
	};
	double ContrastFunction(const pts2D r) {
		if (Qchoice == 0)
			return gaussian(r);
		else if (Qchoice == 1)
			return multipleGaussians(r);
		else if (Qchoice == 2)
			return plasma(r);
		else if (Qchoice == 3)
			return flatBump(r);
		else if (Qchoice == 4)
			return cavity(r);
		else
			return lens(r);
	};
	std::complex<double> IncidenceFunction(const pts2D r) {
		std::complex<double> q = exp(I*kappa*r.x);
		return q;
	};
	std::complex<double> RHSFunction(const pts2D r) {
		std::complex<double> q = -1.0*kappa*kappa*ContrastFunction(r)*IncidenceFunction(r);
		return q;
	};
	~userkernel() {};
};

double arctan(double y, double x) {//returns atan2 in range (0,2*PI)
	double temp = atan2(y, x);
	if (temp < 0.0)
		return temp + 2*PI;
	else
		return temp;
}

#include "G_FMM2DTree.hpp"
