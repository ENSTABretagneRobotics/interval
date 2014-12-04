// Simple interval library from Luc JAULIN, with minor modifications from Fabrice LE BARS.

#ifndef __INTERVAL__
#define __INTERVAL__

#ifdef _MSC_VER
// Disable some Visual Studio warnings.
#	ifndef CRT_SECURE_NO_DEPRECATE
#		define CRT_SECURE_NO_DEPRECATE
#	endif // CRT_SECURE_NO_DEPRECATE
#	ifndef _CRT_SECURE_NO_WARNINGS
#		define _CRT_SECURE_NO_WARNINGS
#	endif // _CRT_SECURE_NO_WARNINGS
//#	ifndef _CRT_NONSTDC_NO_WARNINGS
//#		define _CRT_NONSTDC_NO_WARNINGS
//#	endif // _CRT_NONSTDC_NO_WARNINGS
#endif // _MSC_VER

#include <vector>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include "iboolean.h"

#ifdef _MSC_VER
#ifndef UNREFERENCED_PARAMETER
#define UNREFERENCED_PARAMETER(P) (P)
#endif // UNREFERENCED_PARAMETER
#endif // _MSC_VER

#ifdef __GNUC__
#undef UNREFERENCED_PARAMETER
#define UNREFERENCED_PARAMETER(P) (void)(P)
#endif // __GNUC__

#ifdef __BORLANDC__
#undef UNREFERENCED_PARAMETER
#define UNREFERENCED_PARAMETER(P) 
#endif // __BORLANDC__

// To avoid Visual Studio 2013 warning about overflow in floating-point constant arithmetic 
// each time INFINITY or NAN is used.
#if (_MSC_VER >= 1800)
#pragma warning(disable : 4056)
#endif // (_MSC_VER >= 1800)

#ifndef INFINITY
#define INFINITY HUGE_VAL
#endif // INFINITY

#if defined(_MSC_VER) || defined(__BORLANDC__) 
// Used to define NAN (Not A Number).
#ifndef NAN
extern const unsigned long nan[2];
extern const double nan_double;
#define NAN nan_double
#define NAN_CONSTS_NEEDED
#endif // NAN
#endif // defined(_MSC_VER) || defined(__BORLANDC__) 

// Infinity is denoted by oo.
#ifndef oo
//#define oo 1.0/0.0
//#define oo 1000000000.0
#define oo INFINITY
#endif // oo

// Declaration of the interval class.
class interval;

// Used to define NAI (Not An Interval).
#ifndef NAI
extern const interval nai;
#define NAI nai
#define NAI_CONST_NEEDED
#endif // NAI

using namespace std;

// Include <QDataStream> and <QDebug> before this file to be able to use Qt specific features if you have Qt.
#ifdef QT_VERSION 
class QDataStream;
class QDebug;
#else
#define qDebug() std::cout
#endif // QT_VERSION 

// Deprecated.
typedef double reel;

//----------------------------------------------------------------------
// Useful real-valued functions
//----------------------------------------------------------------------
double Min(vector<double>& x);
double Max(vector<double>& x);
double Sign(const double x);
double Chi(const double a, const double b, const double c);
double Arccossin(const double x, const double y);
double Arg(const double x, const double y);
double Det(double ux, double uy, double vx, double vy);
double DistanceDirSegment(double mx, double my, double theta, double ax, double ay, double bx, double by);
void DistanceDirSegment(double& d, double& phi, double mx, double my, double theta, double ax, double ay, double bx, double by);
double DistanceDirSegments(double mx, double my, double theta, 
						   vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by);
void DistanceDirSegments(double& d, double& phi, double mx, double my, double theta, 
						 vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by);
#define DistanceDirCercle DistanceDirCircle
#define DistanceDirCercles DistanceDirCircles
double DistanceDirCircle(double mx, double my, double theta, double cx, double cy, double r);
void DistanceDirCircle(double& d, double& phi, double mx, double my, double theta, double cx, double cy, double r);
double DistanceDirCircles(double mx, double my, double theta, vector<double> cx, vector<double> cy, vector<double> r);
void DistanceDirCircles(double& d, double& phi, double mx, double my, double theta, 
						vector<double> cx, vector<double> cy, vector<double> r);
double DistanceDirSegmentsOrCircles(double mx, double my, double theta,
									vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by,
									vector<double> cx, vector<double> cy, vector<double> r);
void DistanceDirSegmentsOrCircles(double& d, double& phi, double mx, double my, double theta,
								  vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by,
								  vector<double> cx, vector<double> cy, vector<double> r);

// Used by q_in only.
class borne
{
public:
	double val;
	int ouverture;
	borne();
	borne(const double&, const int&);
	friend bool operator<(const borne &x, const borne &y);
};

class interval
{
public:
	double inf;
	double sup;
	bool isEmpty;

public:
	//----------------------------------------------------------------------
	// Constructors
	//----------------------------------------------------------------------
	interval();
	interval(const double&); // Pas const pour conversion double -> interval
	interval(const double&, const double&);
	interval(const interval&);
	//----------------------------------------------------------------------
	// Operators
	//----------------------------------------------------------------------
	interval& operator=(const interval&);
	friend interval operator+(const interval&, const interval&);
	friend interval operator-(const interval&);
	friend interval operator-(const interval&, const interval&);
	friend interval operator*(const double, const interval&);
	friend interval operator*(const interval&, const double);
	friend interval operator*(const interval&, const interval&);
	friend interval operator/(const interval&, const interval&);
	friend interval operator&(const interval&, const interval&);
	friend interval operator|(const interval&, const interval&);
	friend bool operator==(const interval&, const interval&);
	friend std::ostream& operator<<(std::ostream& os, const interval& a);
#ifdef QT_VERSION 
	inline friend QDataStream& operator<<(QDataStream& s, const interval& i)
	{
		s << i.inf << i.sup << i.isEmpty;
		return s;
	}

	inline friend QDataStream& operator>>(QDataStream& s, interval& i)
	{
		s >> i.inf >> i.sup >> i.isEmpty;
		return s;
	}

	inline friend QDebug operator<<(QDebug os, const interval& a)
	{
		if (a.isEmpty) os.nospace() << "EmptyInterval";
		else if (a.inf != a.sup)
		{ 
			os.nospace() << "[" << a.inf << ", " << a.sup << "] ";
		}
		else os.nospace() << a.inf;
		return os.space();
	}
#endif // QT_VERSION 
	//----------------------------------------------------------------------
	// Interval-valued functions
	//----------------------------------------------------------------------
	friend interval Min(const interval&, const interval&);
	friend interval Min(const interval&, const interval&, const interval&);
	friend interval Max(const interval&, const interval&);
	friend interval Max(const interval&, const interval&, const interval&);
	//Sign?
	friend interval Abs(const interval&);
	friend interval Modulo(const interval& a, double x);
	friend interval Sqr(const interval&);
	friend interval Sqrt(const interval&);
	friend interval InvSqrt(const interval&);
	friend interval Exp(const interval&);
	friend interval Log(const interval&);
	friend interval Pow(const interval&, int);
	friend interval Pow(const interval& x, int num, int den);
	friend interval Power(const interval&, int);
	friend interval PowRoot(const interval& x, int num, int den);
	friend interval Cos(const interval&);
	friend interval Sin(const interval&);
	friend interval Tan(const interval&);
	//Arg/Atan2, param order like double version?
	friend interval Det(interval& ux, interval& uy, interval& vx, interval& vy);
	friend interval Det(interval& ux, interval& uy, double& vx, double& vy);
	friend interval Step(const interval&);
	friend interval Parabole(const interval&, double, double, double);
	friend interval Inter(const interval&, const interval&);
	friend interval Inter(vector<interval> x);
	friend interval Union(const interval&, const interval&);
	friend interval Union(vector<interval> x);
	interval& Intersect(const interval&);
	friend interval InterMin(const interval&, const interval&, char);
	friend interval Inflate(const interval&, double);
#define Enveloppe Envelope
	friend interval Envelope(vector<double>& x);
	//----------------------------------------------------------------------
	// Other functions
	//----------------------------------------------------------------------
	friend double Inf(const interval&);
	friend double Sup(const interval&);
	friend double Center(const interval&);
	friend double Width(const interval&);
	friend double Marge(const interval&, const interval&);
#define ToReel ToReal
#define Todouble ToReal
	friend double ToReal(const interval&);
	friend double Rand(const interval&);
	friend double Eloignement(const interval&, const interval&);
	friend double AbsMax(const interval&);
	friend bool OverLap(const interval&, const interval&);
	friend bool Disjoint(const interval&, const interval&);
	friend bool Subset(const interval&, const interval&);
	friend bool Subset(const interval&, const interval&, double epsilon);
	friend bool SubsetStrict(const interval& a, const interval& b);
	friend iboolean In(const interval&, const interval&);
	friend bool In(double, const interval&);
	//----------------------------------------------------------------------
	// Contractors
	//----------------------------------------------------------------------
#define Cplus Cadd
	friend void Cadd(interval& Z, interval& Y, interval& X, int sens = 0);
	friend void Cadd(interval& Z, double y, interval& X, int sens = 0);
	friend void Cadd(interval& Z, interval& Y, double x, int sens = 0);
	friend void Cadd(double z, interval& Y, interval& X, int sens = 0);
#define Cmoins Csub
	friend void Csub(interval& Z, interval& Y, interval& X, int sens = 0);
	friend void Csub(interval& Z, double y, interval& X, int sens = 0);
	friend void Csub(interval& Z, interval& Y, double x, int sens = 0);
	friend void Csub(double z, interval& Y, interval& X, int sens = 0);
	friend void Csub(interval& Y, interval& X, int sens = 0);
#define Cprod Cmul
	friend void Cmul(interval& Z, interval& Y, interval& X, int sens = 0);
	friend void Cmul(interval& Z, double y, interval& X, int sens = 0);
	friend void Cmul(interval& Z, interval& Y, double x, int sens = 0);
	friend void Cdiv(interval& Z, interval& Y, interval& X, int sens = 0);
#define Cegal Cequal
	friend void Cequal(interval& Y, interval& X, int sens);
	friend void Cequal(interval& Y, interval& X);
	friend void Cmin(interval& a, interval& b, interval& c, int sens = 0);
	friend void Cmin(interval& a, interval& b, interval& c, interval& d, int sens = 0);
	friend void Cmin(interval& a, interval& b, interval& c, interval& d, interval& e, int sens = 0);
	friend int Cmin(interval& a, vector<interval>& x, int sens = 0);
	friend void Cmax(interval& a, interval& b, interval& c, int sens = 0);
	friend void Cabs(interval& Y, interval& X, int sens = 0);
#define Csame_sign Csign
	friend void Csign(interval& Y, interval& X);
	friend void Csign(interval& Y, interval& X, int sens, double a = 0);
	friend void Cchi(interval& F, interval& A, interval& B, interval& C);
	friend void Cgeq(interval& Y, interval& X);
	friend void Cinteger(interval&);
	friend void Cboolean(interval&);
	friend void Csqr(interval& Y, interval& X, int sens = 0);
	friend void Csqrt(interval& Y, interval& X, int sens = 0);
	friend void Cexp(interval& Y, interval& X, int sens = 0);
	friend void Clog(interval& Y, interval& X, int sens = 0);
#define Cpower Cpow
	friend void Cpow(interval& Y, interval& X, int n, int sens = 0);
	friend void Ccos(interval& Y, interval& X, int sens = 0);
	friend void Csin(interval& Y, interval& X, int sens = 0);
	friend void Ctan(interval& Y, interval& X, int sens = 0);
	friend void Catan(interval& Y, interval& X, int sens = 0);
	friend void Csinc(interval& Y, interval& X, int sens = 0);
	//Carg (different from CAngle, with less parameters...)?
	friend int CAngle(interval& X2, interval& Y2, interval& Theta, interval& X1, interval& Y1, bool StrongAngle); // Deprecated.
#define CNorm Cnorm
	friend void Cnorm(interval& N, interval& X, interval& Y);
	friend void Cnorm(interval& N, interval& X, interval& Y, interval& Z, int sens = 0);
#define CScal Cscal
	friend void Cscal(interval& s, interval& ux, interval& uy, interval& vx, interval& vy);
	friend void Cscal(interval& s, double& ux, double& uy, interval& vx, interval& vy);
#define CDet Cdet
	friend void Cdet(interval& det, interval& ux, interval& uy, interval& vx, interval& vy, int sens = 0);
	friend void Cdet(interval& det, double& ux, double& uy, interval& vx, interval& vy, int sens = 0);
	friend void Cdet(interval& det, interval& ux, interval& uy, double& vx, double& vy, int sens = 0);
	friend void Cstep(interval& Y, interval& X);
	friend void Cstep(interval& Y, interval& X, int sens, double a = 0);
	friend void Cramp(interval& Y, interval& X, int sens = 0, double a = 0);
	friend void Cheaviside(interval& Y, interval& X, int sens = 0, double a = 0);
	friend void Crect(interval& Z, interval& X, interval& Y, int sens = 0);
	friend void Crect(interval& Y, interval& X, int sens = 0);
	friend void Ctriangle(interval& Y, interval& X, int sens = 0);
	friend void CDistanceDirLine(interval& dist, interval& mx, interval& my, interval& theta, 
		double& ax, double& ay, double& bx, double& by);
	friend int CDistanceDirSegment(interval& dist, interval& mx, interval& my, interval& theta, 
		double ax, double ay, double bx, double by, int sens = 0);
	friend void CDistanceDirSegments(interval& distmin, interval& mx, interval& my, interval& theta, 
		vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by);
	friend void CPointInLine(interval& mx, interval& my, double& ax, double& ay, double& bx, double& by);
#define CinSegment CPointInSegment
	friend void CPointInSegment(interval& mx, interval& my, double ax, double ay, double bx, double by);
#define CinSegments CPointInSegments
	friend void CPointInSegments(interval& mx, interval& my, vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by);
#define CinCircle CPointInCircle
	friend void CPointInCircle(interval& mx, interval& my, double cx, double cy, double r);
#define CinCircles CPointInCircles
	friend void CPointInCircles(interval& mx, interval& my, vector<double> cx, vector<double> cy, vector<double> r, bool truth = true);
#define CinSegmentsOrCircles CPointInSegmentsOrCircles
	friend void CPointInSegmentsOrCircles(interval& mx, interval& my, 
		vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by, 
		vector<double> cx, vector<double> cy, vector<double> r);
	friend void CPointOutsideSegment(interval& mx, interval& my, double& ax, double& ay, double& bx, double& by, bool outer);
	friend void CPointOutsideSegments(interval& mx, interval& my, 
		vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by, bool outer);

	friend void CPoseInSegment(interval& mx, interval& my, interval& phi, double& ax, double& ay, double& bx, double& by);
	friend void CPoseInSegments(interval& mx, interval& my, interval& phi, 
		vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by);
	friend void CPoseInCircle(interval& mx, interval& my, interval& phi, double& cx, double& cy, double& r);
	friend void CPoseInCircles(interval& mx, interval& my, interval& phi, vector<double> cx, vector<double> cy, vector<double> r);
	friend void CPoseInSegmentsOrCircles(interval& mx, interval& my, interval& malpha, 
		vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by, 
		vector<double> cx, vector<double> cy, vector<double> r);

	friend void CPoseTrans(interval& qx, interval& qy, interval& d, interval& px, interval& py, interval& theta);  //Go straight
	friend void CPoseRotTrans(interval& qx, interval& qy, interval& beta, interval& phi, interval& d, interval& px, interval& py, interval& alpha);
	friend void CPoseTransInWallsOrCircles(interval& px, interval& py, interval& alpha, interval& d, 
		vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by, 
		vector<double> cx,vector<double> cy, vector<double> r);
	friend void CPoseTransRotInWallsOrCircles(interval& px, interval& py, interval& alpha, interval& d, interval& psi, 
		vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by, 
		vector<double> cx, vector<double> cy, vector<double> r);
	friend void CPoseRotTransRotInWallsOrCircles(interval& px, interval& py, interval& alpha, interval& phi, interval& d, interval& psi, 
		vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by, 
		vector<double> cx, vector<double> cy, vector<double> r);
	friend void CPoseRotTransPointInWallsOrCircles(interval& px, interval& py, interval& alpha, interval& phi, interval& d, 
		vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by, 
		vector<double> cx, vector<double> cy, vector<double> r);
	friend void CPoseTransPointInWall(interval& px,interval& py, interval& alpha, interval& d, 
		double ax, double ay, double bx, double by, bool truth = true);
	friend void CPoseTransPointInWalls(interval& px,interval& py, interval& alpha, interval& d0, 
		vector<double>& ax, vector<double>& ay, vector<double>& bx, vector<double>& by, bool truth = true);
	friend void CPoseTransPointInWallsOrCircles(interval& px,interval& py, interval& alpha, interval& d0, 
		vector<double> ax,vector<double> ay,vector<double> bx,vector<double> by, 
		vector<double> cx, vector<double> cy, vector<double> r, bool truth = true);
	friend void CPoseTowardSegment(interval& mx, interval& my, interval& theta, 
		double& ax, double& ay, double& bx, double& by, bool truth = true);

#define Ccroisepas Cnocross
	friend void Cnocross(interval& px, interval& py, interval& mx, interval& my, double& ax, double& ay, double& bx, double& by);
#define CPatteCroiseAucunSegment CLegCrossNoSegment
	friend void CLegCrossNoSegment(interval& dist, interval& px, interval& py, interval& theta, 
		vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by);
	friend void CLegOnWalls(interval& dist, interval& px, interval& py, interval& theta, 
		vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by);
	friend void CLegOnWallsOrCircles(interval& dist, interval& px, interval& py, interval& theta, 
		vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by, 
		vector<double> cx, vector<double> cy, vector<double> r);

	//------- Procedure de reduction elementaires sur les intervalles ----------
	friend void Contract0(char, interval&, interval&, int);
	//friend void Contract0 (char, interval&, interval&, int, int);
	//friend void Contract0 (char, interval&, interval&, int, int n=0);
	friend void Contract0(char, interval&, interval&, interval&, int);
	//friend void Contract0 (char, interval&, double&, interval&, int); //Luc
	friend void Contract0(char, interval&);
	//modifs???
	friend void ShowContraction(interval&, interval&, interval&, interval&);
	friend void IntButterfly(interval& Y, interval Yo, interval dY, interval& X, interval Xo, int sens);
	//modifs???
	friend void Inter1(interval&, interval&, const interval&, const interval&, const interval&);
	friend void Sucre(interval&, const interval&);

	friend void Cnotin(interval& X, interval& Y);
	friend void C_q_in(interval& x, int q, vector<interval>& y);
	//----------------------------------------------------------------------
	// Other
	//----------------------------------------------------------------------
	friend void diffI(interval &x0, interval &x1, interval &c0, interval &c1);
	// Primitive inclusion tests
	friend iboolean TestDiskExists(const interval& X, const interval& Y, const interval& P1, const interval& P2, const interval& P3);
	friend iboolean TestDiskForall(const interval& X, const interval& Y, const interval& P1, const interval& P2, const interval& P3);
};

#endif // __INTERVAL__
