// Simple interval library from Luc JAULIN, with minor modifications from Fabrice LE BARS and Jeremy NICOLA.

#ifndef __BOX__
#define __BOX__

#include "interval.h"

class box
{
public:
	interval* data;
	int dim;

public:
	//----------------------------------------------------------------------
	// Constructors/destructors
	//----------------------------------------------------------------------
	box();
	box(int);
	box(interval x);
	box(interval x, interval y);
	box(interval x, interval y, interval z);
	box(interval x, int n);
	box(const box&);
	~box();
	//----------------------------------------------------------------------
	// Operators
	//----------------------------------------------------------------------
	box& operator=(const box&);
	friend box operator+(const box&, const box&);
	friend box operator-(const box&);
	friend box operator-(const box&, const box&);
	friend box operator*(const interval&, const box&);
	friend box operator*(const box&, const interval&);
	friend box operator*(const double, const box&);
	friend box operator*(const box&, const double);
	friend box operator&(const box&, const box&);
	friend box operator|(const box&, const box&);
	friend bool operator==(const box&, const box&);
	interval& operator[](int) const;
	friend std::ostream& operator<<(std::ostream& os, const box& X);
#ifdef QT_VERSION 
	friend QDebug operator<<(QDebug os, const box&X)
	{   
		os.nospace() << "box :" << "\t dim=" << X.dim << "\n";
		if (X.IsEmpty()) os.nospace() << "EmptyBox";
		for (int i = 1; i <= X.dim; i++)
			os.nospace() << "  " << i << ": "<< X[i] << "\n";
		return (os.space());
	}
#endif // QT_VERSION 
	//----------------------------------------------------------------------
	// Member functions
	//----------------------------------------------------------------------
	box& Intersect(const box& Y);
	double Width(void);
	double SumWidth(void);
	bool IsEmpty(void) const;
	void Resize(int);
};

//----------------------------------------------------------------------
// Box-valued functions
//----------------------------------------------------------------------
box Inf(box);
box Sup(box);
box Center(const box&);
box Center(const box&, std::vector<int>&);
box Zeros(int);
box Empty(int);
box EmptyBox(int);
box EmptyBox(const box&);
box Infinity(int);
box Rand(const box& X);
box Concat(const box&, const box&);
box Proj(const box&, int, int);
box Inter(const box&, const box&);
box Inter(std::vector<box>&);
box Union(const box&, const box&);
box Union(std::vector<box>&);
box Inflate(const box&, double);
//----------------------------------------------------------------------
// Other functions
//----------------------------------------------------------------------
int Size(const box&);
double Width(const box&);
double Width(const box&, std::vector<int>&);
double Width(const box&, const box&);
double Volume(const box&);
double Marge(const box&, const box&);
bool IsBox(const box&);
void Update(const box&);
double Angle(const box&, const box&); // Il faut des vecteurs de dim 2
interval Norm(const box&);
interval NormEuclid(const box&, const box&);
interval NormInf(const box&, const box&);
#define ProduitScalaire Scal
interval Scal(const box&, const box&);
#define Determinant Det
interval Det(const box&, const box&);
double Eloignement(const box&, const box&);
double Eloignement2(const box&, const box&);
double EloignementRelatif2(const box&, const box&);
bool Disjoint(const box&, const box&);
bool Subset(const box&, const box&);
bool Subset(const box& X, const box& Y, double epsilon);
bool SubsetStrict(const box&, const box&);
iboolean In(const box&, const box&);
bool Prop(const box&, const box&);
int AxePrincipal(const box&);
int AxePrincipal(const box&, const box&);
int AxePrincipal(const box&, std::vector<int>&);
double decrease(const box&, const box&);
double decrease(const box&, const box&, std::vector<int>);
//----------------------------------------------------------------------
// Contractors
//----------------------------------------------------------------------
void Cadd(box& Z, box& X, box& Y, int dir = 0);
void Csub(box& Z, box& X, box& Y, int dir = 0);
#define CProd Cmul
void Cmul(box& Y, interval& a, box& X, int dir = 0);
void Cnorm(interval& R, box& X);
#define Cdistance Cdist
void Cdist(interval& R, box& X, box& Y);
#define CProdScalaire Cscal
void Cscal(interval& R, box& X, box& Y);
#define COrtho Cortho
void Cortho(box& X, box& Y);
void Cnotin(box& X, const box& Y);
void C_q_in(box&, int, std::vector<box>&);
//----------------------------------------------------------------------
// Other
//----------------------------------------------------------------------
void Bisect(box& X, box& X1, box& X2);
void Bisect(box& X, box& X1, box& X2, box& V);
void Bisect(box& X, box& X1, box& X2, std::vector<int>& v);
void BisectAlong(box& X, box& X1, box& X2, int i);
void BisectHere(box& X, box& X1, box& X2, int i, double here);
void Trisect(box&, box&, box&, box&);
void TrisectAlong(box&, box&, box&, box&, int);
void Decoup(box&, box&, box&);
void DecoupAlong(box&, box&, box&, int);
//void CheckRange(box&,box&);
void Sucre(box&, box&);
// Operation sur les boites
std::vector<box>* diff(box x, box y);

#endif // __BOX__
