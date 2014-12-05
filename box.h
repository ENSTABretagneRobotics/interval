// Simple interval library from Luc JAULIN, with minor modifications from Fabrice LE BARS and Jeremy NICOLA.

#ifndef __BOX__
#define __BOX__

#include "interval.h"

using namespace std;

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
	//----------------------------------------------------------------------
#ifdef QT_VERSION 
	inline friend QDebug operator<<(QDebug os, const box&X)
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
box Center(const box&, vector<int>&);
box Zeros(int);
box Empty(int);
box EmptyBox(int);
box EmptyBox(const box&);
box Infinity(int);
box Rand(const box& X);
box Concat(const box&, const box&);
box Proj(const box&, int, int);
box Inter(const box&, const box&);
box Inter(vector<box>&);
box Union(const box&, const box&);
box Union(vector<box>&);
box Inflate(box&, double);
//----------------------------------------------------------------------
// Other functions
//----------------------------------------------------------------------
double Width(box&);
double Width(box&, vector<int>&);
double Width(box&, box&);
int Size(const box&);
double Volume(box&);
double Marge(box, box);
bool IsBox(box);
void Update(box&);
double Angle(box&, box&); // Il faut des vecteurs de dim 2
interval Norm(box);
interval NormEuclid(box, box);
interval NormInf(box, box);
interval ProduitScalaire(box&, box&);
#define Determinant Det
interval Det(box&, box&);
double Eloignement(box&, box&);
double Eloignement2(box&, box&);
double EloignementRelatif2(box&, box&);
bool Disjoint(const box&, const box&);
bool Subset(box&, box&);
bool Subset(box& X, box& Y, double epsilon);
bool SubsetStrict(box&, box&);
iboolean In(box, box);
//void Inter1(box&,box&,const box&,const box&,const box&);
bool Prop(box&, box&);
int AxePrincipal(box&);
int AxePrincipal(box&, box&);
int AxePrincipal(box&, vector<int>&);
double decrease(const box&, const box&);
double decrease(const box&, const box&, vector<int>);
//----------------------------------------------------------------------
// Contractors
//----------------------------------------------------------------------
void Cadd(box& Z, box& X, box& Y, int sens = 0);
void Csub(box& Z, box& X, box& Y, int sens = 0);
#define CProd Cmul
void Cmul(box& Y, interval& a, box& X, int sens = 0);
void Cnorm(interval& R, box& X);
void Cdistance(interval& R, box& X, box& Y);
#define CProdScalaire Cscal
void Cscal(interval& R, box& X, box& Y);
#define COrtho Cortho
void Cortho(box& X, box& Y);
void Cnotin(box& X, const box& Y);
void C_q_in(box&, int, vector<box>&);
//----------------------------------------------------------------------
// Other
//----------------------------------------------------------------------
void Bisect(box& X, box& X1, box& X2);
void Bisect(box& X, box& X1, box& X2, box& V);
void Bisect(box& X, box& X1, box& X2, vector<int>& v);
void BisectAlong(box& X, box& X1, box& X2, int i);
void BisectHere(box& X, box& X1, box& X2, int i, double here);
void Trisect(box&, box&, box&, box&);
void TrisectAlong(box&, box&, box&, box&, int);
void Decoup(box&, box&, box&);
void DecoupAlong(box&, box&, box&, int);
//void CheckRange(box&,box&);
void Sucre(box&, box&);
// Operation sur les boites
vector<box>* diff(box x, box y);

#endif // __BOX__
