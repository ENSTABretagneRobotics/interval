// Simple interval library from Luc JAULIN, with minor modifications from Fabrice LE BARS.

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
	// Constructors
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
		os.nospace() << "box :" << "\t dim=" << Size(X) << "\n";
		if (X.IsEmpty()) os.nospace() << "EmptyBox";
		for (int i = 1; i <= Size(X); i++)
			os.nospace() << "  " << i << ": "<< X[i] << "\n";
		return (os.space());
	}
#endif // QT_VERSION 
	//----------------------------------------------------------------------
	// Box-valued functions
	//----------------------------------------------------------------------
	friend box Inf(box);
	friend box Sup(box);
	friend box Center(const box&);
	friend box Center(const box&, vector<int>&);
	friend box Zeros(int);
	friend box Empty(int);
	friend box EmptyBox(int);
	friend box EmptyBox(const box&);
	friend box Infinity(int);
	friend box Rand(const box& X);
	friend box Concat(const box&, const box&);
	friend box Proj(const box&, int, int);
	friend box Inter(const box&, const box&);
	friend box Inter(vector<box>&);
	friend box Union(const box&, const box&);
	friend box Union(vector<box>&);
	box& Intersect(const box& Y);
	friend box Inflate(box&, double);
	//----------------------------------------------------------------------
	// Other functions
	//----------------------------------------------------------------------
	double Width(void);
	friend double Width(box&);
	friend double Width(box&, vector<int>&);
	friend double Width(box&, box&);
	friend int Size(const box&);
	friend double Volume(box&);
	friend double Marge(box, box);
	double SumWidth(void);
	bool IsEmpty(void) const;
	friend bool IsBox(box);
	friend void Update(box&);
	void Resize(int);
	friend double Angle(box&, box&); // Il faut des vecteurs de dim 2
	friend interval Norm(box);
	friend interval NormEuclid(box, box);
	friend interval NormInf(box, box);
	friend interval ProduitScalaire(box&, box&);
#define Determinant Det
	friend interval Det(box&, box&);
	friend double Eloignement(box&, box&);
	friend double Eloignement2(box&, box&);
	friend double EloignementRelatif2(box&, box&);
	friend bool Disjoint(const box&, const box&);
	friend bool Subset(box&, box&);
	friend bool Subset(box& X, box& Y, double epsilon);
	friend bool SubsetStrict(box&, box&);
	friend iboolean In(box, box);
	//friend void Inter1(box&,box&,const box&,const box&,const box&);
	friend bool Prop(box&, box&);
	friend int AxePrincipal(box&);
	friend int AxePrincipal(box&, box&);
	friend int AxePrincipal(box&, vector<int>&);
	friend double decrease(const box&, const box&);
	friend double decrease(const box&, const box&, vector<int>);
	//----------------------------------------------------------------------
	// Contractors
	//----------------------------------------------------------------------
	friend void Cadd(box&, box&, box&, int sens = 0);
	friend void Csub(box&, box&, box&, int sens = 0);
#define CProd Cmul
	friend void Cmul(box& Y, interval& a, box& X, int sens = 0);
	friend void Cnorm(interval& R, box& X);
	friend void Cdistance(interval& R, box& X, box& Y);
#define CProdScalaire Cscal
	friend void Cscal(interval& R, box& X, box& Y);
	friend void COrtho(box& X, box& Y);
	friend void Cnotin(box& X, const box& Y);
	friend void C_q_in(box&, int, vector<box>&);
	//----------------------------------------------------------------------
	// Other
	//----------------------------------------------------------------------
	friend void Bisect(box&, box&, box&);
	friend void Bisect(box&, box&, box&, box&);
	friend void Bisect(box&, box&, box&, vector<int>&);
	friend void BisectAlong(box&, box&, box&, int);
	friend void BisectHere(box&, box&, box&, int, double);
	friend void Trisect(box&, box&, box&, box&);
	friend void TrisectAlong(box&, box&, box&, box&, int);
	friend void Decoup(box&, box&, box&);
	friend void DecoupAlong(box&, box&, box&, int);
	//friend void CheckRange(box&,box&);
	friend void Sucre(box&, box&);
	// Operation sur les boites
	friend vector<box>* diff(box x, box y);
};

//box Empty(int);

#endif // __BOX__
