// Simple interval library from Luc JAULIN, with minor modifications from Fabrice LE BARS and Jeremy NICOLA.

#ifndef __IMATRIX__
#define __IMATRIX__

#include "rmatrix.h"

class imatrix
{
public:
	TNT::Array2D <interval> data;

public:
	//----------------------------------------------------------------------
	// Constructors/destructors
	//----------------------------------------------------------------------
	imatrix();
	imatrix(int, int);
	imatrix(const rmatrix&);
	imatrix(const box&);
	imatrix(const imatrix&);
	~imatrix();
	//----------------------------------------------------------------------
	// Operators
	//----------------------------------------------------------------------
	imatrix& operator=(const imatrix&);
	friend imatrix operator+(const imatrix&, const imatrix&);
	friend imatrix operator-(const imatrix&);
	friend imatrix operator-(const imatrix&, const imatrix&);
	friend imatrix operator*(const imatrix&, const imatrix&);
	friend imatrix operator*(const rmatrix& X, const imatrix& Y);
	friend imatrix operator*(const double a, const imatrix& X);
	friend box operator*(const imatrix& A, const box& x);
	//interval& imatrix::operator()(int i, int j)  const;
	interval operator()(int i, int j) const;
	friend std::ostream& operator<<(std::ostream&, const imatrix&);
	//----------------------------------------------------------------------
	// Member functions
	//----------------------------------------------------------------------
	interval GetVal(int, int) const;
	void SetVal(int, int, interval);
	int dim1(void) const;
	int dim2(void) const;
};

//----------------------------------------------------------------------
// imatrix-valued functions
//----------------------------------------------------------------------
imatrix iZeros(int n, int m);
imatrix iEye(int n);
imatrix Transpose(const imatrix& X);
imatrix RotationPhiThetaPsi(interval& phi, interval& theta, interval& psi);
//----------------------------------------------------------------------
// Other functions
//----------------------------------------------------------------------
rmatrix Center(const imatrix&);
box Row(const imatrix& B, int i);
box Column(const imatrix& B, int j);
box ToBox(const imatrix&);
//----------------------------------------------------------------------
// Contractors
//----------------------------------------------------------------------
#define Cmult Cmul
void Cmul(imatrix& C, imatrix& A, imatrix& B);
void Cmul(box& c, imatrix& A, box& b);
void Crot(imatrix&);
void Cantisym(imatrix& A);

/*void         Update            (imatrix&);
imatrix          Rand              (const imatrix& X);
interval     Determinant       (imatrix&, imatrix&);
bool         Emptyimatrix          (const imatrix&);
bool	    Disjoint	      (const imatrix&,const imatrix&);
iboolean	    In		      (imatrix,imatrix);
imatrix	    Inf 	      (imatrix);
imatrix	    Inter 	      (const imatrix&,const imatrix&);
imatrix          Concat            (const imatrix&, const imatrix&);
interval     Norm              (imatrix);
interval     NormEuclid        (imatrix, imatrix);
interval     NormInf           (imatrix, imatrix);
bool	    Subset	      (imatrix&,imatrix&);
bool	    Subset	      (imatrix&,imatrix&,double);
imatrix	    Sup 	      (imatrix);
imatrix	    Union	      (imatrix&,imatrix&);
double	    Width	      (imatrix&);
double	    Width	      (imatrix&,vector<int>&);
double	    Width	      (imatrix&,imatrix&);
imatrix          Zeros             (int,int);*/

#endif // __IMATRIX__
