// Simple interval library from Luc JAULIN, with minor modifications from Fabrice LE BARS.

#ifndef __IMATRIX__
#define __IMATRIX__

#include "rmatrix.h"

class imatrix
{
public:
	TNT::Array2D <interval> data;

public:
	//----------------------------------------------------------------------
	// Constructors
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
	// imatrix-valued functions
	//----------------------------------------------------------------------
	friend imatrix iZeros(int n, int m);
	friend imatrix iEye(int n);
	friend imatrix Transpose(const imatrix& X);
	friend imatrix RotationPhiThetaPsi(interval& phi, interval& theta, interval& psi);
	//----------------------------------------------------------------------
	// Other functions
	//----------------------------------------------------------------------
	interval GetVal(int,int) const;
	void SetVal(int,int,interval);
	int dim1(void) const;
	int dim2(void) const;
	friend rmatrix Center(const imatrix&);
	friend box Row(const imatrix& B, int i);
	friend box Column(const imatrix& B, int j);
	friend box ToBox(const imatrix&);
	//----------------------------------------------------------------------
	// Contractors
	//----------------------------------------------------------------------
#define Cmult Cmul
	friend void Cmul(imatrix& C, imatrix& A, imatrix& B);
	friend void Cmul(box& c, imatrix& A, box& b);
	friend void Crot(imatrix&);
	friend void Cantisym(imatrix& A);

	/*friend void         Update            (imatrix&);
	friend imatrix          Rand              (const imatrix& X);
	friend interval     Determinant       (imatrix&, imatrix&);
	friend bool         Emptyimatrix          (const imatrix&);
	friend bool	    Disjoint	      (const imatrix&,const imatrix&);
	friend iboolean	    In		      (imatrix,imatrix);
	friend imatrix	    Inf 	      (imatrix);
	friend imatrix	    Inter 	      (const imatrix&,const imatrix&);
	friend imatrix          Concat            (const imatrix&, const imatrix&);
	friend interval     Norm              (imatrix);
	friend interval     NormEuclid        (imatrix, imatrix);
	friend interval     NormInf           (imatrix, imatrix);
	friend bool	    Subset	      (imatrix&,imatrix&);
	friend bool	    Subset	      (imatrix&,imatrix&,double);
	friend imatrix	    Sup 	      (imatrix);
	friend imatrix	    Union	      (imatrix&,imatrix&);
	friend double	    Width	      (imatrix&);
	friend double	    Width	      (imatrix&,vector<int>&);
	friend double	    Width	      (imatrix&,imatrix&);
	friend imatrix          Zeros             (int,int);*/
};
#endif // __IMATRIX__
