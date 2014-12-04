// Simple interval library from Luc JAULIN, with minor modifications from Fabrice LE BARS and Jeremy NICOLA.

#ifndef __RMATRIX__
#define __RMATRIX__

#include "box.h"
#include "matrix_lib/tnt.h"
#include "matrix_lib/jama_lu.h"

class rmatrix
{
public:
	TNT::Array2D <double> data;

public:
	//----------------------------------------------------------------------
	// Constructors/destructors
	//----------------------------------------------------------------------
	rmatrix();
	rmatrix(int, int);
	rmatrix(const box &);
	rmatrix(const rmatrix&);
	~rmatrix();
	//----------------------------------------------------------------------
	// Operators
	//----------------------------------------------------------------------
	rmatrix& operator=(const rmatrix&);
	friend rmatrix operator+(const rmatrix&, const rmatrix&);
	friend rmatrix operator-(const rmatrix&);
	friend rmatrix operator-(const rmatrix&, const rmatrix&);
	friend rmatrix operator*(const rmatrix&, const rmatrix&);
	friend rmatrix operator*(const double, const rmatrix&);
	//friend rmatrix operator*(const rmatrix&, const interval&);
	//friend std::ostream& operator<<(std::ostream&, const rmatrix&);
	//----------------------------------------------------------------------
	// Member functions
	//----------------------------------------------------------------------
	double GetVal(int, int) const;
	void SetVal(int, int, double);
	int dim1(void) const;
	int dim2(void) const;
};

//----------------------------------------------------------------------
// rmatrix-valued functions
//----------------------------------------------------------------------
rmatrix Zeros(int, int);
rmatrix Eye(int);
rmatrix Inv(rmatrix&);
rmatrix RotationPhiThetaPsi(double phi, double theta, double psi);
//----------------------------------------------------------------------
// Other functions
//----------------------------------------------------------------------
box ToBox(const rmatrix& B);
double Det(rmatrix&);

/*double       Angle             (rmatrix& ,rmatrix&); // Il faut des vecteurs de dim 2
int          Size              (const rmatrix&);
int	    AxePrincipal      (rmatrix&);
int	    AxePrincipal      (rmatrix&, vector<int>&);
int	    AxePrincipal      (rmatrix&, rmatrix&);
void         Update            (rmatrix&);
rmatrix          Rand              (const rmatrix& X);
rmatrix	    Center	      (const rmatrix&);
rmatrix          Center            (const rmatrix&, vector<int>&);
//void	    CheckRange	      (rmatrix&,rmatrix&);
interval     Determinant       (rmatrix&, rmatrix&);
bool         Emptyrmatrix          (const rmatrix&);
bool	    Disjoint	      (const rmatrix&,const rmatrix&);
double         decrease          (const rmatrix&, const rmatrix&);
double         decrease          (const rmatrix&, const rmatrix&, vector<int>);
double         Eloignement       (rmatrix&,rmatrix&);
double         Eloignement2      (rmatrix&,rmatrix&);
double         EloignementRelatif2(rmatrix&,rmatrix&);
double         Marge             (rmatrix,rmatrix);
iboolean	    In		      (rmatrix,rmatrix);
rmatrix	    Inf 	      (rmatrix);
rmatrix	    Inflate 	      (rmatrix&,double);
rmatrix	    Inter 	      (const rmatrix&,const rmatrix&);
rmatrix          Concat            (const rmatrix&, const rmatrix&);
rmatrix          Proj              (const rmatrix&, int, int);
//void       Inter1            (rmatrix&,rmatrix&,const rmatrix&,const rmatrix&,const rmatrix&);
interval     Norm              (rmatrix);
interval     NormEuclid        (rmatrix, rmatrix);
interval     NormInf           (rmatrix, rmatrix);
void         Phi0              (rmatrix&,rmatrix*);             //0-intersection (voir qminimax)
void         Phi1              (rmatrix&,rmatrix&,rmatrix*);       //1-intersection (voir qminimax)
void         Phi2              (rmatrix&,rmatrix&,rmatrix&,rmatrix*);  //2-intersection
interval     ProduitScalaire   (rmatrix&,rmatrix&);
bool         Prop              (rmatrix&,rmatrix&);
bool	    Subset	      (rmatrix&,rmatrix&);
bool	    Subset	      (rmatrix&,rmatrix&,double);
bool         Isrmatrix             (rmatrix);
void         Sucre             (rmatrix&,rmatrix&);
rmatrix	    Sup 	      (rmatrix);
rmatrix	    Union	      (rmatrix&,rmatrix&);
double         Volume            (rmatrix&);
double	    Width	      (rmatrix&);
double	    Width	      (rmatrix&,vector<int>&);
double	    Width	      (rmatrix&,rmatrix&);
//rmatrix          Empty             (int);
//rmatrix          Infinity          (int);*/

#endif // __RMATRIX__
