// Simple interval library from Luc JAULIN, with minor modifications from Fabrice LE BARS and Jeremy NICOLA.

#include "imatrix.h"

using namespace TNT;
using namespace JAMA;

//----------------------------------------------------------------------
// Constructors/destructors
//----------------------------------------------------------------------
imatrix::imatrix()
{   
	Array2D <interval> A(1,1);
	data = A;
}
//----------------------------------------------------------------------
imatrix::imatrix(int i, int j)
{   
	Array2D <interval> A(j,i);
	data = A;
}
//----------------------------------------------------------------------
imatrix::imatrix(const rmatrix &X)
{      
	imatrix Res(X.dim1(),X.dim2());
	for (int i = 1; i <= X.dim1(); i++)
		for (int j = 1; j <= X.dim2(); j++)
			Res.SetVal(i,j,interval(X.GetVal(i,j)));
	data = Res.data;
}
//----------------------------------------------------------------------
imatrix::imatrix(const box &X)
{       
	imatrix Z(X.dim,1);
	for (int i = 1; i <= X.dim ; i++)
		Z.SetVal(i,1,X[i]);
	data = Z.data;
}
//----------------------------------------------------------------------
imatrix::imatrix(const imatrix &X)
{       
	if (&X == this) return;
	data = X.data;
}
//----------------------------------------------------------------------
imatrix::~imatrix () {}
//----------------------------------------------------------------------
// Operators
//----------------------------------------------------------------------
imatrix& imatrix::operator= (const imatrix& X)
{  	
	data = X.data;
	return *this;
}
//----------------------------------------------------------------------
imatrix operator+(const imatrix& X, const imatrix& Y)
{  
	imatrix Z(X);
	Z.data = X.data+Y.data;
	return Z;
}
//----------------------------------------------------------------------
imatrix operator-(const imatrix& X)
{    
	imatrix Z(X.dim1(),X.dim2());
	for (int i = 1; i <= X.dim1(); i++)
		for (int j = 1; j <= X.dim2(); j++)
			//Z.SetVal(i,j,-X.GetVal(i,j));
				Z.SetVal(i,j,-X(i,j));
	return Z;
}
//----------------------------------------------------------------------
imatrix operator-(const imatrix& X, const imatrix& Y)
{  
	imatrix Z(X);
	Z.data=X.data-Y.data;
	return Z;
}
//----------------------------------------------------------------------
imatrix operator*(const imatrix& X, const imatrix& Y)
{  
	imatrix Z(X.dim1(),Y.dim2());
	for (int i = 1; i <= X.dim1(); i++)
		for (int j = 1; j <= Y.dim2(); j++)
		{ 
			interval s = 0;
			for (int k = 1; k <= X.dim2(); k++)
				s = s+X(i,k)*Y(k,j);
			Z.SetVal(i,j,s);
			//Z(i,j)=s;
		}
		return Z;
}
//----------------------------------------------------------------------
box operator*(const imatrix& A, const box& x)
{ 
	box y(A.dim1());
	for (int i = 1; i <= A.dim1(); i++)
	{ 
		interval s = 0;
		for (int j = 1; j <= A.dim2(); j++)
			s = s+A(i,j)*x[j];
		//s = s+A(i,j)*x[j];
		y[i] = s;
	}
	return y;
}
//----------------------------------------------------------------------
imatrix operator*(const rmatrix& X, const imatrix& Y)
{  
	imatrix Z(X.dim1(),Y.dim2());
	for (int i = 1; i <= X.dim1(); i++)
		for (int j = 1; j <= Y.dim2(); j++)
		{ 
			interval s = 0;
			for (int k = 1; k <= X.dim2(); k++)
				s = s+X.GetVal(i,k)*Y.GetVal(k,j);
			Z.SetVal(i,j,s);
		}
		return Z;
}
//----------------------------------------------------------------------
imatrix operator*(const double a, const imatrix& X)
{  
	imatrix Z(X.dim1(),X.dim2());
	for (int i = 1; i <= X.dim1(); i++)
		for (int j = 1; j <= X.dim2(); j++)
			Z.SetVal(i,j,a*X(i,j));
	return Z;
}
//----------------------------------------------------------------------
/*interval& imatrix::operator() (int i, int j)  const
{ 
interval& I = data[j-1][i-1];
return I;
};*/
interval imatrix::operator() (int i, int j) const
{ 
	interval I = data[j-1][i-1];
	return I;
}
//----------------------------------------------------------------------
std::ostream& operator<<(std::ostream& os, const imatrix& X)
{  
	os << "matrix(" << X.dim1() << "," << X.dim2() << ")\n";
	for (int i=1; i<=X.dim1() ; i++)
	{
		for (int j = 1; j <= X.dim2(); j++)
			os << "  " << ": " << X.GetVal(i,j);
		os<<"\n";
	}
	return (os);
}
//----------------------------------------------------------------------
// Member functions
//----------------------------------------------------------------------
interval imatrix::GetVal(int i,int j) const
{ return data[j-1][i-1]; };
//----------------------------------------------------------------------
void imatrix::SetVal(int i,int j,interval X)
{ data[j-1][i-1] = X; };
//----------------------------------------------------------------------
int imatrix::dim1(void) const
{ return data.dim2(); };
//----------------------------------------------------------------------
int imatrix::dim2(void) const
{ return data.dim1(); };
//----------------------------------------------------------------------
// imatrix-valued functions
//----------------------------------------------------------------------
imatrix	iZeros(int n, int m)
{  
	rmatrix Z(n,m);
	for (int i = 1; i <= Z.dim1(); i++)
		for (int j = 1; j <= Z.dim2(); j++)
			Z.SetVal(i,j,0);
	return imatrix(Z);
}
//-----------------------------------------------------------------------
imatrix	iEye(int n)
{   
	rmatrix I(n,n);
	for (int i = 1; i <= I.dim1(); i++)
		for (int j = 1; j <= I.dim2(); j++)
			I.SetVal(i,j,0);
	for (int i = 1; i <= n; i++)
		I.SetVal(i,i,1);
	return imatrix(I);
}
//----------------------------------------------------------------------
imatrix Transpose(const imatrix& X)
{    
	imatrix Y(X.dim2(),X.dim1());
	for (int i = 1; i <= X.dim1(); i++)
		for (int j = 1; j <= X.dim2(); j++)
			Y.SetVal(j,i,X.GetVal(i,j));
	return Y;
}
//----------------------------------------------------------------------
imatrix RotationPhiThetaPsi(interval& phi, interval& theta, interval& psi)
{     
	interval Cphi,Sphi,Ctheta,Stheta,Cpsi,Spsi;
	interval a11,a12,a13,a21,a22,a23,a31,a32,a33;
	Cphi   = Cos(phi);     Sphi   = Sin(phi);
	Ctheta = Cos(theta);   Stheta = Sin(theta);
	Cpsi   = Cos(psi);     Spsi   = Sin(psi);
	a11=Ctheta*Cpsi;     a12=-Cphi*Spsi+Stheta*Cpsi*Sphi;  a13=Spsi*Sphi+Stheta*Cpsi*Cphi;
	a21=Ctheta*Spsi;     a22=Cpsi*Cphi+Stheta*Spsi*Sphi;   a23=-Cpsi*Sphi+Stheta*Cphi*Spsi;
	a31=-Stheta;         a32=Ctheta*Sphi;                  a33=Ctheta*Cphi;
	imatrix A(3,3);
	A.SetVal(1,1,a11); A.SetVal(1,2,a12); A.SetVal(1,3,a13);
	A.SetVal(2,1,a21); A.SetVal(2,2,a22); A.SetVal(2,3,a23);
	A.SetVal(3,1,a31); A.SetVal(3,2,a32); A.SetVal(3,3,a33);
	return(A);
}
//----------------------------------------------------------------------
// Other functions
//----------------------------------------------------------------------
rmatrix Center(const imatrix& X)
{  
	rmatrix Z(X.dim1(),X.dim2());
	for (int i = 1; i <= X.dim1(); i++)
		for (int j = 1; j <= X.dim2(); j++)
			Z.SetVal(i,j,Center(X.GetVal(i,j)));
	return Z;
}
//-----------------------------------------------------------------------
box Row(const imatrix& B, int i)
{    
	box Z(B.dim2());
	for (int j = 1; j <= B.dim2(); j++)
		//Z[j] = B.GetVal(i,j);
			Z[j] = B(i,j);
	return Z;
}
//----------------------------------------------------------------------
box Column(const imatrix& B, int j)
{   
	box Z(B.dim1());
	for (int i = 1; i <= B.dim1(); i++)
		//Z[i]=B.GetVal(i,j);
			Z[i] = B(i,j);
	return Z;
}
//----------------------------------------------------------------------
box ToBox(const imatrix& B)
{    
	box Z(B.dim1());
	for (int i = 1; i <= B.dim1(); i++)
		//Z[i] = B.GetVal(i,1);
			Z[i] = B(i,1);
	return Z;
}
//----------------------------------------------------------------------
// Contractors
//----------------------------------------------------------------------
void Cmul(imatrix& C, imatrix& A, imatrix& B)
{   
	for (int i = 1; i <= C.dim1(); i++)
		for (int j = 1; j <= C.dim2(); j++)
		{ 
			box a = Row(A,i);
			box b = Column(B,j);
			interval c = C.GetVal(i,j);
			Cscal(c, a, b);
			C.SetVal(i,j,c);
			for (int k = 1; k <= A.dim2(); k++)
			{
				A.SetVal(i,k,a[k]);
				B.SetVal(k,j,b[k]);
			}
		}
}
//----------------------------------------------------------------------
void Cmul(box& c, imatrix& A, box& b)
{   
	for (int i = 1; i <= c.dim; i++)
	{ 
		box a = Row(A,i);
		Cscal(c[i],a,b);
		for (int k = 1; k <= b.dim; k++)
			A.SetVal(i,k,a[k]);
	}
}
//----------------------------------------------------------------------
void Crot(imatrix& R)
{   
	imatrix Rt = Transpose(R);
	imatrix I = iEye(R.dim1());
	Cmul(I,R,Rt);
}
//----------------------------------------------------------------------
void Cantisym(imatrix& A)
{   
	for (int i = 1; i <= A.dim1(); i++)
	{ 
		A.SetVal(i,i,interval(-0,0));
		for (int j = i+1; j <= A.dim1(); j++)
		{  
			A.SetVal(j,i,Inter(-A(i,j),A(j,i)));
			A.SetVal(i,j,Inter(A(i,j),-A(j,i)));
		}
	}  
}
//----------------------------------------------------------------------
