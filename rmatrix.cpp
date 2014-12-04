// Simple interval library from Luc JAULIN, with minor modifications from Fabrice LE BARS and Jeremy NICOLA.

#include "rmatrix.h"

using namespace TNT;
using namespace JAMA;

//----------------------------------------------------------------------
// Constructors/destructors
//----------------------------------------------------------------------
rmatrix::rmatrix()
{   
	Array2D <double> A(1,1);
	data = A;
}
//----------------------------------------------------------------------
rmatrix::rmatrix(int i, int j)
{   
	Array2D <double> A(j,i);
	data = A;
}
//----------------------------------------------------------------------
rmatrix::rmatrix (const box &X)
{       
	rmatrix Z(X.dim,1);
	for (int i = 1; i <= X.dim; i++)
		Z.SetVal(i, 1, Center(X[i]));
	data = Z.data;
}
//----------------------------------------------------------------------
rmatrix::rmatrix(const rmatrix &X)
{      
	if (&X == this) return;
	data = X.data;
}
//----------------------------------------------------------------------
rmatrix::~rmatrix() {}
//----------------------------------------------------------------------
// Operators
//----------------------------------------------------------------------
rmatrix& rmatrix::operator=(const rmatrix& X)
{  	
	data = X.data;
	return *this;
}
//----------------------------------------------------------------------
rmatrix operator+(const rmatrix& X, const rmatrix& Y)
{  
	rmatrix Z(X);
	Z.data = X.data+Y.data;
	return Z;
}
//----------------------------------------------------------------------
rmatrix operator-(const rmatrix& X)
{    
	rmatrix Z(X.dim1(),X.dim2());
	for (int i = 1; i <= X.dim1(); i++)
		for (int j = 1; j <= X.dim2(); j++)
			Z.SetVal(i,j,-X.GetVal(i,j));
	return Z;
}
//----------------------------------------------------------------------
rmatrix operator-(const rmatrix& X, const rmatrix& Y)
{  
	rmatrix Z(X);
	Z.data = X.data-Y.data;
	return Z;
}
//----------------------------------------------------------------------
rmatrix operator*(const rmatrix& X, const rmatrix& Y)
{  
	rmatrix Z(X.dim1(),Y.dim2());
	for (int i = 1; i <= X.dim1(); i++)
		for (int j = 1; j <= Y.dim2(); j++)
		{ 
			double s = 0;
			for (int k = 1; k <= X.dim2(); k++)
				s = s+X.GetVal(i,k)*Y.GetVal(k,j);
			Z.SetVal(i,j,s);
		}
		return Z;
}
//----------------------------------------------------------------------
rmatrix operator*(const double a, const rmatrix& X)
{  
	rmatrix Z(X.dim1(),X.dim2());
	for (int i = 1; i <= X.dim1(); i++)
		for (int j = 1; j <= X.dim2(); j++)
			Z.SetVal(i,j,a*X.GetVal(i,j));
	return Z;
}
//----------------------------------------------------------------------
/*
rmatrix operator*(const rmatrix& X, const double& a)
{ return (a*X); }
//----------------------------------------------------------------------
std::ostream& operator<<(std::ostream& os, const rmatrix& X)
{   
cout << "rmatrix :"<<"\t dim="<<Size(X)<<"\n";
if (Emptyrmatrix(X)) os<<"Emptyrmatrix";
for (int i=1; i<=Size(X); i++)
os << "  " <<i<< ": "<<X[i]<< "\n";
return (os);
}*/
//----------------------------------------------------------------------
// Member functions
//----------------------------------------------------------------------
double rmatrix::GetVal(int i,int j) const
{ return data[j-1][i-1]; };
//----------------------------------------------------------------------
void rmatrix::SetVal(int i,int j,double X)
{ data[j-1][i-1]=X; };
//----------------------------------------------------------------------
int rmatrix::dim1(void) const
{ return data.dim2(); };
//----------------------------------------------------------------------
int rmatrix::dim2(void) const
{ return data.dim1(); };
//----------------------------------------------------------------------
// rmatrix-valued functions
//----------------------------------------------------------------------
rmatrix	Zeros(int n, int m)
{   
	rmatrix Z(n,m);
	for (int i=1; i<=Z.dim1() ; i++)
		for (int j=1; j<=Z.dim2() ; j++)
			Z.SetVal(i,j,0);
	return Z;
}
//----------------------------------------------------------------------
rmatrix	Eye(int n)
{   
	rmatrix Z=Zeros(n,n);
	for (int i=1; i<=n ; i++)
		Z.SetVal(i,i,1);
	return Z;
}
//----------------------------------------------------------------------
rmatrix	Inv(rmatrix& A)
{ 
	rmatrix Z;
	LU <double> A1(A.data);
	rmatrix Id=Eye(A.dim1());
	Z.data = A1.solve(Id.data);
	return Z;
}
//----------------------------------------------------------------------
rmatrix RotationPhiThetaPsi(double phi,double theta,double psi)
{      
	/* r: angle de roulis (roll) en radians */
	/* p: angle de tangage (pitch) en radians */
	/* h: angle de cap (heading) en radians */
	/* Les rotations sont appliquée dans l'ordre suivant: */
	/*      1: rotation autour de l'axe de roulis */
	/*      2: rotation autour de l'axe de tangage */
	/*      3: rotation autour de l'axe de lacet */
	/*      {x,y,z} = Ch Cp Cr {x,y,z} */
	rmatrix A(3,3);
	double cphi,sphi,ctheta,stheta,cpsi,spsi;
	cphi   = cos(phi);        sphi   = sin(phi);
	ctheta = cos(theta);      stheta = sin(theta);
	cpsi   = cos(psi);        spsi   = sin(psi);
	A.SetVal(1,1,ctheta*cpsi); A.SetVal(1,2,-cphi*spsi+stheta*cpsi*sphi); A.SetVal (1,3,spsi*sphi+stheta*cpsi*cphi);
	A.SetVal(2,1,ctheta*spsi); A.SetVal(2,2,cpsi*cphi+stheta*spsi*sphi);  A.SetVal (2,3,-cpsi*sphi+stheta*cphi*spsi);
	A.SetVal (3,1,-stheta);    A.SetVal(3,2,ctheta*sphi);            A.SetVal (3,3,ctheta*cphi);
	return(A);
}
//----------------------------------------------------------------------
// Other functions
//----------------------------------------------------------------------
box ToBox(const rmatrix& B)
{    
	box Z(B.dim1());
	for (int i=1; i<=B.dim1() ; i++)
		Z[i]=B.GetVal(i,1);
	return Z;
}
double Det(rmatrix& A)
{ 
	LU <double> A1(A.data);
	return A1.det();
}
//----------------------------------------------------------------------
//double determinant_de_A =  A1.det();
//Array1D < double > b(2);
//b[0] = 1;
//b[1] = 1;

//Array1D < double > x(2);
//x = A1.solve(b);

/*
void Update(rmatrix& X)
{ for (int i=1; i<=Size(X); i++)
{ if (X[i].isEmpty)
{ for (int j=1; j<=Size(X); j++)
X[j]=double(); return; }
}
}
//----------------------------------------------------------------------
bool Emptyrmatrix(const rmatrix& X)
{  if (X.dim==0) return true;
for (int i=1; i<=X.dim; i++)
{ if (X[i].isEmpty)  return true;}
return false;
}
//----------------------------------------------------------------------
int Size(const rmatrix& X) {return (X.dim);}
//----------------------------------------------------------------------
rmatrix Inf(rmatrix X)
{   rmatrix Ans(Size(X));
for (int k=1; k<=Size(Ans); k++)  Ans[k]=X[k].inf;
return Ans;
}
//----------------------------------------------------------------------
rmatrix Sup(rmatrix X)
{   rmatrix Ans(Size(X));
for (int k=1; k<=Size(Ans); k++) Ans[k]=X[k].sup;
return Ans;
}
//----------------------------------------------------------------------
rmatrix Concat(const rmatrix& X, const rmatrix& Y)
//  Produit Cartesien ou Concaténation de deuX pavés X et y :
//      Ans=[X,Y]     =>     Ans=Concat(X,Y);
{  double dim=X.dim+Y.dim; rmatrix Ans(dim);
if ((!Emptyrmatrix(X))&&(!Emptyrmatrix(Y)))
{ for (int i=1; i<=dim; i++)
{ if (i<=Size(X)) Ans[i]=X[i]; else Ans[i]=Y[i-X.dim]; }
}
return Ans;
}
//----------------------------------------------------------------------
rmatrix Proj(const rmatrix& X, int i, int j)
// Projection du pavé X dans un espace de dimension dim=(j-i)+1;
//   X=[[X1],[X2],..,[Xi],..,[Xj],..[Xn]]
//   =>  Proj(X,i,j)=[[Xi],..,[Xj]] et Proj(X,i,i)=[Xi]
{  int dim=abs(j-i)+1; rmatrix Ans(dim);
if (!Emptyrmatrix(X))
{ int lb=min(i,j);
for (int k=1; k<=dim; k++) Ans[k]=X[k+lb-1]; }
return Ans;
}
//----------------------------------------------------------------------
int AxePrincipal(rmatrix& X)
{ int kmax=1;
double widthmax=Width(X[kmax]);
for (int k=2; k<=Size(X); k++)
{ if (Width(X[k])>widthmax)
{kmax=k; widthmax=Width(X[k]);}
}
return kmax;
}
//----------------------------------------------------------------------
int AxePrincipal(rmatrix& X, rmatrix& Y)
{ int kmax=1;
double widthmax=Width(X[kmax])*fabs(Y[kmax].inf);
for (int k=2; k<=Size(X); k++)
{	if (Width(X[k])*fabs(Y[k].inf)>widthmax)
{ kmax=k; widthmax=Width(X[k])*fabs(Y[k].inf);}
}
return kmax;
}
//----------------------------------------------------------------------
int AxePrincipal(rmatrix& X, vector<int>& v)
{ int kmax=v[1], SizeV=v.size()-1;
double widthmax=Width(X[kmax]);
for (int k=2; k<=SizeV; k++)
{ if (Width(X[v[k]])>widthmax) {kmax=v[k]; widthmax=Width(X[v[k]]);}}
return kmax;
}
//----------------------------------------------------------------------
void Bisect(rmatrix& X, rmatrix& X1, rmatrix& X2)
{ BisectAlong(X,X1,X2,AxePrincipal(X));	}
//----------------------------------------------------------------------
void Bisect(rmatrix& X, rmatrix& X1, rmatrix& X2, vector<int>& v)
{ BisectAlong(X,X1,X2,AxePrincipal(X,v)); }
//----------------------------------------------------------------------
void Bisect(rmatrix& X, rmatrix& X1, rmatrix& X2, rmatrix& V)
{   BisectAlong(X,X1,X2,AxePrincipal(X,V));  }
//----------------------------------------------------------------------
void BisectAlong(rmatrix& X, rmatrix& X1, rmatrix& X2, int i)
{  X1=X2=X; X1[i].sup=X2[i].inf=Center(X[i]); }
//----------------------------------------------------------------------
void BisectHere(rmatrix& X, rmatrix& X1, rmatrix& X2, int i, double here) //Used by the 3-B consistency
{	X1=X2=X;
X1[i].sup=here;
X2[i].inf=here;
}
//----------------------------------------------------------------------
void TrisectAlong(rmatrix& X, rmatrix& X1, rmatrix& X2, rmatrix& X3, int i)
{  X1=X2=X3=X;
X1[i].sup=X[i].inf+(Width(X[i])/3);
X2[i].inf=X1[i].sup;
X2[i].sup=X2[i].inf+(Width(X[i])/3);
X3[i].inf=X2[i].sup;
}
//-----------------------------------------------------------------------
void Trisect(rmatrix& X, rmatrix& X1, rmatrix& X2, rmatrix& X3)
{ TrisectAlong(X,X1,X2,X3,AxePrincipal(X)); }
//-----------------------------------------------------------------------
void DecoupAlong(rmatrix& X, rmatrix& X1, rmatrix& X2, int i)
{  X1=X2=X;
X1[i].sup=X[i].inf+(Width(X[i])/3);
X2[i].inf=X1[i].sup;
}
//------------------------------------------------------------------------
void Decoup(rmatrix& X, rmatrix& X1, rmatrix& X2)
{ DecoupAlong(X,X1,X2,AxePrincipal(X));}
//-----------------------------------------------------------------------
rmatrix Rand(const rmatrix& X)
{ int sizeX=Size(X); rmatrix Ans(sizeX);
if (Emptyrmatrix(X)) Ans=Empty(sizeX);
else { for (int k=1; k<=sizeX; k++) Ans[k]=Rand(X[k]); }
return Ans;
}
//-----------------------------------------------------------------------
rmatrix Center(const rmatrix& X)
{ int sizeX=Size(X); rmatrix Ans(sizeX);
if (Emptyrmatrix(X)) Ans=Empty(sizeX);
else { for (int k=1; k<=sizeX; k++) Ans[k]=Center(X[k]); }
return Ans;
}
//-----------------------------------------------------------------------
rmatrix Center(const rmatrix& X, vector<int>& v)
{ int sizev=v.size()-1; rmatrix Ans=X;
for (int k=1; k<=sizev; k++) Ans[v[k]]=Center(X[v[k]]);
return Ans;
}
//-----------------------------------------------------------------------
bool Disjoint (const rmatrix& X, const rmatrix& Y)
{   if (Emptyrmatrix(X)||Emptyrmatrix(Y)) return true;
for (int i=1; i<=Size(X); i++)
if (Disjoint(X[i],Y[i])) return true;
return false;
}
//----------------------------------------------------------------------
iboolean In(rmatrix X, rmatrix Y)
{  if (Emptyrmatrix(X)||Emptyrmatrix(Y)) return ifalse;
iboolean r=itrue;
for (int k=1; k<=X.dim; k++)
{ iboolean in1=In(X[k],Y[k]);
if (in1==false) return ifalse;
else if (in1.value==iperhaps) r=iperhaps; }
return r;
}
//------------------------------------------------------------------------------
double decrease(const rmatrix& X, const rmatrix& Y)
{    double e=0;
for (int k=1; k<=X.dim; k++)
{  if ((X[k].isEmpty)||(Y[k].isEmpty)) return (-1);
double e1=0;
double Xinf=X[k].inf, Xsup=X[k].sup;
double Yinf=Y[k].inf, Ysup=Y[k].sup;
if (Xsup>=Ysup) e1=max(e1,Xsup-Ysup);
if (Xinf<=Yinf) e1=max(e1,Yinf-Xinf);
if (Xsup<Ysup) e1=max(e1,Ysup-Xsup);
if (Xinf>Yinf) e1=max(e1,Xinf-Yinf);
e=max(e,e1);
}
return e;
}
//----------------------------------------------------------------------
double decrease(const rmatrix& X, const rmatrix& Y, vector<int> Tab)
{    double e=0; int size=Tab.size()-1;
for (int k=1; k<=size; k++)
{  if ((X[Tab[k]].isEmpty)||(Y[Tab[k]].isEmpty)) return (-1);
double e1=0;
double Xinf=X[Tab[k]].inf, Xsup=X[Tab[k]].sup;
double Yinf=Y[Tab[k]].inf, Ysup=Y[Tab[k]].sup;
if (Xsup>=Ysup) e1=max(e1,Xsup-Ysup);
if (Xinf<=Yinf) e1=max(e1,Yinf-Xinf);
if (Xsup<Ysup) e1=max(e1,Ysup-Xsup);
if (Xinf>Yinf) e1=max(e1,Xinf-Yinf);
e=max(e,e1);
}
return e;
}
//----------------------------------------------------------------------
double Eloignement(rmatrix& X, rmatrix& Y)
{  if ((Emptyrmatrix(X))||(Emptyrmatrix(Y))) return oo;
double e=0;
for (int k=1; k<=Size(X); k++) {e=max(e,Eloignement(X[k],Y[k]));}
return e;
}
//------------------------------------------------------------------------------
double Eloignement2(rmatrix& X, rmatrix& Y)
{// prend le point X1 de X qui est le plus eloigne de [Y] et renvoie la
// distance de X1 avec [Y]
Update(X); Update(Y);
if (Subset(X,Y)) return 0;
double e=0;
for (int k=1; k<=X.dim; k++)
{double e1=0;
double Xinf=X[k].inf;
double Xsup=X[k].sup;
double Yinf=Y[k].inf;
double Ysup=Y[k].sup;
if (Xsup>Ysup) e1=max(e1,Xsup-Ysup);
if (Xinf<Yinf) e1=max(e1,Yinf-Xinf);
e=max(e,e1);
}
return e;
}
//------------------------------------------------------------------------------
double EloignementRelatif2(rmatrix& X,rmatrix& Y)
{// prend le point X1 de X qui est le plus eloigne de [Y] et renvoie la
// distance de X1 avec [Y]
Update(X); Update(Y);
if (Subset(X,Y)) return 0;
if (Emptyrmatrix(Y)) return oo;
double e=0;
for (int k=1; k<=X.dim; k++)
{double e1=0;
double Xinf=X[k].inf;
double Xsup=X[k].sup;
double Yinf=Y[k].inf;
double Ysup=Y[k].sup;
if (Xsup>Ysup) e1=max(e1,fabs((Xsup-Ysup)/Ysup));
if (Xinf<Yinf) e1=max(e1,fabs((Yinf-Xinf)/Yinf));
e=max(e,e1);
}
return e;
}
//----------------------------------------------------------------------
rmatrix Inter(const rmatrix& X, const rmatrix& Y)
{   rmatrix Ans(Size(X));
if ((Emptyrmatrix(X))||(Emptyrmatrix(Y))) {return Ans;}
for (int k=1; k<=Size(Ans); k++)
{ Ans[k]=Inter(X[k],Y[k]);
if (Ans[k].isEmpty) {Update(Ans); return Ans;}
}
return Ans;
}
//----------------------------------------------------------------------
void Phi0(rmatrix& B0,rmatrix* pP) //0-intersection (voir qminimax)
{   if (pP==NULL) {return;}
B0=Inter(B0,*pP);
Phi0(B0,pP->next);
}
//----------------------------------------------------------------------
void Phi1(rmatrix& B0, rmatrix& B1, rmatrix* pP) //1-intersection (voir qminimax)
{// Hyp: B1 contient B0
if (pP==NULL) {return;}
rmatrix B11=Inter(*pP,B1);
B1=Union(B0,B11);
B0=Inter(B0,*pP);
Phi1(B0,B1,pP->next);
}
//----------------------------------------------------------------------
void Phi2(rmatrix& B0, rmatrix& B1,rmatrix& B2,rmatrix* pP) //2-intersection
{// Hyp: B2 contient B1 qui contient B0
if (pP==NULL) {return;}
rmatrix B21=Inter(*pP,B2);   B2=Union(B1,B21);
rmatrix B11=Inter(*pP,B1);   B1=Union(B0,B11);
B0=Inter(B0,*pP);
Phi2(B0,B1,B2,pP->next);
}
//----------------------------------------------------------------------
rmatrix Inflate(rmatrix& X, double eps)
{   Update(X);
rmatrix Ans(Size(X));
for (int k=1; k<=Size(Ans); k++)  {Ans[k]=Inflate(X[k],eps);}
return Ans;
}
//----------------------------------------------------------------------
bool Prop(rmatrix& X, rmatrix& Y)  // Normalement X is a subset of y (used in SIVEX)
{ if (Emptyrmatrix(X)) return false;
if (Emptyrmatrix(Y)) return false;
for (int k=1; k<=Size(X); k++)
if ((X[k].inf==Y[k].inf)||(X[k].sup==Y[k].sup)) return (true);
return false;
}
//----------------------------------------------------------------------
rmatrix Union(rmatrix& X, rmatrix& Y)
{   rmatrix Ans(max(Size(X),Size(Y)));
if (Emptyrmatrix(X)) return (Y);
if (Emptyrmatrix(Y)) return (X);
for (int k=1; k<=Size(Ans); k++)  Ans[k]=Union(X[k],Y[k]);
return Ans;
}
//----------------------------------------------------------------------
bool Subset(rmatrix& X, rmatrix& Y)
{   if (Emptyrmatrix(X)) return true;
if (Emptyrmatrix(Y)) return false;
bool b=true;
for (int k=1; k<=Size(X); k++) b = b && Subset(X[k],Y[k]);
return (b);
}
//----------------------------------------------------------------------
bool Subset(rmatrix& X, rmatrix& Y,double epsilon)
{   if (Subset(X,Y)==false) return false;
for (int k=1; k<=Size(X); k++)
if (Subset(X[k],Y[k],epsilon)) return true;
return (false);
}
//----------------------------------------------------------------------
void Sucre(rmatrix& P, rmatrix& S)    //Sucre le bout de P qui est dans S
{   int sum=0;
int j=1;
for (int i=1; i<=P.dim; i++)
{if (Subset(P[i],S[i])) sum++; else j=i;}
if (sum==(P.dim)-1)
if (!Subset(S[j],P[j])) Sucre(P[j],S[j]);
};
//----------------------------------------------------------------------
double Volume(rmatrix& X)
{   double vol=1;
for (int i=1;i<=Size(X);i++) vol=vol*Width(X[i]);
return vol;
};
//----------------------------------------------------------------------
double Width(rmatrix& X) {
return Width(X[AxePrincipal(X)]);}
//----------------------------------------------------------------------
double Width(rmatrix& X, vector<int>& v)  {
return Width(X[AxePrincipal(X,v)]);}
//----------------------------------------------------------------------
double Width(rmatrix& X, rmatrix& Y)
{ int i=AxePrincipal(X,Y);
return (Width(X[i])*fabs(Y[i].inf));	}
//----------------------------------------------------------------------
double Marge(rmatrix X, rmatrix Y)
{  if ((Emptyrmatrix(X))||(Emptyrmatrix(Y))) return -oo;
double ans=Marge(X[1],Y[1]);
for (int i=2; i<=Size(X); i++)
ans=min(ans,Marge(X[i],Y[i]));
return ans;
}
//-----------------------------------------------------------------------
double ProduitScalaire(rmatrix& U, rmatrix& V)
{  double sum=0;
for (int i=1; i<=Size(U); i++)  sum=sum+U[i]*V[i];
return (sum);
}
//----------------------------------------------------------------------
double Determinant(rmatrix& U, rmatrix& V)
{   double u1=U[1];
double v2=V[2];
double v1=V[1];
double u2=U[2];
double r=u1*v2-v1*u2;
return u1*v2-v1*u2; }
// Calcule l'angle entre 2 vecteur de dimension 2
double Angle (rmatrix& V, rmatrix& W)
// Attention, il faut des vecteurs et non des rmatrix et les vecteurs doivent etre
// de dimension 2
{   if ((Size(V)!=2)||(Size(W)!=2))
{
cout<<"error";
}
double n2,costeta,sinteta;
double nv=Norm(V);
double nw=Norm(W);
double nvw=nv*nw;
n2=ToReel(nvw);
costeta=ToReel(ProduitScalaire(V,W)/n2);
sinteta=ToReel(Determinant(V,W))/n2;
return(Arccossin(costeta,sinteta));
}
//----------------------------------------------------------------------
double Norm(rmatrix X)
{ if (Emptyrmatrix(X)) return double();
double r=0;
for (int i=1; i<=Size(X); i++) r=r+Sqr(X[i]);
return (Sqrt(r));
}
//----------------------------------------------------------------------
double NormEuclid(rmatrix X, rmatrix Y)
{ if (Size(X)!=Size(Y)) return double();
if (Emptyrmatrix(X)||Emptyrmatrix(Y)) return double();
double r=0;
for (int i=1; i<=Size(X); i++) r=r+Sqr(Y[i]-X[i]);
return (Sqrt(r));
}
//----------------------------------------------------------------------
double NormInf(rmatrix X, rmatrix Y)
{ if (Size(X)!=Size(Y)) return double();
if (Emptyrmatrix(X)||Emptyrmatrix(Y)) return double();
double ans=Abs(Y[1]-X[1]);
for (int i=1; i<=Size(X); i++) ans=Max(ans,Abs(Y[i]-X[i]));
return ans;
}
//----------------------------------------------------------------------
rmatrix Zeros(int d)
{   rmatrix Ans(d);
for (int k=1; k<=d; k++) Ans[k]=0;
return Ans;
}
//----------------------------------------------------------------------
bool Isrmatrix(rmatrix X)
{ if (Emptyrmatrix(X)) return false;
for (int i=1; i<=Size(X); i++)
{ if (Width(X[i])==0) return false; }
return true;
}
//---------------------------------------------------------------------*/

