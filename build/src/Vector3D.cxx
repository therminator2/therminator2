/********************************************************************************
 *                                                                              *
 *             THERMINATOR 2: THERMal heavy-IoN generATOR 2                     *
 *                                                                              *
 * Version:                                                                     *
 *      Release, 2.0.3, 1 February 2011                                         *
 *                                                                              *
 * Authors:                                                                     *
 *      Mikolaj Chojnacki   (Mikolaj.Chojnacki@ifj.edu.pl)                      *
 *      Adam Kisiel         (kisiel@if.pw.edu.pl)                               *
 *      Wojciech Broniowski (Wojciech.Broniowski@ifj.edu.pl)                    *
 *      Wojciech Florkowski (Wojciech.Florkowski@ifj.edu.pl)                    *
 *                                                                              *
 * Project homepage:                                                            *
 *      http://therminator2.ifj.edu.pl/                                         *
 *                                                                              *
 * For the detailed description of the program and further references           *
 * to the description of the model please refer to                              *
 * http://arxiv.org/abs/1102.0273                                               *
 *                                                                              *
 * This code can be freely used and redistributed. However if you decide to     *
 * make modifications to the code, please, inform the authors.                  *
 * Any publication of results obtained using this code must include the         *
 * reference to arXiv:1102.0273 and the published version of it, when           *
 * available.                                                                   *
 *                                                                              *
 ********************************************************************************/

#include "THGlobal.h"
#include "Vector3D.h"

// STATIC VARIABLES

double Vector3D::mXin;
double Vector3D::mXi0;
double Vector3D::mXip;
double Vector3D::mX0;
double Vector3D::mXd;

// OPERATORS

double& Vector3D::operator()(unsigned int i, unsigned int j, unsigned int k)
{
  return mVec[i][j][k];
}

// CLASS FUNCTIONS

Vector3D::Vector3D()
: mVecName(""), mVec(0),
mXmin(0.0), mXmax(1.0), mXpts(1),
mYmin(0.0), mYmax(1.0), mYpts(1),
mZmin(0.0), mZmax(1.0), mZpts(1)
{
  mVec		= new double** [1];
  mVec[0]	= new double*  [1];
  mVec[0][0]	= new double   [1];
  mVec[0][0][0]	= 0.0;
}

Vector3D::Vector3D(const char* aName,
		   double aXmin, double aXmax, int aXpts,
		   double aYmin, double aYmax, int aYpts,
		   double aZmin, double aZmax, int aZpts)
: mVecName(aName),
mXmin(aXmin), mXmax(aXmax), mXpts(aXpts),
mYmin(aYmin), mYmax(aYmax), mYpts(aYpts),
mZmin(aZmin), mZmax(aZmax), mZpts(aZpts)
{
  if(mXpts <= 0) mXpts = 1;
  if(mYpts <= 0) mYpts = 1;
  if(mZpts <= 0) mZpts = 1;
  
  mDi = (mXpts - 1) / (mXmax - mXmin);
  mDj = (mYpts - 1) / (mYmax - mYmin);
  mDk = (mZpts - 1) / (mZmax - mZmin);

  mVec = new double** [mXpts];
  for (int i=0; i<mXpts; i++) {
    mVec[i] = new double* [mYpts];
    for (int j=0; j<mYpts; j++) {
      mVec[i][j] = new double [mZpts];
      for (int k=0; k<mZpts; k++)
	mVec[i][j][k] = 0.0;
    }
  }
}

Vector3D::Vector3D(const Vector3D& aVector)
{
  mXmin = aVector.mXmin; mXmax = aVector.mXmax; mXpts = aVector.mXpts; mDi = aVector.mDi;
  mYmin = aVector.mYmin; mYmax = aVector.mYmax; mYpts = aVector.mYpts; mDj = aVector.mDj;
  mZmin = aVector.mZmin; mZmax = aVector.mZmax; mZpts = aVector.mZpts; mDk = aVector.mDk;

  mVecName = "";
  mVec = new double** [mXpts];
  for (int i=0; i<mXpts; i++) {
    mVec[i] = new double* [mYpts];
    for (int j=0; j<mYpts; j++) {
      mVec[i][j] = new double [mZpts];
      for (int k=0; k<mZpts; k++)
	mVec[i][j][k] = 0.0;
    }
  }
}

Vector3D::~Vector3D()
{
  if(mVec) {
    for (int i=0; i<mXpts; i++) {
      for (int j=0; j<mYpts; j++)
        delete[] mVec[i][j];
      delete[] mVec[i];
    }
    delete[] mVec;
  }
}

const char*	Vector3D::GetName() const	{ return mVecName.Data(); }
double		Vector3D::GetXMin() const	{ return mXmin; }
double		Vector3D::GetXMax() const	{ return mXmax; }
int		Vector3D::GetXPts() const	{ return mXpts; }
double		Vector3D::GetYMin() const	{ return mYmin; }
double		Vector3D::GetYMax() const	{ return mYmax; }
int		Vector3D::GetYPts() const	{ return mYpts; }
double		Vector3D::GetZMin() const	{ return mZmin; }
double		Vector3D::GetZMax() const	{ return mZmax; }
int		Vector3D::GetZPts() const	{ return mZpts; }

double Vector3D::Interpolate(double aX, double aY, double aZ)
{
  if(mZpts>1)
    return Interpolate3D(aX, aY, aZ);
  else if (mYpts>1)
    return Interpolate2D(aX, aY);
  else if (mXpts>1)
    return Interpolate1D(aX);
  else {
    return mVec[0][0][0];
  }
}

Vector3D* Vector3D::DerivativeX(const char* aName)
{
  int ti;
  Vector3D* tVec;

  mXd = (mXmax - mXmin) / (mXpts - 1);
  tVec = new Vector3D(*this);
  tVec->mVecName = aName;
  for (int i=0; i<mXpts; i++) {
    ti = InitDerivative(i, mXmin, mXmax, mXpts);
    for (int j=0; j<mYpts; j++)
      for (int k=0; k<mZpts; k++)
	(*tVec)(i,j,k) = Derivative(mVec[ti-1][j][k], mVec[ti][j][k], mVec[ti+1][j][k]);
  }
  return tVec;
}

Vector3D* Vector3D::DerivativeY(const char* aName)
{
  int tj;
  Vector3D* tVec;

  mXd = (mYmax - mYmin) / (mYpts - 1);
  tVec = new Vector3D(*this);
  tVec->mVecName = aName;
  for (int i=0; i<mXpts; i++)
    for (int j=0; j<mYpts; j++) {
      tj = InitDerivative(j, mYmin, mYmax, mYpts);
      for (int k=0; k<mZpts; k++)
	(*tVec)(i,j,k) = Derivative(mVec[i][tj-1][k], mVec[i][tj][k], mVec[i][tj+1][k]);
    }
  return tVec;
}

Vector3D* Vector3D::DerivativeZ(const char* aName)
{
  int tk;
  Vector3D* tVec;

  mXd = (mZmax - mZmin) / (mZpts - 1);
  tVec = new Vector3D(*this);
  tVec->mVecName = aName;
  for (int i=0; i<mXpts; i++)
    for (int j=0; j<mYpts; j++)
      for (int k=0; k<mZpts; k++) {
	tk = InitDerivative(k, mZmin, mZmax, mZpts);
	(*tVec)(i,j,k) = Derivative(mVec[i][j][tk-1], mVec[i][j][tk], mVec[i][j][tk+1]);
      }
  return tVec;
}

double Vector3D::Interpolate1D(double aX)
{
  int    i, j, k;
  double ti;

  ti = (aX - mXmin) * mDi;	i = (int) ti;	if(i+1 > mXpts-1) i--;	ti -= i;
  j  = 0;
  k  = 0;
  return
    mVec[i  ][j  ][k  ] * (1-ti) + mVec[i+1][j  ][k  ] * ti;
}

double Vector3D::Interpolate2D(double aX, double aY)
{
  int    i, j, k;
  double ti,tj;

  ti = (aX - mXmin) * mDi;	i = (int) ti;	if(i+1 > mXpts-1) i--;	ti -= i;
  tj = (aY - mYmin) * mDj;	j = (int) tj;	if(j+1 > mYpts-1) j--;	tj -= j;
  k  = 0;
  return
    (mVec[i  ][j  ][k  ] * (1-ti) + mVec[i+1][j  ][k  ] * ti) * (1-tj) +
    (mVec[i  ][j+1][k  ] * (1-ti) + mVec[i+1][j+1][k  ] * ti) *    tj;
}

double Vector3D::Interpolate3D(double aX, double aY, double aZ)
{
  int    i, j, k;
  double ti,tj,tk;

  ti = (aX - mXmin) * mDi;	i = (int) ti;	if(i+1 > mXpts-1) i--;	ti -= i;
  tj = (aY - mYmin) * mDj;	j = (int) tj;	if(j+1 > mYpts-1) j--;	tj -= j;
  tk = (aZ - mZmin) * mDk;	k = (int) tk;	if(k+1 > mZpts-1) k--;	tk -= k;
  return
    (
      (mVec[i  ][j  ][k  ] * (1-ti) + mVec[i+1][j  ][k  ] * ti) * (1-tj) +
      (mVec[i  ][j+1][k  ] * (1-ti) + mVec[i+1][j+1][k  ] * ti) *    tj
    ) * (1-tk) + (
      (mVec[i  ][j  ][k+1] * (1-ti) + mVec[i+1][j  ][k+1] * ti) * (1-tj) +
      (mVec[i  ][j+1][k+1] * (1-ti) + mVec[i+1][j+1][k+1] * ti) *    tj
    ) *    tk;
}

inline int Vector3D::InitDerivative(int aIdx, double aAMin, double aAMax, int aAPts)
{
  if(aIdx == 0) {
    mXin = aAMin;
    mXi0 = aAMin + mXd;
    mXip = aAMin + mXd * 2.0;
    mX0  = aAMin;
    return aIdx + 1;
  } else if (aIdx == aAPts - 1) {
    mXin = aAMax - mXd * 2.0;
    mXi0 = aAMax - mXd;
    mXip = aAMax;
    mX0  = aAMax;
    return aIdx - 1;
  } else {
    mXin = aAMin + mXd * (aIdx - 1);
    mXi0 = aAMin + mXd *  aIdx;
    mXip = aAMin + mXd * (aIdx + 1);
    mX0  = mXi0;
    return aIdx;
  }
}

inline double Vector3D::Derivative(double aFin, double aFi0, double aFip)
{
  return (
    mXin * (mXin - 2.0 * mX0) * (aFip - aFi0) +
    mXi0 * (mXi0 - 2.0 * mX0) * (aFin - aFip) +
    mXip * (mXip - 2.0 * mX0) * (aFi0 - aFin)
  ) / (2.0 * mXd * mXd * mXd);
}
