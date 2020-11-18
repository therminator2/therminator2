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

#ifndef _TH2_VECTOR_3D_H_
  #define _TH2_VECTOR_3D_H_

#include <TString.h>

class Vector3D
{
  public:
// OPERATORS
    double& operator()(unsigned int i, unsigned int j, unsigned int k);

// CLASS FUNCTIONS
    Vector3D();
    Vector3D(const char* aName,
	     double aXmin, double aXmax, int aXpts,
	     double aYmin, double aYmax, int aYpts,
	     double aZmin, double aZmax, int aZpts);
    Vector3D(const Vector3D& aVector);
    ~Vector3D();

    const char* GetName() const;
    double      GetXMin() const;
    double      GetXMax() const;
    int         GetXPts() const;
    double      GetYMin() const;
    double      GetYMax() const;
    int         GetYPts() const;
    double      GetZMin() const;
    double      GetZMax() const;
    int         GetZPts() const;

    double	Interpolate(double aX, double aY, double aZ);
    Vector3D*	DerivativeX(const char* aName);
    Vector3D*	DerivativeY(const char* aName);
    Vector3D*	DerivativeZ(const char* aName);

  private:
    double	Interpolate1D(double aX);
    double	Interpolate2D(double aX, double aY);
    double	Interpolate3D(double aX, double aY, double aZ);

    inline int	  InitDerivative(int aIdx, double aAMin, double aAMax, int aAPts);
    inline double Derivative(double aFin, double aFi, double aFip);
    
    TString	mVecName;
    double***	mVec;

    double	mXmin, mXmax;
    int		mXpts;
    double	mDi;

    double	mYmin, mYmax;
    int		mYpts;
    double	mDj;

    double	mZmin, mZmax;
    int		mZpts;
    double	mDk;
    
// used by InitDerivative() and Derivative()
    static double mXin, mXi0, mXip, mX0, mXd;
};

#endif

/*! @file Vector3D.h
 * @brief Definition of Vector3D class.
 */
/*! @class Vector3D
 * @brief Class representing a three dimensional vector field. 
 *
 * The table of numbers can be interpolated with a cubic method. Class also provides a simple calculation of the 
 * derivative along all three axes.
 *
 * @fn Vector3D::Vector3D()
 * @brief Default constructor.
 *
 * @fn Vector3D::Vector3D(const char* aName, double aXmin, double aXmax, int aXpts, double aYmin, double aYmax, int aYpts, double aZmin, double aZmax, int aZpts)
 * @brief Initiates a three dimensional vector field of a given name and limiting range and number of points for each axis.
 * @param [in] aName name of the vector field
 * @param [in] aXmin lower limit of X-axis
 * @param [in] aXmax higher limit of X-axis
 * @param [in] aXpts number of points in X direction
 * @param [in] aYmin lower limit of Y-axis
 * @param [in] aYmax higher limit of Y-axis
 * @param [in] aYpts number of points in Y direction
 * @param [in] aZmin lower limit of Z-axis
 * @param [in] aZmax higher limit of Z-axis
 * @param [in] aZpts number of points in Z direction
 *
 * @fn Vector3D::Vector3D(const Vector3D& aVector)
 * @brief Copying constructor.
 * @param [in] aVector referance to Vector3D
 *
 * @fn Vector3D::~Vector3D()
 * @brief Destructor.
 *
 * @fn double& Vector3D::operator()(unsigned int i, unsigned int j, unsigned int k)
 * @brief Operator (i,j,k) providing access to Vector3D array element
 * @param [in] i array index for the x-coordinate
 * @param [in] j array index for the y-coordinate
 * @param [in] k array index for the z-coordinate
 *
 * @fn const char* Vector3D::GetName() const
 * @brief Returns the vector field name
 *
 * @fn double Vector3D::GetXMin() const
 * @brief Returns the X-axis lower limit
 *
 * @fn double Vector3D::GetXMax() const
 * @brief Returns the X-axis higher limit
 *
 * @fn double Vector3D::GetXPts() const
 * @brief Returns the number of X-axis points
 *
 * @fn double Vector3D::GetYMin() const
 * @brief Returns the Y-axis lower limit
 *
 * @fn double Vector3D::GetYMax() const
 * @brief Returns the Y-axis higher limit
 *
 * @fn double Vector3D::GetYPts() const
 * @brief Returns the number of Y-axis points
 *
 * @fn double Vector3D::GetZMin() const
 * @brief Returns the Z-axis lower limit
 *
 * @fn double Vector3D::GetZMax() const
 * @brief Returns the Z-axis higher limit
 *
 * @fn double Vector3D::GetZPts() const
 * @brief Returns the number of Z-axis points
 * 
 * @fn double Vector3D::Interpolate(double aX, double aY, double aZ)
 * @brief Returns an interpolated value of the vector field.
 * @param [in] aX x-coordinate
 * @param [in] aY y-coordinate
 * @param [in] aZ z-coordinate
 *
 * @fn Vector3D* Vector3D::DerivativeX(const char* aName)
 * @brief Calculates the Derivative of the vector field along the X-axis
 * @param [in] aName name of the derived vector field
 * @retval Vector3D pointer to a created vector field
 *
 * @fn Vector3D* Vector3D::DerivativeY(const char* aName)
 * @brief Calculates the Derivative of the vector field along the Y-axis
 * @param [in] aName name of the derived vector field
 * @retval Vector3D pointer to a created vector field
 *
 * @fn Vector3D* Vector3D::DerivativeZ(const char* aName)
 * @brief Calculates the Derivative of the vector field along the Z-axis
 * @param [in] aName name of the derived vector field
 * @retval Vector3D pointer to a created vector field
 */
