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

#ifndef _TH2_HYPERSURFACE_H_
  #define _TH2_HYPERSURFACE_H_

#include <TRandom2.h>
#include "Thermodynamics.h"
#include "Vector3D.h"

class Hypersurface_Lhyquid3D
{
  public:
    Hypersurface_Lhyquid3D();
    ~Hypersurface_Lhyquid3D();

    double GetHyperCubeSpatialVolume();
    double GetDSigmaP(double aMt, double aPt, double aPhiP, double aRapP);
    double GetPdotU(  double aMt, double aPt, double aPhiP, double aRapP);  
    void   SetPositionOnHypersurface(TRandom2* aRandom, double* aXt, double* aXx, double* aXy, double* aXz);
    
  protected:
    double mRapPRange;
    double mLambda;
    double mTauI;
    double mCollidingEnergy;
    double mCentralityMin;
    double mCentralityMax;
    double mImpactParameter;
    double mTempI;
    char   mDeviceName[30];
    char   mCollidingSystem[30];
    Thermodynamics* mThermo;

  private:
    void ReadParameters();

// hypersurface coordinates
    double Zeta, PhiS, Theta;
    double Tau,  Rho,  RapS;
// value of the distance to the hypersurface element and the derivatives
    double Dhs, dDdZeta, dDdPhi, dDdTheta;
// value of the fluid 4-velocity u_x and u_y and rapidity on the hypersurface
    double Ux, Uy, RapF;    
// hypersurface distance and derivatives
    Vector3D* mDistance;
    Vector3D* mDistanceDPhi;
    Vector3D* mDistanceDZeta;
    Vector3D* mDistanceDTheta;
// four-velocity on the hypersurface
    Vector3D* mFluidUx;
    Vector3D* mFluidUy;
    Vector3D* mFluidRapidity;
};

#endif

/*! @file Hypersurface_Lhyquid3D.h
 * @brief Definition of Hypersurface_Lhyquid3D class. 3+1 freeze-out hypersurface.
 */
/*! @class Hypersurface_Lhyquid3D
 * @brief 3+1 freeze-out hypersurface.
 *
 * Output form the <b>Lhyquid</b> hydrodynamic code.<br />
 * Point on the hypersurface is given by a four-vector originating from the space-time point
 * @f[ x_1^\mu = (t_1,x_1,y_1,z_1) = (\tau_i, 0, 0, 0) @f]
 * The four-vector's end point @f$ x_2^\mu = (t_2,x_2,y_2,z_2) @f$ is defined by function @f$ d @f$ 
 * parametrized by three angles:
 * <ul>
 *   <li> @f$ \zeta \in <0, \pi/2> @f$ - angle in the @f$ \rho-\tau @f$ plane (Hypersurface_Lhyquid3D::Zeta)
 *   <li> @f$ \phi \in <0, 2\pi> @f$ - azimuthal angle in @f$ x-y @f$ plane (Hypersurface_Lhyquid3D::PhiS)
 *   <li> @f$ \theta \in <0, \pi> @f$ - angle in the @f$ \Lambda \eta_\parallel-xy\tau @f$ plane (Hypersurface_Lhyquid3D::Theta)
 * </ul>
 * The space-time coordinates are calculated in the following way:
 * @f{eqnarray}{
 *   t_2 &=& \left[ \tau_i + d(\zeta,\phi) \sin(\theta) \sin(\zeta) \right] \cosh(d(\zeta,\phi) \cos(\theta) / \Lambda) \nonumber \\
 *   x_2 &=& d(\zeta,\phi) \sin(\theta) \cos(\zeta) \cos(\phi) \nonumber \\
 *   y_2 &=& d(\zeta,\phi) \sin(\theta) \cos(\zeta) \sin(\phi) \nonumber \\
 *   z_2 &=& \left[ \tau_i + d(\zeta,\phi) \sin(\theta) \sin(\zeta) \right] \sinh(d(\zeta,\phi) \cos(\theta) / \Lambda) \nonumber
 * @f}
 * Other useful variables:
 * <ul>
 *   <li> @f$ \rho = d(\zeta,\phi) \sin(\theta) \cos(\zeta) @f$
 *   <li> @f$ \tau = \tau_i + d(\zeta,\phi) \sin(\theta) \sin(\zeta) @f$
 *   <li> @f$ \eta_\parallel = d(\zeta,\phi) \cos(\theta) / \Lambda @f$
 * </ul>
 *
 * @fn Hypersurface_Lhyquid3D::Hypersurface_Lhyquid3D()
 * @brief Default constructor.
 *
 * @fn Hypersurface_Lhyquid3D::~Hypersurface_Lhyquid3D()
 * @brief Destructor.
 *
 * @fn double Hypersurface_Lhyquid3D::GetHyperCubeSpatialVolume()
 * @brief Returns the spatial-type volume of the Monte-Carlo hypercube @f$ V_{hc} = \eta_\parallel^{rng} (\zeta_{max}-\zeta_{min}) (\phi_{max} - \phi_{min})@f$
 *
 * @fn double Hypersurface_Lhyquid3D::GetDSigmaP(double aMt, double aPt, double aPhiP, double aRapP);
 * @brief Return @f$ d\Sigma_\mu p^\mu@f$ for given four-momentum of the particle.
 * @param [in] aMt particle transverse mass
 * @param [in] aPt particle transverse momentum
 * @param [in] aPhiP particle azimuthal angle of the momentum 
 * @param [in] aRapP particle rapidity
 *
 * @fn double Hypersurface_Lhyquid3D::GetPdotU(  double aMt, double aPt, double aPhiP, double aRapP);
 * @brief Return @f$ p_\mu u^\mu@f$ for given four-momentum of the particle.
 * @param [in] aMt particle transverse mass
 * @param [in] aPt particle transverse momentum
 * @param [in] aPhiP particle azimuthal angle of the momentum 
 * @param [in] aRapP particle rapidity
 *
 * @fn void Hypersurface_Lhyquid3D::SetPositionOnHypersurface(TRandom2* aRandom, double* aXt, double* aXx, double* aXy, double* aXz);
 * @brief Sets a random position on the hypersurface and returns the space-time coordinates.
 * @param [in] aRandom pointer to <a href="http://root.cern.ch/root/html/TRandom2.html">TRandom2</a> random number generator
 * @param [out] aXt particle time
 * @param [out] aXx particle x-coordinate
 * @param [out] aXy particle y-coordinate
 * @param [out] aXz particle z-coordinate
 *
 * @var double Hypersurface_Lhyquid3D::mRapPRange
 * @brief @f$ Y^{rng} @f$, rapidity range
 *
 * @var double Hypersurface_Lhyquid3D::mLambda
 * @brief @f$ \Lambda \,[GeV^{-1}]@f$, conversion constant to spatial rapidity @f$ \eta_\parallel^{max} = d(\zeta,\phi,\theta = 0) / \Lambda @f$
 *
 * @var double Hypersurface_Lhyquid3D::mTauI
 * @brief @f$ \tau_i \,[GeV^{-1}]@f$, initial proper time of the hydrodynamic evolution
 *
 * @var double Hypersurface_Lhyquid3D::mCollidingEnergy
 * @brief (optional) @f$ \sqrt{s_{NN}} \,[GeV] @f$, center of mass energy
 *
 * @var double Hypersurface_Lhyquid3D::mCentralityMin
 * @brief (optional) @f$ c_{min} \,[\%] @f$, centrality lower limit
 *
 * @var double Hypersurface_Lhyquid3D::mCentralityMax
 * @brief (optional) @f$ c_{max} \,[\%] @f$, centrality higher limit
 *
 * @var double Hypersurface_Lhyquid3D::mImpactParameter
 * @brief (optional) @f$ b \,[fm] @f$, impact parameter
 *
 * @var double Hypersurface_Lhyquid3D::mTempI
 * @brief (optional) @f$ T_i \,[MeV] @f$, temperature in the center of the fireball at time @f$ \tau_i @f$
 *
 * @var char Hypersurface_Lhyquid3D::mDeviceName[30]
 * @brief (optional) name of the experimental device e.g. SPS, RHIC, LHC, etc.
 *
 * @var char Hypersurface_Lhyquid3D::mCollidingSystem[30]
 * @brief (optional) colliding particles e.g. pp, dCu, AuAu, PbPb, etc.
 * 
 * @var Thermodynamics* Hypersurface_Lhyquid3D::mThermo
 * @brief Pointer to the Thermodynamics class.
 */
