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

#ifndef _TH2_PARTICLE_COORDINATES_H_
  #define _TH2_PARTICLE_COORDINATES_H_

#include <Rtypes.h>

#define _PARTICLES_TREE_  "particles"
#define _PARTICLE_BRANCH_ "particle"
#define _PARTICLE_FORMAT_ "mass/F:t:x:y:z:e:px:py:pz:decayed/I:pid:fatherpid:rootpid:eid:fathereid:eventid/i"

class ParticleCoor {
  public:
    ParticleCoor()
    : mass(0.0), t(0.0), x(0.0), y(0.0), z(0.0), e(0.0), px(0.0), py(0.0), pz(0.0),
      decayed(0), pid(0), fatherpid(0), rootpid(0), eid(0), fathereid(0), eventid(0)
    { }
    ~ParticleCoor()
    { }

    ParticleCoor& operator= (ParticleCoor& aPartCoor);
    ParticleCoor& operator= (ParticleCoor* aPartCoor);

    void    GetMomentum(double* aPe, double* aPx, double* aPy, double* aPz) const;
    void    GetPosition(double* aXt, double* aXx, double* aXy, double* aXz) const;

    Float_t GetTau() const;
    Float_t GetR() const;
    Float_t GetRho() const;
    Float_t GetPhiS() const;
    Float_t GetRapidityS() const;

    Float_t GetP() const;
    Float_t GetPt() const;
    Float_t GetMt() const;
    Float_t GetPhiP() const;
    Float_t GetRapidityP() const;
    Float_t GetEtaP() const;

    Int_t   GetDecayed() const;
    void    SetDecayed();

    void    TransformToLCMS(Float_t aBetaZ);
    void    TransformRotateZ(Float_t aPhi);
    void    TransformToPRF(Float_t aBetaT);

    Float_t mass;
    Float_t t, x,  y,  z;
    Float_t e, px, py, pz;
    Int_t   decayed;
    Int_t   pid;
    Int_t   fatherpid;
    Int_t   rootpid;
    Int_t   eid;
    Int_t   fathereid;
    UInt_t  eventid;
};

#endif

/*! @file ParticleCoor.h
 * @brief Definition of Particle class. %Particle entry to the <b>eventXXX.root</b> file.
 * 
 * Information copied form Particle object.
 */
/*! @def _PARTICLES_TREE_
 * @brief Macro with the name of the <a href="http://root.cern.ch/root/html/TTree.html">TTree</a> object.
 *
 * @def _PARTICLE_BRANCH_
 * @brief Macro with the name of a <a href="http://root.cern.ch/root/html/TBranch.html">TBranch</a> in the <a href="http://root.cern.ch/root/html/TTree.html">TTree</a>. 
 *
 * @def _PARTICLE_FORMAT_
 * @brief Macro with the format in which the information is stored in the ROOT event file.
 */
/*! @class ParticleCoor
 * @brief Class with particle coordinates.
 *
 * @fn ParticleCoor::ParticleCoor()
 * @brief Default constructor.
 * 
 * @fn ParticleCoor::~ParticleCoor()
 * @brief Destructor.
 * 
 * @fn ParticleCoor& ParticleCoor::operator= (ParticleCoor& aPartCoor)
 * @brief Assignment operator.
 * @param [in,out] aPartCoor reference to ParticleCoor
 * @retval ParticleCoor& reference to ParticleCoor
 * 
 * @fn ParticleCoor& ParticleCoor::operator= (ParticleCoor* aPartCoor)
 * @brief Assignment operator.
 * @param [in,out] aPartCoor pointer to ParticleCoor
 * @retval ParticleCoor& reference to ParticleCoor
 * 
 * @fn void ParticleCoor::GetMomentum(double* aPe, double* aPx, double* aPy, double* aPz) const
 * @brief Returns particle energy and momentum.
 * @param [out] aPe particle energy
 * @param [out] aPx particle px-coordinate
 * @param [out] aPy particle py-coordinate
 * @param [out] aPz particle pz-coordinate
 * 
 * @fn void ParticleCoor::GetPosition(double* aXt, double* aXx, double* aXy, double* aXz) const
 * @brief Returns particle space-time position.
 * @param [out] aXt particle time
 * @param [out] aXx particle x-coordinate
 * @param [out] aXy particle y-coordinate
 * @param [out] aXz particle z-coordinate
 * 
 * @fn Float_t ParticleCoor::GetTau() const
 * @brief Returns particle proper time.
 *
 * @fn Float_t ParticleCoor::GetR() const;
 * @brief Returns particle distance (spherical coordinates)
 * 
 * @fn Float_t ParticleCoor::GetRho() const
 * @brief Returns particle transverse distance.
 *
 * @fn Float_t ParticleCoor::GetPhiS() const
 * @brief Returns particle azimuthal angle.
 *
 * @fn Float_t ParticleCoor::GetRapidityS() const
 * @brief Returns particle spatial rapidity
 *
 * @fn Float_t ParticleCoor::GetP() const;
 * @brief Returns particle momentum.
 * 
 * @fn Float_t ParticleCoor::GetPt() const
 * @brief Returns particle transverse momentum.
 *
 * @fn Float_t ParticleCoor::GetMt() const
 * @brief Returns particle transverse mass.
 *
 * @fn Float_t ParticleCoor::GetPhiP() const
 * @brief Returns particle transverse momentum azimuthal angle.
 *
 * @fn Float_t ParticleCoor::GetRapidityP() const
 * @brief Returns particle rapidity.
 *
 * @fn Float_t ParticleCoor::GetEtaP() const
 * @brief Returns particle pseudo-rapidity.
 *
 * @fn Int_t ParticleCoor::GetDecayed() const
 * @brief If particle has decayed.
 * @retval 0 particle did not decay
 * @retval 1 particle decayed
 *
 * @fn void ParticleCoor::SetDecayed()
 * @brief Particle has decayed.
 * 
 * @fn void ParticleCoor::TransformToLCMS(Float_t aBetaZ)
 * @brief Transformes the particle spatial and momentum coordinates to LCMS.
 * @param [in] aBetaZ longnitudal velocity (@f$ \beta_z = \frac{v_z}{c} @f$)
 * 
 * @fn void ParticleCoor::TransformRotateZ(Float_t aPhi)
 * @brief Rotates the particle spatial and momentum coordinates in the transverse plane.
 * @param [in] aPhi azimuthal angle (@f$ \phi @f$)
 * 
 * @fn void ParticleCoor::TransformToPRF(Float_t aBetaT)
 * @brief Transformes the particle spatial and momentum coordinates to PRF.
 * @param [in] aBetaT transverse velocity (@f$ \beta_T = \frac{v_T}{c} @f$)
 * 
 * @var Float_t ParticleCoor::mass
 * @brief @f$ m @f$ mass
 *
 * @var Float_t ParticleCoor::t
 * @brief @f$ t @f$ time coordinate
 *
 * @var Float_t ParticleCoor::x
 * @brief @f$ x @f$ coordinate
 *
 * @var Float_t ParticleCoor::y
 * @brief @f$ y @f$ coordinate
 *
 * @var Float_t ParticleCoor::z
 * @brief @f$ z @f$ coordinate
 *
 * @var Float_t ParticleCoor::e
 * @brief @f$ E @f$ energy
 *
 * @var Float_t ParticleCoor::px
 * @brief @f$ p_x @f$ coordinate
 *
 * @var Float_t ParticleCoor::py
 * @brief @f$ p_y @f$ coordinate
 *
 * @var Float_t ParticleCoor::pz
 * @brief @f$ p_z @f$ coordinate
 *
 * @var Int_t ParticleCoor::decayed
 * @brief information if particle has decayed
 *
 * @var Int_t ParticleCoor::pid
 * @brief particle ID, form PDG
 *
 * @var Int_t ParticleCoor::fatherpid
 * @brief father particle ID, form PDG
 *
 * @var Int_t ParticleCoor::rootpid
 * @brief root (primordial) particle ID, form PDG
 *
 * @var Int_t ParticleCoor::eid
 * @brief event unique particle ID
 *
 * @var Int_t ParticleCoor::fathereid
 * @brief event unique father particle ID
 *
 * @var UInt_t ParticleCoor::eventid
 * @brief event number
 */
