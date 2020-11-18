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

#ifndef _TH2_CHEMISTRY_H_
  #define _TH2_CHEMISTRY_H_

#include "Vector3D.h"
#include "ParticleType.h"

class Chemistry
{
  public:
    Chemistry();
    ~Chemistry();

    int    GetChemistryType() const;  
    double GetChemicalPotential(ParticleType* aPartType) const;
    double GetChemicalPotential(ParticleType* aPartType, double aX, double aY, double aZ) const;
    double GetFugacity(ParticleType* aPartType) const;
    double GetFugacity(ParticleType* aPartType, double aX, double aY, double aZ) const;
    double GetMuB() const;
    double GetMuI() const;
    double GetMuS() const;
    double GetMuC() const;
    double GetLambdaQ() const;
    double GetLambdaI() const;
    double GetLambdaS() const;
    double GetLambdaC() const;
    double GetGammaQ() const;
    double GetGammaS() const;
    double GetGammaC() const;
    void   SetChemistry(double    aMuB,     double    aMuI,     double    aMuS,     double    aMuC);
    void   SetChemistry(Vector3D* aMuB,     Vector3D* aMuI,     Vector3D* aMuS,     Vector3D* aMuC);
    void   SetChemistry(double    aLambdaQ, double    aLambdaI, double    aLambdaS, double    aLambdaC,
			double    aGammaQ,                      double    aGammaS,  double    aGammaC);
    void   SetChemistry(Vector3D* aLambdaQ, Vector3D* aLambdaI, Vector3D* aLambdaS, Vector3D* aLambdaC,
			Vector3D* aGammaQ,                      Vector3D* aGammaS,  Vector3D* aGammaC);  
    void   SetGammas(double    aGammaQ, double aGammaS, double aGammaC);

  private:
    int		mChemistryType;

    double	mMuBConst;
    double	mMuIConst;
    double	mMuSConst;
    double	mMuCConst;
    Vector3D*	mMuBVar;
    Vector3D*	mMuIVar;
    Vector3D*	mMuSVar;
    Vector3D*	mMuCVar;

    double	mGammaQConst;
    double	mGammaSConst;
    double	mGammaCConst;
    double	mLambdaIConst;
    double	mLambdaQConst;
    double	mLambdaSConst;
    double	mLambdaCConst;
    Vector3D*	mGammaQVar;
    Vector3D*	mGammaSVar;
    Vector3D*	mGammaCVar;
    Vector3D*	mLambdaIVar;
    Vector3D*	mLambdaQVar;
    Vector3D*	mLambdaSVar;
    Vector3D*	mLambdaCVar;
};

#endif

/*! @file Chemistry.h
 * @brief Definition of Chemistry class. Gives the local chemical potential or fugacity.
 */
/*! @class Chemistry
 * @brief Class gives the local chemical potential or fugacity for the given ParticleType.
 *
 * Depending on the number and type of arguments given to the SetChemistry() function, the class 
 * will return chemical potential or fugacity that is constant or as a function of <b>x</b>,<b>y</b> and <b>z</b>.
 * 
 * ParticleType class holds information on the particle's quark content @f$ N_i @f$ (with index 
 * @f$ i = q, \bar q, s, \bar s, c, \bar c @f$) and quantum numbers @f$ B @f$, @f$ I_3 @f$, @f$ S @f$ and @f$ C @f$. <br />
 * Chemical potential is calculated as follows
 * @f[ \mu = \mu_B B + \mu_{I_3} I_3 + \mu_S S + \mu_C C@f]
 * and the fugacity
 * @f[ \Upsilon = \lambda_{I_3} \,\gamma_q^{N_q + N_{\bar q}} \,\lambda_q^{N_q - N_{\bar q}} \,\gamma_s^{N_s + N_{\bar s}} \,\lambda_s^{N_s - N_{\bar s}} \,\gamma_c^{N_c + N_{\bar c}} \,\lambda_c^{N_c - N_{\bar c}} @f]
 *
 * @fn Chemistry::Chemistry()
 * @brief Default constructor.
 *
 * @fn Chemistry::~Chemistry()
 * @brief Destructor.
 *
 * @fn int Chemistry::GetChemistryType() const
 * @brief Get the type of chemistry stored.
 * @retval 0 constant chemical potentials @f$ \mu_i @f$
 * @retval 1 chemical potentials @f$ \mu_i(x, y, z) @f$
 * @retval 2 constant fugacities @f$ \lambda_j @f$ and phase space occupancies @f$ \gamma_k @f$
 * @retval 3 fugacities @f$ \lambda_j(x, y, z) @f$ and phase space occupancies @f$ \gamma_k(x, y, z) @f$
 * @retval -1 unknown type
 * 
 * where @f$ i = B, I_3, S, C @f$, @f$ j = I_3, q, s, c @f$, @f$ k = q, s, c @f$.
 *
 * @fn double Chemistry::GetChemicalPotential(ParticleType* aPartType) const
 * @brief Returns the constant chemical potential @f$ \mu @f$ for a given particle type.
 * @param [in] aPartType type of the particle
 *
 * @fn double Chemistry::GetChemicalPotential(ParticleType* aPartType, double aX, double aY, double aZ) const
 * @brief Returns the chemical potential @f$ \mu(x, y, z) @f$ for a given particle type
 * @param [in] aPartType type of the particle
 * @param [in] aX x-coordinate
 * @param [in] aY y-coordinate
 * @param [in] aZ z-coordinate
 *
 * @fn double Chemistry::GetFugacity(ParticleType* aPartType) const
 * @brief Returns the constant fugacity @f$ \Upsilon @f$ for a given particle type.
 * @param [in] aPartType type of the particle
 *
 * @fn double Chemistry::GetFugacity(ParticleType* aPartType, double aX, double aY, double aZ) const
 * @brief Returns the fugacity @f$ \Upsilon(x, y, z) @f$ for a given particle type
 * @param [in] aPartType type of the particle
 * @param [in] aX x-coordinate
 * @param [in] aY y-coordinate
 * @param [in] aZ z-coordinate
 *
 * @fn double Chemistry::GetMuB() const
 * @brief returns baryochemical potential
 * @warning only constant type
 *
 * @fn double Chemistry::GetMuI() const
 * @brief returns isospin chemical potential
 * @warning only constant type
 *
 * @fn double Chemistry::GetMuS() const
 * @brief returns strange chemical potential
 * @warning only constant type
 *
 * @fn double Chemistry::GetMuC() const
 * @brief returns charm chemical potential
 * @warning only constant type
 *
 * @fn double Chemistry::GetLambdaQ() const
 * @brief returns light quark fugacity factor
 * @warning only constant type
 *
 * @fn double Chemistry::GetLambdaI() const
 * @brief returns isospin fugacity factor
 * @warning only constant type
 *
 * @fn double Chemistry::GetLambdaS() const
 * @brief returns strange quark fugacity factor
 * @warning only constant type
 *
 * @fn double Chemistry::GetLambdaC() const
 * @brief returns charm quark fugacity factor
 * @warning only constant type
 *
 * @fn double Chemistry::GetGammaQ() const
 * @brief returns light quark phase space occupancy
 * @warning only constant type
 *
 * @fn double Chemistry::GetGammaS() const
 * @brief returns strange quark phase space occupancy
 * @warning only constant type
 *
 * @fn double Chemistry::GetGammaC() const
 * @brief returns charm quark phase space occupancy
 * @warning only constant type
 *
 * @fn void Chemistry::SetChemistry(double aMuB, double aMuI, double aMuS, double aMuC)
 * @brief Sets the chemistry type to constant chemical potentials.
 * @param [in] aMuB baryochemical potential
 * @param [in] aMuI isospin chemical potential
 * @param [in] aMuS strange chemical potential
 * @param [in] aMuC strange chemical potential
 *
 * @fn void Chemistry::SetChemistry(Vector3D* aMuB, Vector3D* aMuI, Vector3D* aMuS, Vector3D* aMuC)
 * @brief Sets the chemistry type to variable chemical potentials.
 * @param [in] aMuB baryochemical potential
 * @param [in] aMuI isospin chemical potential
 * @param [in] aMuS strange chemical potential
 * @param [in] aMuC strange chemical potential
 *
 * @fn void Chemistry::SetChemistry(double aLambdaQ, double aLambdaI, double aLambdaS, double aLambdaC, double aGammaQ, double aGammaS, double aGammaC)
 * @brief Sets the chemistry type to constant fugacities and phase space occupancies.
 * @param [in] aLambdaQ light quark fugacity factor
 * @param [in] aLambdaI isospin fugacity factor
 * @param [in] aLambdaS strange quark fugacity factor
 * @param [in] aLambdaC charm quark fugacity factor
 * @param [in] aGammaQ light quark phase space occupancy
 * @param [in] aGammaS strange quark phase space occupancy
 * @param [in] aGammaC charm quark phase space occupancy
 * 
 *
 * @fn void Chemistry::SetChemistry(Vector3D* aLambdaQ, Vector3D* aLambdaI, Vector3D* aLambdaS, Vector3D* aLambdaC, Vector3D* aGammaQ, Vector3D* aGammaS, Vector3D* aGammaC)
 * @brief Sets the chemistry type to variable fugacities and phase space occupancies.
 * @param [in] aLambdaQ light quark fugacity factor
 * @param [in] aLambdaI isospin fugacity factor
 * @param [in] aLambdaS strange quark fugacity factor
 * @param [in] aLambdaC charm quark fugacity factor
 * @param [in] aGammaQ light quark phase space occupancy
 * @param [in] aGammaS strange quark phase space occupancy
 * @param [in] aGammaC charm quark phase space occupancy
 */
