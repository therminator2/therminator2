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

#ifndef _TH2_MODEL_BLASTWAVE_H_
  #define _TH2_MODEL_BLASTWAVE_H_

#include "Model.h"
#include "Thermodynamics.h"

class Model_BlastWave
: public Model
{
  public:
    Model_BlastWave();
    Model_BlastWave(TRandom2* aRandom);
    ~Model_BlastWave();
    
    double GetIntegrand(ParticleType* aPartType);
    void   AddParameterBranch(TTree* aTree);
    
  protected:
    double mRapPRange;
    double mRapSRange;
    double mRhoMax;
    double mTau;
    double mVt;
    Thermodynamics* mThermo;
    
  private:
    void   Description();
    void   ReadParameters();
};

#endif

/*! @file Model_BlastWave.h
 * @brief Definition of Model_BlastWave class. Blast-Wave Model.
 */
/*! @class Model_BlastWave
 * @brief Blast-Wave freeze-out model.
 * 
 * <table>
 *   <tr><th>sModel</th>	<th>Model::mModelName</th>	<th>Description</th></tr>
 *   <tr><td>1</td>		<td>"BlastWave"</td>		<td>Blast-Wave</td></tr>
 * </table>
 * 
 * The integrand has the form
 * @f{eqnarray}{
 *   I			&=& \frac{2 s+1}{\left(2\pi\right)^3} \,d\Sigma_\mu p^\mu \left[ \exp \left( \frac{p_\mu u^\mu - \mu_{ch}}{T}\right) \pm 1\right]^{-1} p_T \,dp_T \,d\phi_p \,dy \nonumber \\
 *   \nonumber\\ \nonumber\\
 *   d\Sigma_\mu p^\mu	&=& m_T \cosh(\eta_\parallel - y) \,\tau \,\rho \,d\rho \,d\phi \,d\eta_\parallel \nonumber \nonumber \\
 *   p_\mu u^\mu	&=& \frac{1}{\sqrt{1-v_T^2}} \left[ m_T \cosh(\eta_\parallel - y) - v_T \,p_T \cos(\phi - \phi_p)\right] \nonumber
 * @f}
 * where (link to all variable tables)
 *
 * @fn Model_BlastWave::Model_BlastWave()
 * @brief Default constructor.
 *
 * @fn Model_BlastWave::Model_BlastWave(TRandom2* aRandom)
 * @brief Constructor that links a random number generator form Integrator.
 * @param [in] aRandom pointer to <a href="http://root.cern.ch/root/html/TRandom2.html">TRandom2</a> random number generator
 *
 * @fn Model_BlastWave::~Model_BlastWave()
 * @brief Destructor.
 *
 * @fn double Model_BlastWave::GetIntegrand(ParticleType* aPartType)
 * @brief Generates random space-time coordinates and four-momentum of a particle of a given type 
 * and returns the value of the integrand of the Cooper-Frye formula.
 * @param [in] aPartType pointer to the ParticleType class holding particle properties.
 *
 * @fn void Model_BlastWave::AddParameterBranch(TTree* aTree)
 * @brief Function called by EventGenerator to add a new <a href="http://root.cern.ch/root/html/TBranch.html">TBranch</a> 
 * structure to an opened event <a href="http://root.cern.ch/root/html/TFile.html">TFile</a>.
 * 
 * This branch contains model parameters and model description.
 * @param [in] aTree pointer to the current <a href="http://root.cern.ch/root/html/TTree.html">TTree</a> object.
 * 
 * @var double Model_BlastWave::mRapPRange
 * @brief @f$ Y^{rng} @f$, rapidity range
 *
 * @var double Model_BlastWave::mRapSRange
 * @brief @f$ \eta_\parallel^{rng} @f$, spatial rapidity range
 *
 * @var double Model_BlastWave::mRhoMax
 * @brief @f$ \rho_{max} \,[GeV^{-1}] @f$, size of the fire-cylinder
 *
 * @var double Model_BlastWave::mTau
 * @brief @f$ \tau \,[GeV^{-1}] @f$, proper time
 *
 * @var double Model_BlastWave::mVt
 * @brief @f$ v_T \,[c] @f$, transverse velocity
 * 
 * @var Thermodynamics* Model_BlastWave::mThermo
 * @brief Pointer to the Thermodynamics class
 */
