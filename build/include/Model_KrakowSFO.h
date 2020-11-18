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

#ifndef _TH2_MODEL_CRACOW_SINGLE_FREEZE_OUT_H_
  #define _TH2_MODEL_CRACOW_SINGLE_FREEZE_OUT_H_

#include "Model.h"
#include "Thermodynamics.h"

class Model_KrakowSFO
: public Model
{
  public:
    Model_KrakowSFO();
    Model_KrakowSFO(TRandom2* aRandom);
    ~Model_KrakowSFO();
    
    double GetIntegrand(ParticleType* aPartType);
    void   AddParameterBranch(TTree* aTree);
    
  protected:
    double mRapPRange;
    double mRapSRange;
    double mTauC;
    double mRhoMax;
    Thermodynamics* mThermo;
        
  private:
    void Description();
    void ReadParameters();
};

#endif

/*! @file Model_KrakowSFO.h
 * @brief Definition of Model_KrakowSFO class. Krakow Single Freeze-out Model.
 */
/*! @class Model_KrakowSFO
 * @brief Krakow Single Freeze-out model.
 * 
 * <table>
 *   <tr><th>sModel</th>	<th>Model::mModelName</th>	<th>Description</th></tr>
 *   <tr><td>0</td>		<td>"KrakowSFO"</td>		<td>Krakow Single Freeze-out model</td></tr>
 * </table>
 * 
 * The integrand has the form
 * @f{eqnarray}{
 *   I			&=& \frac{2 s+1}{\left(2\pi\right)^3} \,d\Sigma_\mu p^\mu \left[ \exp \left( \frac{p_\mu u^\mu - \mu_{ch}}{T}\right) \pm 1\right]^{-1} p_T \,dp_T \,d\phi_p \,dy \nonumber \\
 *   \nonumber\\ \nonumber\\
 *   d\Sigma_\mu p^\mu	&=& \left[ \sqrt{\tau_C^2 + \rho^2} \cosh(\eta_\parallel - y) - \rho \,p_T \cos(\phi - \phi_p)\right] \,\rho \,d\rho \,d\phi \,d\eta_\parallel \nonumber \nonumber \\
 *   p_\mu u^\mu	&=& \frac{1}{\tau_C} \left[ \sqrt{\tau_C^2 + \rho^2} \cosh(\eta_\parallel - y) - \rho \,p_T \cos(\phi - \phi_p)\right] \nonumber
 * @f}
 * where (link to all variable tables)
 *
 * @fn Model_KrakowSFO::Model_KrakowSFO()
 * @brief Default constructor.
 *
 * @fn Model_KrakowSFO::Model_KrakowSFO(TRandom2* aRandom)
 * @brief Constructor that links a random number generator form Integrator.
 * @param [in] aRandom pointer to <a href="http://root.cern.ch/root/html/TRandom2.html">TRandom2</a> random number generator
 *
 * @fn Model_KrakowSFO::~Model_KrakowSFO()
 * @brief Destructor.
 *
 * @fn double Model_KrakowSFO::GetIntegrand(ParticleType* aPartType)
 * @brief Generates random space-time coordinates and four-momentum of a particle of a given type 
 * and returns the value of the integrand of the Cooper-Frye formula.
 * @param [in] aPartType pointer to the ParticleType class holding particle properties.
 *
 * @fn void Model_KrakowSFO::AddParameterBranch(TTree* aTree)
 * @brief Function called by EventGenerator to add a new <a href="http://root.cern.ch/root/html/TBranch.html">TBranch</a> 
 * structure to an opened event <a href="http://root.cern.ch/root/html/TFile.html">TFile</a>.
 * 
 * This branch contains model parameters and model description.
 * @param [in] aTree pointer to the current <a href="http://root.cern.ch/root/html/TTree.html">TTree</a> object.
 * 
 * @var double Model_KrakowSFO::mRapPRange
 * @brief @f$ Y^{rng} @f$, rapidity range
 *
 * @var double Model_KrakowSFO::mRapSRange
 * @brief @f$ \eta_\parallel^{rng} @f$, spatial rapidity range
 *
 * @var double Model_KrakowSFO::mRhoMax
 * @brief @f$ \rho_{max} \,[GeV^{-1}] @f$, size of the fire-cylinder
 *
 * @var double Model_KrakowSFO::mTauC
 * @brief @f$ \tau_C \sqrt{t^2 - x^2 - y^2 - z^2} \,[GeV^{-1}] @f$, Krakow-type proper time
 * 
 * The well known proper time can be evaluated as @f$ \tau = \sqrt{\tau_C^2 + \rho^2} @f$.
 * 
 * @var Thermodynamics* Model_KrakowSFO::mThermo
 * @brief Pointer to the Thermodynamics class.
 */
