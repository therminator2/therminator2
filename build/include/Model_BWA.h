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

#ifndef _TH2_MODEL_BWA_H_
  #define _TH2_MODEL_BWA_H_

#include "Model.h"
#include "Thermodynamics.h"

class Model_BWA
: public Model
{
  public:
    Model_BWA();
    Model_BWA(TRandom2* aRandom);
    ~Model_BWA();
    
<<<<<<< HEAD
    double GetIntegrand(ParticleType* aPartType, bool);
=======
    double GetIntegrand(ParticleType* aPartType);
>>>>>>> 6eef3d1fc1f3ac17aedbca11e4a34ae9c4017af9
    void   AddParameterBranch(TTree* aTree);

  protected:
    double mRapPRange;
    double mRapSRange;
    double mRhoMax;
    double mTau;
    double mVt;
    double mA;
    double mDelay;
    Thermodynamics* mThermo;
        
  private:
    void   Description();
    void   ReadParameters();
};

#endif

/*! @file Model_BWA.h
 * @brief Definition of Model_BWA class. Modified versions of Blast-Wave Model.
 */
/*! @class Model_BWA
 * @brief Modified Blast-Wave freeze-out models. Parameter <em>Model_BWA::mA</em> tilts the hypersurface in the @f$ \rho - \tau@f$ plane.
 *
 * <table>
 *   <tr><th>sModel</th>	<th>Model::mModelName</th>	<th>Description</th></tr>
 *   <tr><td>2</td>		<td>"BWAVT"</td>		<td>Blast-Wave A-type</td></tr>
 *   <tr><td>3</td>		<td>"BWAVTDelay"</td>		<td>Blast-Wave A-type with delayed emission time</td></tr>
 *   <tr><td>4</td>		<td>"BWAVLinear"</td>		<td>Blast-Wave A-type with linearized transverse flow</td></tr>
 *   <tr><td>5</td>		<td>"BWAVLinearDelay"</td>	<td>Blast-Wave A-type with linearized transverse flow and delayed emission time</td></tr>
 *   <tr><td>6</td>		<td>"BWAVLinearFormation"</td>	<td>Blast-Wave A-type with linearized transverse flow, delayed emission time and space formation</td></tr>
 * </table>
 * <ul>
 *   <li> linearized transverse flow: @f$ v^{lin}_{T}(\rho) = \frac{\rho / \rho_{max}}{v_T + \rho / \rho_{max}} @f$
 *   <li> delayed emission time: @f$ t(X) = t_0 + \lambda * \log(X) @f$
 *        where random variable @f$ X \in (0,1) @f$ and @f$ t_0 @f$ is the particle time on the hypersurface </li>
 *   <li> space formation: @f$ x_{i}(t(X)) = x^{0}_{i} + \frac{p_{i}}{E} \,t @f$ 
 *        where @f$ x^{0}_{i} @f$ - particle position on the hypersurface, 
 *        @f$ p_{i} @f$ - particle momentum and @f$ E @f$ - energy. Index @f$ i @f$ takes values x,y or z.</li>
 * </ul>
 * 
 * The integrand has the form
 * @f{eqnarray}{
 *   I			&=& \frac{2 s+1}{\left(2\pi\right)^3} \,d\Sigma_\mu p^\mu \left[ \exp \left( \frac{p_\mu u^\mu - \mu_{ch}}{T}\right) \pm 1\right]^{-1} p_T \,dp_T \,d\phi_p \,dy \nonumber \\
 *   \nonumber\\ \nonumber\\
 *   d\Sigma_\mu p^\mu	&=& \left[ m_T \cosh(\eta_\parallel - y) - A  \,p_T \cos(\phi - \phi_p) \right] \left( \tau + A \,\rho \right) \,\rho \,d\rho \,d\phi \,d\eta_\parallel \nonumber \nonumber \\
 *   p_\mu u^\mu	&=& \frac{1}{\sqrt{1-v_T^2}} \left[ m_T \cosh(\eta_\parallel - y) - v_T \,p_T \cos(\phi - \phi_p)\right] \nonumber
 * @f}
 * where (link to all variable tables)
 * 
 * Standard Blast-Wave model (Model_BlastWave) is recreated for <em>Model_BWA::mA = 0</em>.<br />
 * Parameter Model_BWA::mDelay in some model versions delays the particle emittion according to 
 * the exponetial distribution.<br />
 *
 * @fn Model_BWA::Model_BWA()
 * @brief Default constructor.
 *
 * @fn Model_BWA::Model_BWA(TRandom2* aRandom)
 * @brief Constructor that links a random number generator form Integrator.
 * @param [in] aRandom pointer to <a href="http://root.cern.ch/root/html/TRandom2.html">TRandom2</a> random number generator
 *
 * @fn Model_BWA::~Model_BWA()
 * @brief Destructor.
 *
 * @fn double Model_BWA::GetIntegrand(ParticleType* aPartType)
 * @brief Generates random space-time coordinates and four-momentum of a particle of a given type 
 * and returns the value of the integrand of the Cooper-Frye formula.
 * @param [in] aPartType pointer to the ParticleType class holding particle properties.
 *
 * @fn void Model_BWA::AddParameterBranch(TTree* aTree)
 * @brief Function called by EventGenerator to add a new <a href="http://root.cern.ch/root/html/TBranch.html">TBranch</a> 
 * structure to an opened event <a href="http://root.cern.ch/root/html/TFile.html">TFile</a>.
 * 
 * This branch contains model parameters and model description.
 * @param [in] aTree pointer to the current <a href="http://root.cern.ch/root/html/TTree.html">TTree</a> object.
 * 
 * @var double Model_BWA::mRapPRange
 * @brief @f$ Y^{rng} @f$, rapidity range
 *
 * @var double Model_BWA::mRapSRange
 * @brief @f$ \eta_\parallel^{rng} @f$, spatial rapidity range
 *
 * @var double Model_BWA::mRhoMax
 * @brief @f$ \rho_{max} \,[GeV^{-1}] @f$, size of the fire-cylinder
 *
 * @var double Model_BWA::mTau
 * @brief @f$ \tau \,[GeV^{-1}] @f$, proper time
 *
 * @var double Model_BWA::mVt
 * @brief @f$ v_T \,[c] @f$, transverse velocity
 *
 * @var double Model_BWA::mA;
 * @brief @f$ A @f$, the A-parameter
 *
 * @var double Model_BWA::mDelay;
 * @brief @f$ \lambda \,[GeV^{-1}]@f$, particle time emission delay
 * 
 * @var Thermodynamics* Model_BWA::mThermo
 * @brief Pointer to the Thermodynamics class
 */
