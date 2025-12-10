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

#ifndef _TH2_MODEL_EXAMPLE_H_
  #define _TH2_MODEL_EXAMPLE_H_

#include "Model.h"
#include "Thermodynamics.h"

class Model_Example
: public Model
{
  public:
    Model_Example();
    Model_Example(TRandom2* aRandom);
    ~Model_Example();
    
    std::pair<double, double> GetIntegrand(ParticleType* aPartType, bool finiteWidth, bool positionDistribution = false);
    void   AddParameterBranch(TTree* aTree);
    
  protected:
    Thermodynamics* mThermo;
    double mParameter1;	// model parameter #1
    double mParameter2;	// model parameter #2
    double mParameter3;	// model parameter #3
    
  private:
    void Description();
    void ReadParameters();
};

#endif

/*! @file Model_Example.h
 * @brief Definition of Model_Example class. Example freeze-out Model.
 */
/*! @class Model_Example
 * @brief The purpose of this class is to show how to build custom model class.
 * 
 * <table>
 *   <tr><th>sModel</th>	<th>Model::mModelName</th>	<th>Description</th></tr>
 *   <tr><td>99</td>		<td>"Example"</td>		<td>Example model class</td></tr>
 * </table>
 *
 * @fn Model_Example::Model_Example()
 * @brief Default constructor.
 *
 * @fn Model_Example::Model_Example(TRandom2* aRandom)
 * @brief Constructor that links a random number generator form Integrator.
 * @param [in] aRandom pointer to <a href="http://root.cern.ch/root/html/TRandom2.html">TRandom2</a> random number generator
 *
 * @fn Model_Example::~Model_Example()
 * @brief Destructor.
 *
 * @fn double Model_Example::GetIntegrand(ParticleType* aPartType)
 * @brief Generates random space-time coordinates and four-momentum of a particle of a given type 
 * and returns the value of the integrand of the Cooper-Frye formula.
 * @param [in] aPartType pointer to the ParticleType class holding particle properties.
 *
 * @fn void Model_Example::AddParameterBranch(TTree* aTree)
 * @brief Function called by EventGenerator to add a new <a href="http://root.cern.ch/root/html/TBranch.html">TBranch</a> 
 * structure to an opened event <a href="http://root.cern.ch/root/html/TFile.html">TFile</a>.
 * 
 * This branch contains model parameters and model description.
 * @param [in] aTree pointer to the current <a href="http://root.cern.ch/root/html/TTree.html">TTree</a> object.
 * 
 * @var double Model_Example::mParameter1
 * @brief model parameter 1
 *
 * @var double Model_Example::mParameter2
 * @brief model parameter 2
 *
 * @var double Model_Example::mParameter3
 * @brief model parameter 3
 * 
 * @var Thermodynamics* Model_Example::mThermo
 * @brief Pointer to the Thermodynamics class.
 */
