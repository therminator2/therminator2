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

#ifndef _TH2_MODEL_LHYQUID2DBI_H_
  #define _TH2_MODEL_LHYQUID2DBI_H_

#include "Model.h"
#include "Hypersurface_Lhyquid2D.h"

class Model_Lhyquid2DBI
: public Model,
  public Hypersurface_Lhyquid2D
{
  public:
    Model_Lhyquid2DBI();
    Model_Lhyquid2DBI(TRandom2* aRandom);
    ~Model_Lhyquid2DBI();
    
<<<<<<< HEAD
    double GetIntegrand(ParticleType* aPartType, bool);
=======
    double GetIntegrand(ParticleType* aPartType);
>>>>>>> 6eef3d1fc1f3ac17aedbca11e4a34ae9c4017af9
    void   AddParameterBranch(TTree* aTree);
    
  private:
    void Description();
    void ReadParameters();
};

#endif

/*! @file Model_Lhyquid2DBI.h
 * @brief Definition of Model_Lhyquid2DBI class. Hydro based 2+1 freeze-out Model with boost-invariance.
 */
/*! @class Model_Lhyquid2DBI
 * @brief Lhyquid2DBI is a boost-invariant freeze-out model that uses a hypersurface from 
 * hydrodynamic calculation as input.
 *
 * <table>
 *   <tr><th>sModel</th>	<th>Model::mModelName</th>	<th>Description</th></tr>
 *   <tr><td>11</td>		<td>"Lhyquid2DBI"</td>		<td>Hydro based freeze-out model with boost-invariance</td></tr>
 * </table>
 *
 * @fn Model_Lhyquid2DBI::Model_Lhyquid2DBI()
 * @brief Default constructor.
 *
 * @fn Model_Lhyquid2DBI::Model_Lhyquid2DBI(TRandom2* aRandom)
 * @brief Constructor that links a random number generator form Integrator.
 * @param [in] aRandom pointer to <a href="http://root.cern.ch/root/html/TRandom2.html">TRandom2</a> random number generator
 *
 * @fn Model_Lhyquid2DBI::~Model_Lhyquid2DBI()
 * @brief Destructor.
 *
 * @fn double Model_Lhyquid2DBI::GetIntegrand(ParticleType* aPartType)
 * @brief Generates random space-time coordinates and four-momentum of a particle of a given type 
 * and returns the value of the integrand of the Cooper-Frye formula.
 * @param [in] aPartType pointer to the ParticleType class holding particle properties.
 *
 * @fn void Model_Lhyquid2DBI::AddParameterBranch(TTree* aTree)
 * @brief Function called by EventGenerator to add a new <a href="http://root.cern.ch/root/html/TBranch.html">TBranch</a> 
 * structure to an opened event <a href="http://root.cern.ch/root/html/TFile.html">TFile</a>.
 * 
 * This branch contains model parameters and model description.
 * @param [in] aTree pointer to the current <a href="http://root.cern.ch/root/html/TTree.html">TTree</a> object.
 * 
 */
