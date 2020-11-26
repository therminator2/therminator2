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

#ifndef _TH2_MODEL_LHYQUID3D_H_
  #define _TH2_MODEL_LHYQUID3D_H_

#include "Model.h"
#include "Hypersurface_Lhyquid3D.h"

class Model_Lhyquid3D
: public Model,
  public Hypersurface_Lhyquid3D
{
  public:
    Model_Lhyquid3D();
    Model_Lhyquid3D(TRandom2* aRandom);
    ~Model_Lhyquid3D();
    
    double GetIntegrand(ParticleType* aPartType, bool);
    void   AddParameterBranch(TTree* aTree);
    
  private:
    void Description();
    void ReadParameters();
};

#endif

/*! @file Model_Lhyquid3D.h
 * @brief Definition of Model_Lhyquid3D class. Hydro based 3+1 freeze-out Model.
 */
/*! @class Model_Lhyquid3D
 * @brief Lhyquid3D is a 3+1 freeze-out model that uses a hypersurface from hydrodynamic calculation as input.
 *
 * <table>
 *   <tr><th>sModel</th>	<th>Model::mModelName</th>	<th>Description</th></tr>
 *   <tr><td>10</td>		<td>"Lhyquid3D"</td>		<td>Hydro based 3+1 freeze-out model</td></tr>
 * </table>
 *
 * @fn Model_Lhyquid3D::Model_Lhyquid3D()
 * @brief Default constructor.
 *
 * @fn Model_Lhyquid3D::Model_Lhyquid3D(TRandom2* aRandom)
 * @brief Constructor that links a random number generator form Integrator.
 * @param [in] aRandom pointer to <a href="http://root.cern.ch/root/html/TRandom2.html">TRandom2</a> random number generator
 *
 * @fn Model_Lhyquid3D::~Model_Lhyquid3D()
 * @brief Destructor.
 *
 * @fn double Model_Lhyquid3D::GetIntegrand(ParticleType* aPartType)
 * @brief Generates random space-time coordinates and four-momentum of a particle of a given type 
 * and returns the value of the integrand of the Cooper-Frye formula.
 * @param [in] aPartType pointer to the ParticleType class holding particle properties.
 *
 * @fn void Model_Lhyquid3D::AddParameterBranch(TTree* aTree)
 * @brief Function called by EventGenerator to add a new <a href="http://root.cern.ch/root/html/TBranch.html">TBranch</a> 
 * structure to an opened event <a href="http://root.cern.ch/root/html/TFile.html">TFile</a>.
 * 
 * This branch contains model parameters and model description.
 * @param [in] aTree pointer to the current <a href="http://root.cern.ch/root/html/TTree.html">TTree</a> object.
 * 
 */
