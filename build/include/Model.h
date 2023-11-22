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

#ifndef _TH2_MODEL_H_
  #define _TH2_MODEL_H_

#include <iomanip>
#include <TTree.h>
#include <TString.h>
#include <TRandom2.h>
#include "ParticleType.h"
#include "Particle.h"
#include "ParticleDB.h"
#include "StructModel.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"

#define MODEL_PAR_DESC(_name_,_unit_)	setw(8)<<_name_<<" "<<setw(7)<<left<<_unit_<<setw(5)<<right<<"#"<<endl
#define MODEL_NAME(_name_)		"# Model ("<<setw(2)<< sModel <<") "<<left<<setw(36)<<_name_<<right<<"#"<<endl

class Model
{
  public:
    Model();
    Model(TRandom2* aRandom);
    virtual ~Model();
    
    virtual double GetIntegrand(ParticleType* aPartType, bool) = 0;
    virtual void   AddParameterBranch(TTree* aTree);
    
    void  	SetParticlePX(Particle* aParticle);
    void SetBreitWignerMap();
    void SetBranchingMap();
    double	GetHyperCubeVolume();
    const char*	GetHash();
    const char*	GetName();
    const char*	GetDescription();
    bool HasBW(TString);
    void GetParticleMass(ParticleType *,bool finiteWidth, double &M, double &spectralFunctionWeight );
    

  protected:
    void	CreateEventSubDir();
    void	CalculateHash(TString aPreHash);

    double CalcMass(ParticleType *, double &statWeigt);
    double CalcMassEBW(ParticleType *, double &statWeigt);
    double CalcMassDistr(ParticleType *, double &statWeigt);


    double	Xt, Xx, Xy, Xz;
    double	Pe, Px, Py, Pz;
    double	mHyperCube;
    TString	mHash;
    TString	mName;
    TString	mDescription;
    TRandom2*	mRandom;

    std::map<int, TF1*> *mSpectralFncs;
    std::map<int, double> *mSpectralFncIntegrals;
    ParticleDB *mDB;
    std::map<TString,TH1D*> *mBreitWigner;
    


};

#endif

/*! @file Model.h
 * @brief Definition of Model class. Provides an interface to derived <em>Model_*</em> classes.
 */
/*! @class Model
 * @brief Base abstract class that provides an interface to user defined derived Model classes.
 * 
 * Class implements a model dependent Cooper-Frye formula. Currently implemented models are:
 * <ul>
 *   <li><b>Model_Cracow</b>		the Cracow Single Freeze-out model
 *   <li><b>Model_BlastWave</b>		classical blast-wave model
 *   <li><b>Model_BWA</b>		the blast-wave model of <em>A-type</em>
 *   <li><b>Model_Lhyquid3D</b>		3-dimensional hydro
 *   <li><b>Model_Lhyquid2DBI</b>	boost-invariant 2-dimensional hydro
 * </ul>
 *
 * @fn Model::Model()
 * @brief Default constructor.
 *
 * @fn Model::Model(TRandom2* aRandom)
 * @brief Constructor that links a random number generator form Integrator.
 * @param [in] aRandom pointer to <a href="http://root.cern.ch/root/html/TRandom2.html">TRandom2</a> random number generator
 *
 * @fn virtual Model::~Model();
 * @brief Destructor.
 *
 * @fn virtual double Model::GetIntegrand(ParticleType* aPartType) = 0
 * @brief Generates random space-time coordinates and four-momentum of a particle of a given type 
 * and returns the value of the integrand of the Cooper-Frye formula.
 * @warning The user MUST define this function in his own derived class.
 * @param [in] aPartType pointer to the ParticleType class holding particle properties.
 *
 * @fn void Model::AddParameterBranch(TTree* aTree)
 * @brief Function called by EventGenerator to add a new <a href="http://root.cern.ch/root/html/TBranch.html">TBranch</a> 
 * structure to an opened event <a href="http://root.cern.ch/root/html/TFile.html">TFile</a>.
 * 
 * This branch contains model parameters and model description.
 * @param [in] aTree pointer to the current <a href="http://root.cern.ch/root/html/TTree.html">TTree</a> object.
 * 
 * @fn void Model::SetParticlePX(Particle* aParticle)
 * @brief Fills the Particle with the coordinates and four momentum generated by the last call of GetIntegrand.
 *
 * @fn double Model::GetHyperCubeVolume()
 * @brief Returns the volume of the multidimensional cube of the Monte Carlo integration. 
 *
 * @fn const char* Model::GetHash()
 * @brief Returns the CRC32 hash generated form model parameters.
 * 
 * Hash is part of the filename that stores the integrand maximal values and average multiplicities 
 * for all particles types, e.g. <em>fmultiplicity_01234567.txt</em> if hash is 01234567.
 *
 * @fn const char* Model::GetName()
 * @brief Returns a name of the model.
 * 
 * @fn const char* Model::GetDescription()
 * @brief Returns a text description of the model and parameters for output files.
 *
 * @fn void Model::CreateEventSubDir()
 * @brief creates a model subdirectory structure in the main event directory (see sEventDIR)
 * 
 * @fn void Model::CalculateHash(TString aPreHash)
 * @brief calculates a CRC-32 hash from the given string of characters.
 * @param [in] aPreHash character string to hash
 * 
 * @var double Model::Xt
 * @brief particle creation time @f$ t @f$
 *
 * @var double Model::Xx
 * @brief particle @f$ x @f$ coordinate
 *
 * @var double Model::Xy
 * @brief particle @f$ y @f$ coordinate
 *
 * @var double Model::Xz
 * @brief particle @f$ z @f$ coordinate
 *
 * @var double Model::Pe
 * @brief particle energy @f$ E @f$
 *
 * @var double Model::Px
 * @brief particle @f$ p_x @f$ coordinate
 *
 * @var double Model::Py
 * @brief particle @f$ p_y @f$ coordinate
 *
 * @var double Model::Pz
 * @brief particle @f$ p_z @f$ coordinate
 *
 * @var double Model::mHyperCube
 * @brief volume of the Monte-Carlo hypercube
 *
 * @var char Model::mHash[9]
 * @brief 8-character hex-base CRC32 hash generated form model parameters
 *
 * @var char Model::mName[30]
 * @brief name of the model
 *
 * @var char Model::mDescription[2000]
 * @brief text description of model and parameters
 * 
 * @var TRandom2* Model::mRandom;
 * @brief random number generator, <a href="http://root.cern.ch/root/html/TRandom2.html">TRandom2</a> class
 */
/*! @def MODEL_PAR_DESC(_name_,_unit_)
 * @brief Macro used by models in creating the Model description. Writes the parameter <em>_name_</em> and its <em>_unit_</em> with 
 * a minimum number of characters.
 *
 * @def MODEL_NAME(_name_)
 * @brief Macro used by models in creating the Model description. Writes the model ID number and <em>_name_</em>.
 */
