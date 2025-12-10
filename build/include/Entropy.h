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

#ifndef _TH2_ENTROPY_H_
  #define _TH2_ENTROPY_H_

#include "Vector3D.h"

class Entropy
{
  public:
    Entropy();
    ~Entropy();

    int    GetEntropyType() const;
    double GetEntropy() const;
    double GetEntropy(double aX, double aY, double aZ) const;
    void   SetEntropy(double aEntropy);
    void   SetEntropy(Vector3D* aEntropy);

  private:
    double    mEntropyConst;
    Vector3D* mEntropyVar;
    int       mEntropyType;
};

#endif

/*! @file Entropy.h
 * @brief Definition of Entropy class. Gives the local entropy density @f$ s @f$.
 */
/*! @class Entropy
 * @brief Class gives the local entropy density @f$ s @f$.
 *
 * If SetEntropy() is called with a Vector3D object it will return entropy density at @f$ s(x,y,z) @f$ 
 * in other case it will return a constant value of entropy density.
 *
 * @fn Entropy::Entropy()
 * @brief Default constructor.
 *
 * @fn Entropy::~Entropy()
 * @brief Destructor.
 *
 * @fn int Entropy::GetEntropyType() const
 * @brief Get the type of entropy density stored.
 * @retval 0 constant, @f$ s = const @f$
 * @retval 1 variable, @f$ s = s(x,y,z) @f$
 * @retval -1 unknown, @f$ s = 0 @f$
 *
 * @fn double Entropy::GetEntropy() const
 * @brief Returns constant entropy density value
 *
 * @fn double Entropy::GetEntropy(double aX, double aY, double aZ) const
 * @brief Returns entropy density as a function of <b>x</b>,<b>y</b> and <b>z</b>
 * @param [in] aX x-coordinate
 * @param [in] aY y-coordinate
 * @param [in] aZ z-coordinate
 *
 * @fn void Entropy::SetEntropy(double aEntropy)
 * @brief Sets the entropy density type to constant.
 * @param [in] aEntropy constant entropy density value
 *
 * @fn void Entropy::SetEntropy(Vector3D* aEntropy)
 * @brief Sets the entropy density type to variable.
 * @param [in] aEntropy pointer to Vector3D object with the entropy density data
 */
