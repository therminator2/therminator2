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

#include "Entropy.h"

Entropy::Entropy()
: mEntropyConst(0.0), mEntropyVar(0), mEntropyType(-1)
{
}

Entropy::~Entropy()
{
  if(mEntropyVar)
    delete mEntropyVar;
}

double Entropy::GetEntropy() const
{
  return mEntropyConst;
}

double Entropy::GetEntropy(double aX, double aY, double aZ) const
{
  return mEntropyVar->Interpolate(aX, aY, aZ);
}

int Entropy::GetEntropyType() const
{
  return mEntropyType;
}

void Entropy::SetEntropy(double aEntropy)
{
  mEntropyType = 0;
  mEntropyConst = aEntropy;
}

void Entropy::SetEntropy(Vector3D* aEntropy)
{
  mEntropyType = 1;
  mEntropyVar = aEntropy;
}
