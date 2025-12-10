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

#ifndef _TH2_TEMPERATURE_H_
  #define _TH2_TEMPERATURE_H_

#include "Vector3D.h"

class Temperature
{
  public:
    Temperature();
    ~Temperature();

    int    GetTemperatureType() const;
    double GetTemperature() const;
    double GetTemperature(double aX, double aY, double aZ) const;
    void   SetTemperature(double aTemp);
    void   SetTemperature(Vector3D* aTemp);

  private:
    double    mTempConst;
    Vector3D* mTempVar;
    int       mTempType;
};

#endif

/*! @file Temperature.h
 * @brief Definition of Temperature class. Gives the local temperature.
 */
/*! @class Temperature
 * @brief Class gives the local temperature.
 *
 * If SetTemperature() is called with a Vector3D object it will return the temperature at @f$ T(x,y,z) @f$ 
 * in other case it will return a constant value of temperature.
 *
 * @fn Temperature::Temperature()
 * @brief Default constructor.
 *
 * @fn Temperature::~Temperature()
 * @brief Destructor.
 *
 * @fn int Temperature::GetTemperatureType() const
 * @brief Get the type of temperature stored.
 * @retval 0 constant, @f$ T = const @f$
 * @retval 1 variable, @f$ T = T(x,y,z) @f$
 * @retval -1 unknown, @f$ T = 0 @f$
 *
 * @fn double Temperature::GetTemperature() const
 * @brief Returns constant temperature value
 *
 * @fn double Temperature::GetTemperature(double aX, double aY, double aZ) const
 * @brief Returns temperature as a function of <b>x</b>,<b>y</b> and <b>z</b>
 * @param [in] aX x-coordinate
 * @param [in] aY y-coordinate
 * @param [in] aZ z-coordinate
 *
 * @fn void Temperature::SetTemperature(double aTemp)
 * @brief Sets the temperature type to constant type.
 * @param [in] aTemp constant temperature value
 *
 * @fn void Temperature::SetTemperature(Vector3D* aTemp)
 * @brief Sets the temperature type to variable.
 * @param [in] aTemp pointer to Vector3D object with the temperature data
 */
