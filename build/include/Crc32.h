#ifndef _TH2_CRC32_H_
  #define _TH2_CRC32_H_

class Crc32 {
  public:
    Crc32();
    Crc32(const char* aData, unsigned int aSize);
    ~Crc32();

    void	 Update(const char* aData, unsigned int aSize);
    unsigned int Finish() const;
    unsigned int GetValue() const;
    const std::string GetValueHex() const;

  private:
    mutable unsigned int      mCrc32;
    static const unsigned int mCrc32Tab[];
};

#endif

/*! @file Crc32.h
 * @brief Definition of Crc32 class. Generates a 32-bit CRC (Cyclic Redundancy Check) checksum.
 */
/*! @class Crc32
 * @brief Class for generating a 32-bit CRC (Cyclic Redundancy Check) checksum.
 * 
 * for more information see <a href=http://en.wikipedia.org/wiki/Cyclic_redundancy_check>CRC @ Wikipedia</a> <br />
 * CRC-32-IEEE 802.3 polynomial:
 * @f[ X^{32} + X^{26} + X^{23} + X^{22} + X^{16} + X^{12} + X^{11} + X^{10} + X^8 + X^7 + X^5 + X^4 + X^2 + X + 1 @f]
 *
 * @fn Crc32::Crc32()
 * @brief Default constructor.
 *
 * @fn Crc32::Crc32(const char* aData, unsigned int aSize)
 * @brief Initiate the calculation of the CRC checksum from a character string 
 * of a given size.
 * @param [in] aData character string
 * @param [in] aSize length of the string
 *
 * @fn Crc32::~Crc32()
 * @brief Destructor.
 *
 * @fn void Crc32::Update(const char* aData, unsigned int aSize)
 * @brief Updates the value of pre-CRC checksum with a new character string of a given lenght.
 * @param [in] aData character string
 * @param [in] aSize length of the string
 *
 * @fn unsigned int Crc32::Finish() const
 * @brief Finishes and returns the value of the CRC checksum as an unsigned 32-bit integer.
 * @warning After calling this of function further calling of Crc32::Update() will produce 
 * in invalid checksum value.
 *
 * @fn unsigned int Crc32::GetValue() const
 * @brief Returns the value of the CRC checksum as an unsigned 32-bit integer.
 *
 * @fn const char* Crc32::GetValueHex() const
 * @brief returns the value of the CRC checksum converted to as 8-character hex string. 
 */
