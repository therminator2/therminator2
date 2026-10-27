#ifndef COULOMB_UNITS
#define COULOMB_UNITS

// Physical quantities of the Coulomb afterburner with units, using the mp-units library (external/mp-units).
// A quantity carries its unit in its type: formulas that are not dimensionally consistent (adding a length to a time,
// a missing factor c, ...) do not compile, and conversions between units are done automatically. Numbers without
// units appear only at the boundary to the rest of THERMINATOR (Particle, histograms, printouts), through the
// functions below, which state the units explicitly.

#include <cmath>
#include <ostream>
#include <type_traits>

#include <mp-units/math.h>
#include <mp-units/systems/si.h>

#include "Particle.h"

namespace CoulombUnits {

using namespace mp_units;

inline constexpr auto fm = si::femto<si::metre>;
inline constexpr auto c = si::si2019::speed_of_light_in_vacuum;
inline constexpr auto fm_c = fm / c; // unit of time, fm/c
inline constexpr auto MeV = si::mega<si::electronvolt>;
inline constexpr auto GeV = si::giga<si::electronvolt>;

// THERMINATOR's internal unit of positions and times, GeV^-1 = hbar c / GeV (/ c), with hbar c = kHbarC GeV fm
// (THGlobal.h). The ratio must be the same as kHbarC, which CheckInternalUnits() verifies.
inline constexpr struct inverse_gev_length final : named_unit<"hbar_c/GeV", mag_ratio<1973269631, 10000000000> * fm> {
} inverse_gev_length;
inline constexpr auto inverse_gev_time = inverse_gev_length / c;

// Quantities used by the afterburner, with the units in which they are stored
using Dimensionless = quantity<one, double>;
using Length = quantity<fm, double>;
using Area = quantity<square(fm), double>;           // e.g. Minkowski intervals c^2 dt^2 - dx^2
using Time = quantity<fm_c, double>;
using Speed = quantity<c, double>;                   // velocities in units of c
using Acceleration = quantity<c / fm_c, double>;
using Energy = quantity<GeV, double>;
using Mass = quantity<GeV / square(c), double>;
using Momentum = quantity<GeV / c, double>;
using Force = quantity<GeV / fm, double>;
using Power = quantity<GeV / fm_c, double>;
using Coupling = quantity<GeV * fm, double>;         // energy times length, e.g. k e^2 = alpha hbar c

// Coulomb constant with unit charges, k e^2 = e^2 / (4 pi epsilon_0) = alpha hbar c = 1.44 MeV fm,
// from epsilon_0 = 55.26349406 e^2 / (GeV fm)
inline constexpr Coupling kCoulomb = (1. / (4. * 3.14159265358979323846 * 55.26349406)) * (GeV * fm);

// ---------------------------------------------------------------------------------------------------------------------
// Three-vectors of quantities. Products and magnitudes have the units derived from those of the components.

template<typename Q>
struct Vec3 {
    Q x, y, z;

    Vec3() : x(Q::zero()), y(Q::zero()), z(Q::zero()) {}
    Vec3(Q aX, Q aY, Q aZ) : x(aX), y(aY), z(aZ) {}
    // From vectors in other units of the same quantity
    template<typename R>
        requires std::is_convertible_v<R, Q>
    Vec3(const Vec3<R> &aOther) : x(aOther.x), y(aOther.y), z(aOther.z) {}

    template<typename R>
    Vec3 &operator+=(const Vec3<R> &aOther) { x += aOther.x; y += aOther.y; z += aOther.z; return *this; }
    template<typename R>
    Vec3 &operator-=(const Vec3<R> &aOther) { x -= aOther.x; y -= aOther.y; z -= aOther.z; return *this; }
};

template<typename T> struct IsVec3 : std::false_type {};
template<typename Q> struct IsVec3<Vec3<Q>> : std::true_type {};

template<typename A, typename B>
auto operator+(const Vec3<A> &a, const Vec3<B> &b) -> Vec3<decltype(a.x + b.x)> { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
template<typename A, typename B>
auto operator-(const Vec3<A> &a, const Vec3<B> &b) -> Vec3<decltype(a.x - b.x)> { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
template<typename A>
Vec3<A> operator-(const Vec3<A> &a) { return {-a.x, -a.y, -a.z}; }
// Multiplication and division by a scalar: a number, a quantity or a unit
template<typename A, typename S>
    requires(!IsVec3<S>::value)
auto operator*(const Vec3<A> &a, const S &s) -> Vec3<decltype(a.x * s)> { return {a.x * s, a.y * s, a.z * s}; }
template<typename S, typename A>
    requires(!IsVec3<S>::value)
auto operator*(const S &s, const Vec3<A> &a) -> Vec3<decltype(s * a.x)> { return {s * a.x, s * a.y, s * a.z}; }
template<typename A, typename S>
    requires(!IsVec3<S>::value)
auto operator/(const Vec3<A> &a, const S &s) -> Vec3<decltype(a.x / s)> { return {a.x / s, a.y / s, a.z / s}; }

template<typename A, typename B>
auto Dot(const Vec3<A> &a, const Vec3<B> &b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
template<typename A, typename B>
auto Cross(const Vec3<A> &a, const Vec3<B> &b) -> Vec3<decltype(a.y * b.z - a.z * b.y)>
{
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
template<typename A>
auto Mag2(const Vec3<A> &a) { return Dot(a, a); }
template<typename A>
A Mag(const Vec3<A> &a) { return A(sqrt(Mag2(a))); }
// Unit vector; the zero vector for a zero vector (like TVector3::Unit)
template<typename A>
Vec3<Dimensionless> UnitVector(const Vec3<A> &a)
{
    A tMag = Mag(a);
    if (tMag == A::zero()) return Vec3<Dimensionless>();
    return a / tMag;
}

template<typename Q>
bool IsFinite(const Q &q) { return std::isfinite(q.numerical_value_in(q.unit)); }
template<typename Q>
bool IsFinite(const Vec3<Q> &v) { return IsFinite(v.x) && IsFinite(v.y) && IsFinite(v.z); }
template<typename Q>
Q NaN() { return std::numeric_limits<double>::quiet_NaN() * Q::unit; }

template<typename Q>
std::ostream &operator<<(std::ostream &os, const Vec3<Q> &v) { return os << "(" << v.x << ", " << v.y << ", " << v.z << ")"; }

// ---------------------------------------------------------------------------------------------------------------------
// Boundary to Particle, which stores positions and times in GeV^-1, energies, momenta and masses in GeV

inline Time TimeOf(const Particle &aParticle) { return aParticle.T() * inverse_gev_time; }
inline Vec3<Length> PositionOf(const Particle &aParticle)
{
    return {aParticle.X() * inverse_gev_length, aParticle.Y() * inverse_gev_length, aParticle.Z() * inverse_gev_length};
}
inline Energy EnergyOf(const Particle &aParticle) { return aParticle.E() * GeV; }
inline Vec3<Momentum> MomentumOf(const Particle &aParticle)
{
    return {aParticle.Px() * GeV / c, aParticle.Py() * GeV / c, aParticle.Pz() * GeV / c};
}
inline Mass MassOf(const Particle &aParticle) { return aParticle.mass * GeV / square(c); }
// Velocity p c^2 / E; zero for a particle with zero energy
inline Vec3<Speed> VelocityOf(const Vec3<Momentum> &aMomentum, Energy aEnergy)
{
    if (aEnergy == Energy::zero()) return Vec3<Speed>();
    return aMomentum * square(c) / aEnergy;
}
inline Vec3<Speed> VelocityOf(const Particle &aParticle) { return VelocityOf(MomentumOf(aParticle), EnergyOf(aParticle)); }
// Energy sqrt(m^2 c^4 + p^2 c^2)
inline Energy EnergyOf(Mass aMass, const Vec3<Momentum> &aMomentum)
{
    return Energy(sqrt(pow<2>(aMass * square(c)) + Mag2(aMomentum) * square(c)));
}

inline void SetSpaceTime(Particle &aParticle, Time aTime, const Vec3<Length> &aPosition)
{
    aParticle.SetTXYZ(aTime.numerical_value_in(inverse_gev_time), aPosition.x.numerical_value_in(inverse_gev_length),
                      aPosition.y.numerical_value_in(inverse_gev_length), aPosition.z.numerical_value_in(inverse_gev_length));
}
inline void SetFourMomentum(Particle &aParticle, Energy aEnergy, const Vec3<Momentum> &aMomentum)
{
    aParticle.SetEPxPyPz(aEnergy.numerical_value_in(GeV), aMomentum.x.numerical_value_in(GeV / c),
                         aMomentum.y.numerical_value_in(GeV / c), aMomentum.z.numerical_value_in(GeV / c));
}

// Positions and times of THERMINATOR are in GeV^-1 with hbar c = kHbarC: returns false if inverse_gev_length differs
bool CheckInternalUnits();

} // namespace CoulombUnits

#endif
