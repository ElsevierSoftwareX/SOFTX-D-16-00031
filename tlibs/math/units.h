/*
 * wrapper for boost.units
 * @author Tobias Weber
 * @date dec-2015
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIBS_UNITS__
#define __TLIBS_UNITS__


#include <boost/units/unit.hpp>
#include <boost/units/quantity.hpp>
#include <boost/units/dimensionless_quantity.hpp>
#include <boost/units/cmath.hpp>
#include <boost/units/physical_dimensions.hpp>
#include <boost/units/io.hpp>

#include <boost/units/systems/si.hpp>
#include <boost/units/systems/angle/degrees.hpp>
#include <boost/units/systems/si/codata/universal_constants.hpp>
#include <boost/units/systems/si/codata/neutron_constants.hpp>
#include <boost/units/systems/si/codata/electron_constants.hpp>
#include <boost/units/systems/si/codata/electromagnetic_constants.hpp>
#include <boost/units/systems/si/codata/physico-chemical_constants.hpp>

#include <boost/numeric/ublas/vector.hpp>


namespace tl {

namespace units = boost::units;
namespace co = boost::units::si::constants::codata;


// general quantities
template<class Sys, class T=double> using t_length =
	units::quantity<units::unit<units::length_dimension, Sys>, T>;
template<class Sys, class T=double> using t_momentum =
	units::quantity<units::unit<units::momentum_dimension, Sys>, T>;
template<class Sys, class T=double> using t_wavenumber =
	units::quantity<units::unit<units::wavenumber_dimension, Sys>, T>;
template<class Sys, class T=double> using t_velocity =
	units::quantity<units::unit<units::velocity_dimension, Sys>, T>;
template<class Sys, class T=double> using t_frequency =
	units::quantity<units::unit<units::frequency_dimension, Sys>, T>;
template<class Sys, class T=double> using t_energy =
	units::quantity<units::unit<units::energy_dimension, Sys>, T>;
template<class Sys, class T=double> using t_angle =
	units::quantity<units::unit<units::plane_angle_dimension, Sys>, T>;
template<class Sys, class T=double> using t_temperature =
	units::quantity<units::unit<units::temperature_dimension, Sys>, T>;
template<class Sys, class T=double> using t_mass =
	units::quantity<units::unit<units::mass_dimension, Sys>, T>;
template<class Sys, class T=double> using t_time =
	units::quantity<units::unit<units::time_dimension, Sys>, T>;
template<class Sys, class T=double> using t_flux =
	units::quantity<units::unit<units::magnetic_flux_density_dimension, Sys>, T>;
template<class Sys, class T=double> using t_area =
	units::quantity<units::unit<units::area_dimension, Sys>, T>;
template<class Sys, class T=double> using t_volume =
	units::quantity<units::unit<units::volume_dimension, Sys>, T>;

template<class Sys, class T=double> using t_length_inverse =
	units::quantity<units::unit<units::derived_dimension<units::length_base_dimension, -1>::type, Sys>, T>;
template<class Sys, class T=double> using t_length_square =
	units::quantity<units::unit<units::derived_dimension<units::length_base_dimension, 2>::type, Sys>, T>;
template<class Sys, class T=double> using t_momentum_square =
	units::quantity<units::unit<units::derived_dimension<units::momentum_dimension, 2>::type, Sys>, T>;
template<class Sys, class T=double> using t_action =
	units::quantity<units::unit<typename units::derived_dimension
	<units::mass_base_dimension,1, units::length_base_dimension,2, units::time_base_dimension,-1>::type, Sys>, T>;
template<class Sys, class T=double> using t_energy_per_temperature =
	units::quantity<units::unit<typename units::derived_dimension
	<units::mass_base_dimension,1, units::length_base_dimension,2,
	units::time_base_dimension,-2, units::temperature_base_dimension,-1>::type, Sys>, T>;
template<class Sys, class T=double> using t_energy_per_field =
	units::quantity<units::unit<typename units::derived_dimension
	<units::current_base_dimension,1, units::length_base_dimension,2>::type, Sys>, T>;
template<class Sys, class T=double> using t_dimensionless =
	units::quantity<units::unit<units::dimensionless_type, Sys>, T>;


// synonyms
template<class Sys, class T=double> using t_freq = t_frequency<Sys, T>;
template<class Sys, class T=double> using t_temp = t_temperature<Sys, T>;


// si quantities -- partial specialisations
template<class Y=double> using t_length_si = t_length<units::si::system, Y>;
template<class Y=double> using t_length_inverse_si = t_length_inverse<units::si::system, Y>;
template<class Y=double> using t_momentum_si = t_momentum<units::si::system, Y>;
template<class Y=double> using t_wavenumber_si = t_wavenumber<units::si::system, Y>;
template<class Y=double> using t_velocity_si = t_velocity<units::si::system, Y>;
template<class Y=double> using t_frequency_si = t_frequency<units::si::system, Y>;
template<class Y=double> using t_energy_si = t_energy<units::si::system, Y>;
template<class Y=double> using t_angle_si = t_angle<units::si::system, Y>;
template<class Y=double> using t_temperature_si = t_temperature<units::si::system, Y>;
template<class Y=double> using t_mass_si = t_mass<units::si::system, Y>;
template<class Y=double> using t_time_si = t_time<units::si::system, Y>;
template<class Y=double> using t_flux_si = t_flux<units::si::system, Y>;
template<class Y=double> using t_area_si = t_area<units::si::system, Y>;
template<class Y=double> using t_action_si = t_action<units::si::system, Y>;
template<class Y=double> using t_energy_per_temperature_si = t_energy_per_temperature<units::si::system, Y>;


// si quantities -- full specialisations
using length = t_length_si<>;
using inv_length = decltype(1./length());
using momentum = t_momentum_si<>;
using wavenumber = t_wavenumber_si<>;
using velocity = t_velocity_si<>;
using frequency = t_frequency_si<>;
using energy = t_energy_si<>;
using angle = t_angle_si<>;
using temperature = t_temperature_si<>;
using mass = t_mass_si<>;
using time = t_time_si<>;
using flux = t_flux_si<>;
using area = t_area_si<>;
using action = t_action_si<>;


// synonyms
typedef frequency freq;
typedef temperature temp;


// constants
template<class Y=double> t_energy<units::si::system, Y> get_one_meV()
{ return Y(1e-3) * Y(co::e/units::si::coulombs)*units::si::coulombs*units::si::volts; }
template<class Y=double> t_energy<units::si::system, Y> get_one_eV()
{ return Y(1) * Y(co::e/units::si::coulombs)*units::si::coulombs*units::si::volts; }
template<class Y=double> t_energy<units::si::system, Y> get_one_MeV()
{ return Y(1e6) * Y(co::e/units::si::coulombs)*units::si::coulombs*units::si::volts; }
template<class Y=double> t_length<units::si::system, Y> get_one_angstrom()
{ return Y(1e-10) * units::si::meters; }
template<class Y=double> t_length<units::si::system, Y> get_one_meter()
{ return Y(1) * units::si::meters; }
template<class Y=double> t_area<units::si::system, Y> get_one_barn()
{ return Y(1e-28) * units::si::meters*units::si::meters; }
template<class Y=double> t_temperature<units::si::system, Y> get_one_kelvin()
{ return Y(1) * units::si::kelvins; }
template<class Y=double> t_length<units::si::system, Y> get_one_centimeter()
{ return Y(1e-2) * units::si::meters; }
template<class Y=double> t_time<units::si::system, Y> get_one_second()
{ return Y(1) * units::si::seconds; }
template<class Y=double> t_time<units::si::system, Y> get_one_picosecond()
{ return Y(1e-12) * units::si::seconds; }
template<class Y=double> t_angle<units::si::system, Y> get_one_radian()
{ return Y(1) * units::si::radians; }
template<class Y=double> t_flux<units::si::system, Y> get_one_tesla()
{ return Y(1) * units::si::teslas; }
template<class Y=double> t_flux<units::si::system, Y> get_one_kilogauss()
{ return Y(0.1) * units::si::teslas; }

template<class Y=double> t_mass<units::si::system, Y> get_m_n()
{ return Y(co::m_n/units::si::kilograms)*units::si::kilograms; }
template<class Y=double> t_mass<units::si::system, Y> get_amu()
{ return Y(co::m_u/units::si::kilograms)*units::si::kilograms; }
template<class Y=double> t_action<units::si::system, Y> get_hbar()
{ return Y(co::hbar/units::si::joules/units::si::seconds)*units::si::joules*units::si::seconds; }
template<class Y=double> t_action<units::si::system, Y> get_h()
{ return get_hbar<Y>() * Y(2)*get_pi<Y>(); }
template<class Y=double> t_velocity<units::si::system, Y> get_c()
{ return Y(co::c/units::si::meters*units::si::seconds)*units::si::meters/units::si::seconds; }
template<class Y=double> t_energy_per_temperature<units::si::system, Y> get_kB()
{ return Y(co::k_B*units::si::kelvin/units::si::joules)/units::si::kelvin*units::si::joules; }
template<class Y=double> t_energy_per_field<units::si::system, Y> get_muB()
{ return Y(co::mu_B/units::si::joules*units::si::tesla)*units::si::joules/units::si::tesla; }
template<class Y=double> t_energy_per_field<units::si::system, Y> get_mu_n()
{ return Y(co::mu_n/units::si::joules*units::si::tesla)*units::si::joules/units::si::tesla; }
template<class Y=double> t_energy_per_field<units::si::system, Y> get_mu_N()
{ return Y(co::mu_N/units::si::joules*units::si::tesla)*units::si::joules/units::si::tesla; }
template<class Y=double> Y get_g_n() { return Y(co::g_n.value()); }
template<class Y=double> Y get_g_e() { return Y(co::g_e.value()); }


// template constants
#if __cplusplus >= 201402L
	template<class Y=double> const t_length_si<Y> t_meters = get_one_meter<Y>();
	template<class Y=double> const t_flux_si<Y> t_teslas = Y(1)*units::si::teslas;
	template<class Y=double> const t_time_si<Y> t_seconds = Y(1)*units::si::seconds;
	template<class Y=double> const t_angle_si<Y> t_radians = Y(1)*units::si::radians;
	template<class Y=double> const t_temperature_si<Y> t_kelvins = Y(1)*units::si::kelvins;
	template<class Y=double> const t_area_si<Y> t_barns = Y(1e-28)*units::si::meters*units::si::meters;

	template<class Y=double> const t_energy_si<Y> t_meV = get_one_meV<Y>();
	template<class Y=double> const t_length_si<Y> t_angstrom = get_one_angstrom<Y>();
#endif


// constants
static const length meters = 1.*units::si::meters;
static const flux teslas = 1.*units::si::teslas;
static const time seconds = 1.*units::si::seconds;
static const angle radians = 1.*units::si::radians;
static const temp kelvins = 1.*units::si::kelvins;
static const mass amu = co::m_u;
static const area barns = 1e-28 * units::si::meters*units::si::meters;

static const energy one_meV = 1e-3 * co::e * units::si::volts;
static const energy one_eV = co::e * units::si::volts;
static const length angstrom = 1e-10 * meters;
static const length cm = meters/100.;
static const time ps = 1e-12 * seconds;


// synonyms
static const temp kelvin = kelvins;
static const length meter = meters;
static const time second = seconds;
static const energy meV = one_meV;
static const flux tesla = teslas;
static const area barn = barns;



// helper functions
template<class t_quant>
t_quant my_units_sqrt(const decltype(t_quant() * t_quant())& val)
{
	using t_quant_sq = decltype(t_quant() * t_quant());
	using Y = typename t_quant::value_type;

	t_quant one_quant = t_quant::from_value(Y(1));
	t_quant_sq one_quant_sq = t_quant_sq::from_value(Y(1));

	Y valsq = Y(val / one_quant_sq);
	return std::sqrt(valsq) * one_quant;
}

template<class t_quant>
decltype(t_quant()*t_quant()) my_units_pow2(const t_quant& val)
{
	return val*val;
}

template<class t_elem, template<class...> class t_vec=boost::numeric::ublas::vector>
t_elem my_units_norm2(const t_vec<t_elem>& vec)
{
	using t_elem_sq = decltype(t_elem()*t_elem());
	t_elem_sq tRet = t_elem_sq();

	for(std::size_t i=0; i<vec.size(); ++i)
		tRet += vec[i]*vec[i];

	return tl::my_units_sqrt<t_elem>(tRet);
}

}
#endif
