/*
 * Custom type traits
 * @author Tobias Weber
 * @date 19-nov-2014
 * @license GPLv2 or GPLv3
 */

#ifndef __MY_TRAITS_H__
#define __MY_TRAITS_H__

#include <type_traits>


namespace tl {

// -----------------------------------------------------------------------------
template<class T, bool bScalar=std::is_scalar<T>::value>
struct underlying_value_type
{};

template<class T>
struct underlying_value_type<T, 1>
{
	using value_type = T;
};

template<class T>
struct underlying_value_type<T, 0>
{
	using value_type = typename T::value_type;
};

template<class T>
using underlying_value_type_t =
	typename underlying_value_type<T, std::is_scalar<T>::value>::value_type;

// -----------------------------------------------------------------------------



// -----------------------------------------------------------------------------
}
#include <vector>
#include <array>
#include <list>
#include <initializer_list>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/math/quaternion.hpp>

namespace tl {

typedef std::integral_constant<int, 0> dim_0d_type;
typedef std::integral_constant<int, 1> dim_1d_type;
typedef std::integral_constant<int, 2> dim_2d_type;

template<class> struct get_type_dim : dim_0d_type {};

template<class... PARAMS> struct get_type_dim<std::vector<PARAMS...>> : dim_1d_type {};
template<class... PARAMS> struct get_type_dim<std::array<PARAMS...>> : dim_1d_type {};
template<class... PARAMS> struct get_type_dim<std::list<PARAMS...>> : dim_1d_type {};
template<class... PARAMS> struct get_type_dim<boost::numeric::ublas::vector<PARAMS...>> : dim_1d_type {};
template<class... PARAMS> struct get_type_dim<std::initializer_list<PARAMS...>> : dim_1d_type {};

template<class... PARAMS> struct get_type_dim<boost::numeric::ublas::matrix<PARAMS...>> : dim_2d_type {};



enum class LinalgType : short
{
	UNKNOWN,
	VECTOR,
	MATRIX,
	QUATERNION,
	REAL
};

template<LinalgType val>  struct linalg_type { static constexpr LinalgType value = val; };

template<class> struct get_linalg_type : linalg_type<LinalgType::UNKNOWN> {};
template<class... PARAMS> struct get_linalg_type<std::vector<PARAMS...>> : linalg_type<LinalgType::VECTOR> {};
template<class... PARAMS> struct get_linalg_type<boost::numeric::ublas::vector<PARAMS...>> : linalg_type<LinalgType::VECTOR> {};
template<class... PARAMS> struct get_linalg_type<boost::numeric::ublas::matrix<PARAMS...>> : linalg_type<LinalgType::MATRIX> {};
template<class... PARAMS> struct get_linalg_type<boost::math::quaternion<PARAMS...>> : linalg_type<LinalgType::QUATERNION> {};
template<> struct get_linalg_type<double> : linalg_type<LinalgType::REAL> {};
template<> struct get_linalg_type<float> : linalg_type<LinalgType::REAL> {};
// -----------------------------------------------------------------------------



// -----------------------------------------------------------------------------
}

#include <boost/units/dimensionless_quantity.hpp>

namespace tl {
enum class ScalarType : short
{
	TRIVIAL,
	DIMENSIONLESS,
	QUANTITY
};

template<ScalarType val> struct scalar_type { static constexpr ScalarType value = val; };

template<class T> struct get_scalar_type : scalar_type<ScalarType::TRIVIAL> { using value_type = T; };
template<class Sys, class T> struct get_scalar_type<boost::units::dimensionless_quantity<Sys, T>> : scalar_type<ScalarType::DIMENSIONLESS> { using value_type = T; };
template<class Sys, class T> struct get_scalar_type<boost::units::quantity<Sys, T>> : scalar_type<ScalarType::QUANTITY> 
{ using value_type = typename boost::units::quantity<Sys, T>::value_type; };
// -----------------------------------------------------------------------------




// -----------------------------------------------------------------------------
template<class T>
struct remove_constref
{
	typedef typename std::remove_const<
		typename std::remove_reference<T>::type
			>::type type;
};

// like C++14 style
template<class T>
using remove_constref_t = typename remove_constref<T>::type;
// -----------------------------------------------------------------------------


}

#endif
