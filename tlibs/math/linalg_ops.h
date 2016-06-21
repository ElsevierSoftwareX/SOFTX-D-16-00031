/*
 * Linalg operators
 * @author tweber
 * @date apr-2015
 * @license GPLv2 or GPLv3
 */

#ifndef __TL_LINALG_OPS_H__
#define __TL_LINALG_OPS_H__

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "../helper/traits.h"

namespace tl {

namespace ublas = boost::numeric::ublas;

template<class T1, class T2,
	LinalgType ty1=get_linalg_type<T1>::value,
	LinalgType ty2=get_linalg_type<T2>::value>
struct linalg_mult_op_impl
{
	void operator()(const T1&, const T2&) const
	{
		throw Err("No implementation for linalg_mult_op found.");
	}
};

// vec * vec
template<class T1, class T2>
struct linalg_mult_op_impl<T1, T2,
	LinalgType::VECTOR, LinalgType::VECTOR>
{
	typedef typename T1::value_type ret_type;

	ret_type operator()(const T1& vec1, const T2& vec2) const
	{
		return ublas::inner_prod(vec1, vec2);
	}
};

// mat * mat
template<class T1, class T2>
struct linalg_mult_op_impl<T1, T2,
	LinalgType::MATRIX, LinalgType::MATRIX>
{
	typedef T1 ret_type;

	ret_type operator()(const T1& mat1, const T2& mat2) const
	{
		return ublas::prod(mat1, mat2);
	}
};

// mat * vec
template<class T1, class T2>
struct linalg_mult_op_impl<T1, T2,
	LinalgType::MATRIX, LinalgType::VECTOR>
{
	typedef T2 ret_type;

	ret_type operator()(const T1& mat, const T2& vec) const
	{
		return ublas::prod(mat, vec);
	}
};

}

template<class T1, class T2>
typename tl::linalg_mult_op_impl<T1, T2>::ret_type operator*(const T1& t1, const T2& t2)
{
	return tl::linalg_mult_op_impl<T1, T2>()(t1, t2);
}

#endif
