/*
 * basic linalg helpers
 *
 * @author: tweber
 * @date: 30-apr-2013
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIBS_LINALG_H__
#define __TLIBS_LINALG_H__

#include "../helper/flags.h"
#include "../helper/exception.h"
#include "math.h"
#include "../log/log.h"
#include "../log/debug.h"
#include "../helper/traits.h"

#include <initializer_list>
#include <cmath>

#include <boost/algorithm/minmax_element.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/exception.hpp>


namespace tl {

namespace ublas = boost::numeric::ublas;


template<class matrix_type=ublas::matrix<double> >
typename matrix_type::value_type determinant(const matrix_type& mat);


// creates a vector
template<class t_vec=ublas::vector<double>, template<class...> class t_lst=std::initializer_list>
t_vec make_vec(t_lst<typename t_vec::value_type>&& lst)
{
	using T = typename t_vec::value_type;
	using t_iter = typename t_lst<T>::const_iterator;

	t_vec vec(lst.size());

	std::size_t i=0;
	for(t_iter iter = lst.begin(); iter!=lst.end(); ++i, ++iter)
		vec[i] = *iter;

	return vec;
}

// creates a matrix
template<class t_mat=ublas::matrix<double>, template<class...> class t_lst=std::initializer_list>
t_mat make_mat(t_lst<t_lst<typename t_mat::value_type>>&& lst)
{
	using T = typename t_mat::value_type;

	std::size_t I = lst.size();
	std::size_t J = lst.begin()->size();

	t_mat mat(I, J);
	typename t_lst<t_lst<T>>::const_iterator iter = lst.begin();

	for(std::size_t i=0; i<I; ++i, ++iter)
	{
		typename t_lst<T>::const_iterator iterinner = iter->begin();
		for(std::size_t j=0; j<J; ++j, ++iterinner)
		{
			mat(i,j) = *iterinner;
		}
	}

	return mat;
}

template<class matrix_type = ublas::matrix<double>>
matrix_type unit_matrix(std::size_t N)
{
	matrix_type mat(N,N);

	for(std::size_t i=0; i<N; ++i)
		for(std::size_t j=0; j<N; ++j)
			mat(i,j) = (i==j ? 1 : 0);
	return mat;
}


// converts vector
template<class t_from, class t_to, template<class...> class t_vec = ublas::vector>
t_vec<t_to> convert_vec(const t_vec<t_from>& vec)
{
	using t_vec_to = t_vec<t_to>;

	t_vec_to vecRet(vec.size());

	for(std::size_t i=0; i<vec.size(); ++i)
		vecRet[i] = t_to(vec[i]);

	return vecRet;
}


template<class vec_type>
bool vec_equal(const vec_type& vec0, const vec_type& vec1,
	typename vec_type::value_type eps = std::numeric_limits<typename vec_type::value_type>::epsilon())
{
	typedef typename vec_type::value_type T;

	if(vec0.size() != vec1.size())
		return false;

	for(std::size_t i=0; i<vec0.size(); ++i)
		if(!float_equal<T>(vec0[i], vec1[i], eps))
			return false;
	return true;
}

template<class mat_type>
bool mat_equal(const mat_type& mat0, const mat_type& mat1,
	typename mat_type::value_type eps = std::numeric_limits<typename mat_type::value_type>::epsilon())
{
	typedef typename mat_type::value_type T;

	if(mat0.size1() != mat1.size1() || mat0.size2() != mat1.size2())
		return false;

	for(std::size_t i=0; i<mat0.size1(); ++i)
		for(std::size_t j=0; j<mat0.size2(); ++j)
			if(!float_equal<T>(mat0(i,j), mat1(i,j), eps))
				return false;
	return true;
}


template<class vec_type>
typename vec_type::value_type vec_len(const vec_type& vec)
{
	typename vec_type::value_type t = typename vec_type::value_type();

	for(std::size_t i=0; i<vec.size(); ++i)
		t += vec[i]*vec[i];

	t = std::sqrt(t);
	return t;
}



// remove an element from a vector
template<class vector_type>
vector_type remove_elem(const vector_type& vec, std::size_t iIdx)
{
	vector_type vecret(vec.size()-1);

	for(std::size_t i=0, j=0; i<vec.size() && j<vecret.size();)
	{
		vecret[j] = vec[i];

		if(i!=iIdx) ++j;
		++i;
	}

	return vecret;
}

template<class matrix_type>
matrix_type submatrix(const matrix_type& mat, std::size_t iRow, std::size_t iCol)
{
	matrix_type matret(mat.size1()-1, mat.size2()-1);

	for(std::size_t i=0, i0=0; i<mat.size1() && i0<matret.size1();)
	{
		for(std::size_t j=0, j0=0; j<mat.size2() && j0<matret.size2();)
		{
			matret(i0,j0) = mat(i,j);

			if(j!=iCol) ++j0;
			++j;
		}

		if(i!=iRow) ++i0;
		++i;
	}

	return matret;
}

template<class matrix_type>
matrix_type remove_column(const matrix_type& mat, std::size_t iCol)
{
	matrix_type matret(mat.size1(), mat.size2()-1);
	for(std::size_t i=0; i<mat.size1(); ++i)
	{
		for(std::size_t j=0, j0=0; j<mat.size2() && j0<matret.size2(); ++j)
		{
			matret(i,j0) = mat(i,j);
                        if(j!=iCol) ++j0;
		}
	}
	return matret;
}

template<class matrix_type>
void submatrix_copy(matrix_type& mat, const matrix_type& sub,
	std::size_t iRowBegin, std::size_t iColBegin)
{
	for(std::size_t i=0; i<sub.size1(); ++i)
		for(std::size_t j=0; j<sub.size2(); ++j)
			mat(iRowBegin+i, iColBegin+j) = sub(i,j);
}

template<class vec_type>
void subvector_copy(vec_type& vec, const vec_type& sub, std::size_t iRowBegin)
{
	for(std::size_t i=0; i<sub.size(); ++i)
		vec[iRowBegin+i] = sub[i];
}


template<class matrix_type>
matrix_type remove_elems(const matrix_type& mat, std::size_t iIdx)
{
	return submatrix(mat, iIdx, iIdx);
}

template<class t_vec=ublas::vector<double>,
	class t_mat=ublas::matrix<typename t_vec::value_type>>
void set_column(t_mat& M, std::size_t iCol, const t_vec& vec)
{
	std::size_t s = std::min(vec.size(), M.size1());
	for(std::size_t i=0; i<s; ++i)
		M(i, iCol) = vec[i];
}

template<class vector_type=ublas::vector<double>,
	class matrix_type=ublas::matrix<typename vector_type::value_type>>
vector_type get_column(const matrix_type& mat, std::size_t iCol)
{
	vector_type vecret(mat.size1());

	for(std::size_t i=0; i<mat.size1(); ++i)
		vecret[i] = mat(i, iCol);

	return vecret;
}

template<class vector_type = ublas::vector<double>,
	class matrix_type = ublas::matrix<typename vector_type::value_type>,
	class cont_type = std::vector<vector_type>>
cont_type get_columns(const matrix_type& mat)
{
	cont_type vec;
	vec.reserve(mat.size2());

	for(std::size_t i=0; i<mat.size2(); ++i)
		vec.push_back(get_column(mat, i));

	return vec;
}

template<class vector_type=ublas::vector<double>,
	class matrix_type=ublas::matrix<typename vector_type::value_type>>
vector_type get_row(const matrix_type& mat, std::size_t iRow)
{
	vector_type vecret(mat.size2());

	for(std::size_t i=0; i<mat.size2(); ++i)
		vecret[i] = mat(iRow, i);

	return vecret;
}


template<class t_mat=ublas::matrix<double>>
t_mat mirror_matrix(std::size_t iSize, std::size_t iComp)
{
	using T = typename t_mat::value_type;

	t_mat mat = unit_matrix<t_mat>(iSize);
	mat(iComp, iComp) = T(-1);

	return mat;
}


template<class matrix_type=ublas::matrix<double>>
matrix_type rotation_matrix_2d(typename matrix_type::value_type angle)
{
	typedef typename matrix_type::value_type T;

	T s, c;

	if(angle==0.)
	{
		s = T(0);
		c = T(1);
	}
	else
	{
		s = std::sin(angle);
		c = std::cos(angle);
	}

	return make_mat<matrix_type>
	({	{c, -s},
		{s,  c}});
}


/**
 * generates points in an arc defined by vec1 and vec2 at an angle phi around vec1
 */
template<class t_mat=ublas::matrix<double>, class t_vec=ublas::vector<double>>
t_vec arc(const t_vec& vec1, const t_vec& vec2, tl::underlying_value_type_t<t_vec> phi)
{
	return std::cos(phi)*vec1 + std::sin(phi)*vec2;
}


template<class matrix_type=ublas::matrix<double>>
matrix_type rotation_matrix_3d_x(typename matrix_type::value_type angle)
{
	typedef typename matrix_type::value_type T;

	T s, c;
	if(angle==0.)
	{
		s = T(0);
		c = T(1);
	}
	else
	{
		s = std::sin(angle);
		c = std::cos(angle);
	}

	return make_mat<matrix_type>
	({	{1, 0,  0},
		{0, c, -s},
		{0, s,  c}});
}

template<class matrix_type=ublas::matrix<double>>
matrix_type rotation_matrix_3d_y(typename matrix_type::value_type angle)
{
	typedef typename matrix_type::value_type T;

	T s, c;
	if(angle==0.)
	{
		s = T(0);
		c = T(1);
	}
	else
	{
		s = std::sin(angle);
		c = std::cos(angle);
	}

	return make_mat<matrix_type>
	({	{c,  0, s},
		{0,  1, 0},
		{-s, 0, c}});
}

template<class matrix_type=ublas::matrix<double>>
matrix_type rotation_matrix_3d_z(typename matrix_type::value_type angle)
{
	typedef typename matrix_type::value_type T;

	T s, c;
	if(angle==0.)
	{
		s = T(0);
		c = T(1);
	}
	else
	{
		s = std::sin(angle);
		c = std::cos(angle);
	}

	return make_mat<matrix_type>
	({	{c, -s, 0},
		{s,  c, 0},
		{0,  0, 1}});
}

template<class matrix_type = ublas::matrix<double>,
	class vector_type = ublas::vector<typename matrix_type::value_type>>
matrix_type skew(const vector_type& vec)
{
	if(vec.size() == 3)
	{
		return make_mat<matrix_type>
		({	{       0, -vec[2],  vec[1]},
			{  vec[2],       0, -vec[0]},
			{ -vec[1],  vec[0],       0}});
	}
	else
		throw Err("Skew only defined for three dimensions.");
}


template<class matrix_type = ublas::matrix<double>,
		class cont_type = std::initializer_list<typename matrix_type::value_type>>
matrix_type diag_matrix(const cont_type& lst)
{
	matrix_type mat(lst.size(), lst.size());

	for(std::size_t i=0; i<lst.size(); ++i)
		for(std::size_t j=0; j<lst.size(); ++j)
			mat(i,j) = 0.;

	std::size_t i = 0;
	for(typename cont_type::const_iterator iter=lst.begin(); iter!=lst.end(); ++iter, ++i)
		mat(i,i) = *iter;

	return mat;
}


// Euler-Rodrigues formula
template<class mat_type=ublas::matrix<double>,
	class vec_type=ublas::vector<typename mat_type::value_type>,
	typename T = typename mat_type::value_type>
mat_type rotation_matrix(const vec_type& _vec, T angle)
{
	const vec_type vec = _vec/ublas::norm_2(_vec);

	T s, c;
	if(angle==0.)
	{
		s = T(0);
		c = T(1);
	}
	else
	{
		s = std::sin(angle);
		c = std::cos(angle);
	}

	return (T(1) - c) * ublas::outer_prod(vec,vec) +
		c * unit_matrix(vec.size()) +
		s * skew(vec);
}

template<class matrix_type=ublas::matrix<double>>
typename matrix_type::value_type trace(const matrix_type& mat)
{
	typedef typename matrix_type::value_type T;

	if(mat.size1() != mat.size2())
		return T(0);

	T tr = T(0.);
	for(std::size_t i=0; i<mat.size1(); ++i)
		tr += mat(i,i);
	return tr;
}

// see: https://www.opengl.org/sdk/docs/man2/xhtml/gluPerspective.xml
// also see: similar gnomonic projection of spherical coordinates onto a plane
template<class matrix_type=ublas::matrix<double,ublas::row_major, ublas::bounded_array<double,4*4>>,
	class T=typename matrix_type::value_type>
matrix_type perspective_matrix(T yfov, T asp, T n, T f)
{
	const T y = std::tan(T(0.5)*get_pi<T>() - T(0.5)*yfov);
	const T x = y/asp;
	const T dsgn = -1.;

	return make_mat<matrix_type>
	({
		{    x, T(0),             T(0),              T(0) },
		{ T(0),    y,             T(0),              T(0) },
		{ T(0), T(0), dsgn*(f+n)/(f-n), (-T(2)*f*n)/(f-n) },
		{ T(0), T(0),             dsgn,              T(0) }
	});
}

// see: https://www.opengl.org/sdk/docs/man2/xhtml/glOrtho.xml
template<class matrix_type=ublas::matrix<double, ublas::row_major, ublas::bounded_array<double,4*4>>,
	class T=typename matrix_type::value_type>
matrix_type ortho_matrix(T l, T r, T b, T t, T n, T f)
{
	return make_mat<matrix_type>
	({
		{ T(2)/(r-l),       T(0),       T(0), (l+r)/(l-r) },
		{       T(0), T(2)/(t-b),       T(0), (b+t)/(b-t) },
		{       T(0),       T(0), T(2)/(n-f), (n+f)/(n-f) },
		{       T(0),       T(0),       T(0), T(1)        }
	});
}

// -----------------------------------------------------------------------------
template<typename T, class FKT, const int iDim=get_type_dim<T>::value>
struct is_nan_or_inf_impl
{
	is_nan_or_inf_impl(const FKT&) {}
	bool operator()(T) const { throw Err("No implementation of is_nan_or_inf!"); }
};

template<typename real_type, class FKT>
struct is_nan_or_inf_impl<real_type, FKT, 0>	// scalar impl.
{
	const FKT& m_fkt;
	is_nan_or_inf_impl(const FKT& fkt) : m_fkt(fkt) {}
	bool operator()(real_type d) const { return m_fkt(d); }
};

template<typename vec_type, class FKT>
struct is_nan_or_inf_impl<vec_type, FKT, 1>		// vector impl.
{
	const FKT& m_fkt;
	is_nan_or_inf_impl(const FKT& fkt) : m_fkt(fkt) {}

	bool operator()(const vec_type& vec) const
	{
		for(std::size_t i=0; i<vec.size(); ++i)
			if(m_fkt(vec[i]))
				return true;
		return false;
	}
};

template<typename mat_type, class FKT>
struct is_nan_or_inf_impl<mat_type, FKT, 2>		// matrix impl.
{
	const FKT& m_fkt;
	is_nan_or_inf_impl(const FKT& fkt) : m_fkt(fkt) {}

	bool operator()(const mat_type& mat) const
	{
		for(std::size_t i=0; i<mat.size1(); ++i)
			for(std::size_t j=0; j<mat.size2(); ++j)
				if(m_fkt(mat(i,j)))
					return true;
		return false;
	}
};

template<class T=ublas::matrix<double>>
bool isnan(const T& mat)
{
	typedef underlying_value_type_t<T> real_type;

	using fkt = std::function<bool(real_type)>;
	fkt stdisnan = (bool(*)(real_type))std::isnan;
	is_nan_or_inf_impl<T, fkt> _isnan(stdisnan);
	return _isnan(mat);
}

template<class T=ublas::matrix<double>>
bool isinf(const T& mat)
{
	typedef underlying_value_type_t<T> real_type;

	using fkt = std::function<bool(real_type)>;
	fkt stdisinf = (bool(*)(real_type))std::isinf;
	is_nan_or_inf_impl<T, fkt> _isinf(stdisinf);
	return _isinf(mat);
}

template<class T=ublas::matrix<double>>
bool is_nan_or_inf(const T& mat)
{
	typedef underlying_value_type_t<T> real_type;

	using fkt = std::function<bool(real_type)>;
	fkt stdisnaninf = [](real_type d)->bool { return std::isnan(d) || std::isinf(d); };
	is_nan_or_inf_impl<T, fkt> _isnaninf(stdisnaninf);
	return _isnaninf(mat);
}
// -----------------------------------------------------------------------------


/**
 * calculates the matrix inverse
 *
 * @desc code for inverse based on boost/libs/numeric/ublas/test/test_lu.cpp
 * @desc Boost's test_lu.cpp is (c) 2008 by G. Winkler
 */
template<class mat_type=ublas::matrix<double>>
bool inverse(const mat_type& mat, mat_type& inv)
{
	using T = typename mat_type::value_type;
	const typename mat_type::size_type N = mat.size1();
	if(N != mat.size2())
		return false;
	//if(isnan(mat) || isinf(mat))
	//	return false;

	try
	{
		mat_type lu = mat;
		ublas::permutation_matrix<typename mat_type::size_type> perm(N);

		if(ublas::lu_factorize(lu, perm) != 0)
			return false;

		inv = ublas::identity_matrix<T>(N);
		ublas::lu_substitute(lu, perm, inv);
	}
	catch(const std::exception& ex)
	{
		log_err("Matrix inversion failed with exception: ", ex.what(), ".", "\n",
				"Matrix to be inverted was: ", mat, ".");
		//log_backtrace();
		return false;
	}
	return true;
}


// R = T^(-1) M T
template<class mat_type=ublas::matrix<double>>
mat_type transform(const mat_type& mat, const mat_type& matTrafo, bool bOrtho=0)
{
	mat_type matTrafoInv;
	if(bOrtho)
		matTrafoInv = ublas::trans(matTrafo);
	else
		inverse(matTrafo, matTrafoInv);

	mat_type MT = ublas::prod(mat, matTrafo);
	mat_type TinvMT = ublas::prod(matTrafoInv, MT);

	return TinvMT;
}

// R = T M T^(-1)
template<class mat_type=ublas::matrix<double>>
mat_type transform_inv(const mat_type& mat, const mat_type& matTrafo, bool bOrtho=0)
{
	mat_type matTrafoInv;
	if(bOrtho)
		matTrafoInv = ublas::trans(matTrafo);
	else
		inverse(matTrafo, matTrafoInv);

	mat_type MT = ublas::prod(mat, matTrafoInv);
	mat_type TinvMT = ublas::prod(matTrafo, MT);

	return TinvMT;
}


template<typename T>
bool solve_linear_approx(const ublas::matrix<T>& M, const ublas::vector<T>& v,
				ublas::vector<T>& x);

// solve Mx = v for x
template<typename T=double>
bool solve_linear(const ublas::matrix<T>& M,
	const ublas::vector<T>& v, ublas::vector<T>& x)
{
	if(M.size1() == M.size2())		// determined, TODO: check rank
	{
		try
		{
			const std::size_t N = M.size1();

			ublas::matrix<T> lu = M;
			ublas::permutation_matrix<typename ublas::matrix<T>::size_type> perm(N);

			typename ublas::matrix<T>::size_type sing = ublas::lu_factorize(lu, perm);
			if(sing != 0)
				return false;

			x = v;
			ublas::lu_substitute(lu, perm, x);
		}
		catch(const std::exception& ex)
		{
			log_err("Linear equation solver failed with exception: ", ex.what(), ".");
			return false;
		}
	}
	else if(M.size1() < M.size2())	// underdetermined
	{
		ublas::matrix<T> Q, R;
		if(!qr(M, Q, R))
			return false;
		typedef typename ublas::vector<T>::size_type t_int;

		// M x = v
		// QR x = v
		// R x = Q^T v

		ublas::vector<T> vnew = ublas::prod(ublas::trans(Q), v);

		x = ublas::zero_vector<T>(M.size2());
		ublas::vector<T> xnew(R.size1());
		bool bOk = 0;

		// find non-singular right-upper submatrix
		std::vector<t_int> vecDelCols;
		std::size_t iNumToDel = R.size2()-R.size1();
		if(iNumToDel != 1)
		{
			log_err(__func__, " not yet implemented.");
			return false;
		}

		bool bFoundNonSingular = 0;
		ublas::matrix<T> Rsub;
		for(std::ptrdiff_t iCol=std::ptrdiff_t(R.size2()-1); iCol>=0; --iCol)
		{
			Rsub = remove_column(R, (std::size_t)iCol);

			T det = determinant<ublas::matrix<T>>(Rsub);
			if(!float_equal(det, 0.))
			{
				bFoundNonSingular = 1;
				vecDelCols.push_back(iCol);
				break;
			}
		}

		if(!bFoundNonSingular)
		{
			log_err("No non-singluar submatrix found in linear equation solver.");
			return false;
		}

		bOk = solve_linear(Rsub, vnew, xnew);

		for(t_int i=0, i0=0; i<xnew.size() && i0<x.size(); ++i, ++i0)
		{
			while(std::find(vecDelCols.begin(), vecDelCols.end(), i0) != vecDelCols.end())
				++i0;
			x[i0] = xnew[i];
		}

		return bOk;
	}
	else if(M.size1() > M.size2())	// overdetermined
		return solve_linear_approx<T>(M,v,x);
	else
		return false;

	return true;
}

// solve M^T M x = M^T v for x
template<typename T=double>
bool solve_linear_approx(const ublas::matrix<T>& M,
	const ublas::vector<T>& v, ublas::vector<T>& x)
{
	if(M.size1() <= M.size2())
	{
		//std::cerr << "Error: Matrix has to be overdetermined." << std::endl;
		return false;
	}

	ublas::matrix<T> Q, R;
	if(!qr(M, Q, R))
		return false;

	// M^T M x = M^T v
	// R^T Q^T Q R x = R^T Q^T v
	// R^T R x = R^T Q^T v

	const ublas::matrix<T> RT = ublas::trans(R);
	const ublas::matrix<T> QT = ublas::trans(Q);
	const ublas::matrix<T> RTR = ublas::prod(RT, R);
	const ublas::matrix<T> RTQT = ublas::prod(RT, QT);

	const ublas::vector<T> vnew = ublas::prod(RTQT, v);
	return solve_linear<T>(RTR, vnew, x);
}


template<class matrix_type=ublas::matrix<double>>
bool is_diag_matrix(const matrix_type& mat)
{
	typedef typename matrix_type::value_type T;

	for(std::size_t i=0; i<mat.size1(); ++i)
		for(std::size_t j=0; j<mat.size2(); ++j)
		{
			if(i==j) continue;

			if(!float_equal(mat(i,j), T(0.)))
				return false;
		}

	return true;
}


template<class matrix_type=ublas::matrix<double>,
	class vec_type=ublas::vector<typename matrix_type::value_type>,
	class container_type=std::initializer_list<vec_type>, const bool bRowMat>
inline matrix_type row_col_matrix(const container_type& vecs)
{
	if(vecs.size() == 0)
		return matrix_type(0,0);

	const std::size_t N = vecs.size();
	const std::size_t M = vecs.begin()->size();

	matrix_type mat(bRowMat?N:M, bRowMat?M:N);
	std::size_t j=0;
	for(typename container_type::const_iterator iter=vecs.begin(); iter!=vecs.end(); ++iter)
	{
		const vec_type& vec = *iter;

		for(std::size_t i=0; i<vec.size(); ++i)
		{
			if(bRowMat)
				mat(j,i) = vec[i];
			else
				mat(i,j) = vec[i];
		}

		++j;
	}

	return mat;
}

// vectors form rows of matrix
template<class matrix_type=ublas::matrix<double>,
	class vec_type=ublas::vector<typename matrix_type::value_type>,
	class container_type=std::initializer_list<vec_type>>
matrix_type row_matrix(const container_type& vecs)
{
	return row_col_matrix<matrix_type, vec_type, container_type, true>(vecs);
}

// vectors form columns of matrix
template<class matrix_type=ublas::matrix<double>,
	class vec_type=ublas::vector<typename matrix_type::value_type>,
	class container_type=std::initializer_list<vec_type>>
matrix_type column_matrix(const container_type& vecs)
{
	return row_col_matrix<matrix_type, vec_type, container_type, false>(vecs);
}

template<typename vector_type = ublas::vector<double>>
vector_type cross_3(const vector_type& vec0, const vector_type& vec1)
{
	return make_vec<vector_type>
		({
			vec0[1]*vec1[2] - vec1[1]*vec0[2],
			vec0[2]*vec1[0] - vec1[2]*vec0[0],
			vec0[0]*vec1[1] - vec1[0]*vec0[1]
		});
}


template<class t_mat/*=ublas::matrix<double>*/>
typename t_mat::value_type determinant(const t_mat& mat)
{
	typedef typename t_mat::value_type T;
	typedef typename t_mat::size_type t_size;

	if(mat.size1() != mat.size2())
		return T(0);

	if(mat.size1()==0)
		return T(0);
	else if(mat.size1()==1)
		return mat(0,0);
	else if(mat.size1()==2)
		return mat(0,0)*mat(1,1) - mat(1,0)*mat(0,1);
	else if(mat.size1()==3)
	{
		ublas::vector<T> vec0 = get_column(mat, 0);
		ublas::vector<T> vec1 = get_column(mat, 1);
		ublas::vector<T> vec2 = get_column(mat, 2);

		ublas::vector<T> vecCross = cross_3<ublas::vector<T> >(vec1, vec2);
		return ublas::inner_prod(vec0, vecCross);
	}
	else if(mat.size1()>3 && mat.size1()<6)		// recursive expansion, complexity: O(n!)
	{
		const t_size i = 0;
		T val = T(0);

		for(t_size j=0; j<mat.size2(); ++j)
		{
			if(float_equal<T>(mat(i,j), 0.))
				continue;

			T dSign = 1.;
			if(is_odd<std::size_t>(i+j))
				dSign = -1.;

			t_mat matSub = submatrix(mat, i, j);
			val += dSign * mat(i,j) * determinant<t_mat>(matSub);
		}

		return val;
	}
	else if(mat.size1()>=6)				// LU decomposition, complexity: O(n^3)
	{
		t_mat lu = mat;
		t_size N = mat.size1();
		ublas::permutation_matrix<typename t_mat::size_type> perm(N);

		ublas::lu_factorize(lu, perm);

		t_mat L = ublas::triangular_adaptor<t_mat, ublas::unit_lower>(lu);
		t_mat U = ublas::triangular_adaptor<t_mat, ublas::upper>(lu);

		T dDet = T(1.);
		for(t_size i=0; i<mat.size1(); ++i)
			dDet *= L(i,i)*U(i,i);

		std::size_t iNumSwaps=0;
		for(t_size iSwap=0; iSwap<perm.size(); ++iSwap)
			if(iSwap != perm(iSwap))
				++iNumSwaps;

		if(is_odd<std::size_t>(iNumSwaps))
			dDet *= T(-1.);

		return dDet;
	}

	return T(0);
}

template<class matrix_type=ublas::matrix<double>>
typename matrix_type::value_type get_volume(const matrix_type& mat)
{
	return determinant<matrix_type>(mat);
}


template<class matrix_type=ublas::matrix<double>>
typename matrix_type::value_type get_ellipsoid_volume(const matrix_type& mat)
{
	typedef typename matrix_type::value_type T;
	T tDet = determinant<matrix_type>(mat);

	return T(4./3.) * get_pi<T>() * std::sqrt(T(1)/tDet);
}



// calculate fractional coordinate basis vectors from angles
// see: http://www.bmsc.washington.edu/CrystaLinks/man/pdb/part_75.html
// for the reciprocal lattice this is equal to the B matrix from Acta Cryst. (1967), 22, 457
template<class t_vec>
bool fractional_basis_from_angles(typename t_vec::value_type a,
	typename t_vec::value_type b,
	typename t_vec::value_type c,
	typename t_vec::value_type alpha,
	typename t_vec::value_type beta,
	typename t_vec::value_type gamma,
	t_vec& veca, t_vec& vecb, t_vec& vecc)
{
	typedef typename t_vec::value_type T;

	const T dSG = std::sin(gamma);
	const T dCG = std::cos(gamma);
	const T dCA = std::cos(alpha);
	const T dCB = std::cos(beta);

	const T dCA2 = dCA*dCA;
	const T dCB2 = dCB*dCB;
	const T dCG2 = dCG*dCG;

	const T dVol =a*b*c*std::sqrt(1.- dCA2 - dCB2 - dCG2 + 2.*dCA*dCB*dCG);
	//std::cout << "vol = " <<  dVol << std::endl;

	if(std::isinf(dVol) || std::isnan(dVol))
		return false;

	veca[0] = a;
	veca[1] = 0.;
	veca[2] = 0.;

	vecb[0] = b*dCG;
	vecb[1] = b*dSG;
	vecb[2] = 0.;

	vecc[0] = c*dCB;
	vecc[1] = c*(dCA - dCB*dCG) / dSG;
	vecc[2] = dVol / (a*b*dSG);

	return true;
}


// signed angle wrt basis
template<typename vec_type>
typename vec_type::value_type vec_angle(const vec_type& vec)
{
	if(vec.size() == 2)
		return std::atan2(vec[1], vec[0]);

	throw Err("vec_angle not yet implemented for size != 2.");
}


// -----------------------------------------------------------------------------
template<typename T> void set_eps_0(T& d, underlying_value_type_t<T> eps=-1.);

template<typename T, bool bScalar=std::is_scalar<T>::value>
struct set_eps_0_impl
{
	void operator()(T&) const { throw Err("No implementation of set_eps_0!"); }
};

template<typename real_type>
struct set_eps_0_impl<real_type, 1>
{
	real_type eps = std::numeric_limits<real_type>::epsilon();

	void operator()(real_type& d) const
	{
		if(std::abs(d) < eps)
			d = real_type(0);
	}
};

template<typename vec_type>
struct set_eps_0_impl<vec_type, 0>
{
	using real_type = typename vec_type::value_type;
	real_type eps = std::numeric_limits<real_type>::epsilon();

	void operator()(vec_type& vec) const
	{
		for(real_type& d : vec)
			set_eps_0<real_type>(d, eps);
	}
};

template<typename T>
void set_eps_0(T& d, underlying_value_type_t<T> eps)
{
	set_eps_0_impl<T, std::is_scalar<T>::value> op;
	if(eps >= underlying_value_type_t<T>(0))
		op.eps = eps;
	op(d);
}
// -----------------------------------------------------------------------------

template<typename t_vec, typename T = typename t_vec::value_type>
bool vec_is_collinear(const t_vec& _vec1, const t_vec& _vec2, T eps = std::numeric_limits<T>::epsilon())
{
	const t_vec vec1 = _vec1 / ublas::norm_2(_vec1);
	const t_vec vec2 = _vec2 / ublas::norm_2(_vec2);

	T tdot = std::abs(ublas::inner_prod(vec1, vec2));
	return float_equal<T>(tdot, 1, eps);
}

// signed angle between two vectors
template<typename vec_type>
typename vec_type::value_type vec_angle(const vec_type& vec0,
	const vec_type& vec1, const vec_type* pvec_norm=0)
{
	typedef typename vec_type::value_type real_type;

	if(vec0.size() != vec1.size())
		throw Err("In vec_angle: Vector sizes do not match.");

	if(vec0.size() == 2)
	{
		return vec_angle<vec_type>(vec0) - vec_angle<vec_type>(vec1);
	}
	if(vec0.size() == 3)
	{
		real_type dC = ublas::inner_prod(vec0, vec1);
		vec_type veccross = cross_3<vec_type>(vec0, vec1);
		real_type dS = ublas::norm_2(veccross);

		real_type dAngle = std::atan2(dS, dC);

		// get signed angle
		if(pvec_norm)
		{
			if(ublas::inner_prod(veccross, *pvec_norm) < real_type(0))
				dAngle = -dAngle;
		}

		return dAngle;
	}

	throw Err("vec_angle only implemented for size == 2 and size == 3.");
}


template<class T, LinalgType ty=get_linalg_type<T>::value>
struct vec_angle_unsigned_impl
{
	void operator()(const T&, const T&) const { throw Err("No implementation of vec_angle_unsigned!"); }
};

// unsigned angle between two vectors
template<class T>
struct vec_angle_unsigned_impl<T, LinalgType::VECTOR>
{
	typename T::value_type operator()(const T& q1, const T& q2) const
	{
		typedef typename T::value_type REAL;

		if(q1.size() != q2.size())
			return REAL();

		REAL dot = REAL();
		REAL len1 = REAL();
		REAL len2 = REAL();
		for(std::size_t i=0; i<q1.size(); ++i)
		{
			dot += q1[i]*q2[i];

			len1 += q1[i]*q1[i];
			len2 += q2[i]*q2[i];
		}

		len1 = std::sqrt(len1);
		len2 = std::sqrt(len2);

		dot /= len1;
		dot /= len2;

		return std::acos(dot);
	}
};

template<class T>
typename T::value_type vec_angle_unsigned(const T& q1, const T& q2)
{
	return vec_angle_unsigned_impl<T>()(q1, q2);
}

// -----------------------------------------------------------------------------


// see: K. Shoemake, "Animating rotation with quaternion curves":
// http://dx.doi.org/10.1145/325334.325242
template<class T>
T slerp(const T& q1, const T& q2, typename T::value_type t)
{
	typedef typename T::value_type REAL;

	REAL angle = vec_angle_unsigned<T>(q1, q2);

	T q = std::sin((1.-t)*angle)/std::sin(angle) * q1 +
		std::sin(t*angle)/std::sin(angle) * q2;

	return q;
}



// --------------------------------------------------------------------------------


// see e.g.: http://www.itl.nist.gov/div898/handbook/pmc/section5/pmc541.htm
template<typename T=double>
ublas::matrix<T> covariance(const std::vector<ublas::vector<T>>& vecVals,
	const std::vector<T>* pProb = 0)
{
	if(vecVals.size() == 0) return ublas::matrix<T>();

	using t_vecvec = typename std::remove_reference<decltype(vecVals)>::type;
	using t_innervec_org = decltype(vecVals[0]);
	using t_innervec = typename std::remove_const<
						typename std::remove_reference<t_innervec_org>::type>
								::type;

	t_innervec vecMean = mean_value<t_vecvec>(vecVals);
	ublas::matrix<T> matCov(vecVals[0].size(), vecVals[0].size());

	T tSum = T(0);
	const std::size_t N = vecVals.size();
	for(std::size_t i=0; i<N; ++i)
	{
		T tprob = 1.;

		t_innervec vec = vecVals[i] - vecMean;
		ublas::matrix<T> matOuter = ublas::outer_prod(vec, vec);

		if(pProb)
		{
			tprob = (*pProb)[i];
			matOuter *= tprob;
		}

		matCov += matOuter;
		tSum += tprob;
	}
	matCov /= tSum;

	return matCov;
}




// --------------------------------------------------------------------------------
}

#include <boost/math/common_factor_rt.hpp>

namespace tl{

template<class t_vec=ublas::vector<int>>
t_vec get_gcd_vec(const t_vec& vec)
{
	if(vec.size() <= 1)
		return vec;

	typedef typename t_vec::value_type t_int;

	t_int igcd_total = 1;
	for(std::size_t i=0; i<vec.size()-1; ++i)
	{
		t_int i0 = vec[i];
		t_int i1 = vec[i+1];

		t_int igcd = boost::math::gcd<t_int>(i0, i1);

		if(i==0)
			igcd_total = igcd;
		else
			igcd_total = boost::math::gcd<t_int>(igcd, igcd_total);
	}

	if(igcd_total == 0)
		return vec;

	return vec/igcd_total;
}


// --------------------------------------------------------------------------------

// Householder reflection matrix
template<class t_mat = ublas::matrix<double>,
	class t_vec = ublas::vector<typename t_mat::value_type>,
	typename T = typename t_mat::value_type>
	t_mat reflection_matrix(const t_vec& vecNorm)
{
	t_mat mat = -T(2) * ublas::outer_prod(vecNorm, vecNorm);
	mat /= ublas::inner_prod(vecNorm, vecNorm);

	for(std::size_t i=0; i<vecNorm.size(); ++i)
		mat(i,i) += T(1);

	return mat;
}

// Householder reflection
template<class t_vec = ublas::vector<double>,
	class t_mat = ublas::matrix<typename t_vec::value_type>,
	typename T = typename t_mat::value_type>
	t_vec reflection(const t_vec& vec, const t_vec& vecNorm)
{
	t_mat mat = reflection_matrix<t_mat, t_vec, T>(vecNorm);
	return ublas::prod(mat, vec);
}

// add a nxn unit matrix to the upper left of a matrix
template<class t_mat = ublas::matrix<double>,
	typename T = typename t_mat::value_type>
	t_mat insert_unity(const t_mat& M, std::size_t n)
{
	if(M.size1()!=M.size2())
		throw Err("Non-square matrix not yet supported.");

	std::size_t m = M.size1();
	t_mat M2 = t_mat(m+n, m+n);

	for(std::size_t iR=0; iR<m+n; ++iR)
		for(std::size_t jR=0; jR<m+n; ++jR)
		{
			if(iR<n || jR<n)
				M2(iR, jR) = (iR==jR ? 1 : 0);
			else
				M2(iR, jR) = M(iR-n, jR-n);
		}

	return M2;
}

// QR decomposition via householder reflections
template<class t_mat = ublas::matrix<double>,
	class t_vec = ublas::vector<typename t_mat::value_type>,
	typename T = typename t_mat::value_type>
bool qr_decomp(const t_mat& M, t_mat& Q, t_mat& R)
{
	std::size_t m = M.size1();
	std::size_t n = M.size2();

	t_mat A = M;
	std::vector<t_mat> vecRefls;

	for(std::size_t i=0; i<std::min(m-1,n); ++i)
	{
		t_vec vec0 = get_column(A, 0);

		// vector of form [123.4 0 0 0] ?
		t_vec vec0_rest = ublas::subrange(vec0, 1, vec0.size());
		if(vec_equal<t_vec>(vec0_rest, ublas::zero_vector<T>(vec0_rest.size())))
		{
			t_mat matReflM = unit_matrix(m);
			vecRefls.push_back(matReflM);
			continue;
		}

		t_vec vecE0 = ublas::zero_vector<T>(vec0.size());
		vecE0[0] = ublas::norm_2(vec0);

		t_vec vecReflNorm = vec0-vecE0;
		//std::cout << "refl norm: " << vecReflNorm << std::endl;
		t_mat matRefl = reflection_matrix(vecReflNorm);

		A = ublas::prod(matRefl, A);
		A = submatrix(A,0,0);

		t_mat matReflM = insert_unity(matRefl, m-matRefl.size1());
		//std::cout << "refl: " << matReflM << std::endl;
		vecRefls.push_back(matReflM);
	}

	if(vecRefls.size() == 0)
		return false;

	Q = unit_matrix(m);
	for(const t_mat& matRefl : vecRefls)
	{
		t_mat matReflT = ublas::trans(matRefl);
		Q = ublas::prod(Q, matReflT);
	}

	t_mat QT = ublas::trans(Q);
	R = ublas::prod(QT, M);

	return true;
}


template<typename t_vec = ublas::vector<double>,
	typename T = typename t_vec::value_type>
std::vector<t_vec> gram_schmidt(const std::vector<t_vec>& vecs, bool bNorm=true);

// QR decomposition via gram-schmidt orthogonalisation
template<class t_mat = ublas::matrix<double>,
	class t_vec = ublas::vector<typename t_mat::value_type>,
	typename T = typename t_mat::value_type>
bool qr_decomp_gs(const t_mat& M, t_mat& Q, t_mat& R)
{
	Q = column_matrix(gram_schmidt(get_columns(M), 1));

	// M = QR  =>  Q^T M = R
	R = ublas::prod(ublas::trans(Q), M);
	return 1;
}


template<class t_mat = ublas::matrix<double>,
	class t_vec = ublas::vector<typename t_mat::value_type>,
	typename T = typename t_mat::value_type>
t_mat norm_col_vecs(const t_mat& M)
{
	t_mat N(M.size1(), M.size2());

	for(std::size_t i=0; i<M.size2(); ++i)
	{
		t_vec vec0 = get_column(M, i);
		vec0 /= ublas::norm_2(vec0);

		set_column(N, i, vec0);
	}

	return N;
}

template<class t_mat=ublas::matrix<double>, class t_real=underlying_value_type_t<t_mat>>
bool is_symmetric(const t_mat& mat, t_real eps=std::numeric_limits<t_real>::epsilon())
{
	if(mat.size1() != mat.size2())
		return false;

	for(std::size_t i=0; i<mat.size1(); ++i)
		for(std::size_t j=i+1; i<mat.size2(); ++i)
			if(!float_equal(mat(i,j), mat(j,i), eps))
				return false;

	return true;
}


// -----------------------------------------------------------------------------
template<class T, LinalgType ty=get_linalg_type<T>::value>
struct apply_fkt_impl
{
	using value_type = underlying_value_type_t<T>;

	T operator()(T, const std::function<value_type(value_type)>& fkt) const
	{
		throw Err("No implementation of apply_fkt!");
	}
};

template<class T>
struct apply_fkt_impl<T, LinalgType::REAL>
{
	T operator()(T t, const std::function<T(T)>& fkt) const
	{
		return fkt(t);
	}
};

template<class t_vec>
struct apply_fkt_impl<t_vec, LinalgType::VECTOR>
{
	using value_type = underlying_value_type_t<t_vec>;

	t_vec operator()(const t_vec& vec, const std::function<value_type(value_type)>& fkt) const
	{
		t_vec v;
		v.resize(vec.size());

		for(std::size_t i=0; i<vec.size(); ++i)
			v[i] = fkt(vec[i]);

		return v;
	}
};

template<class t_mat>
struct apply_fkt_impl<t_mat, LinalgType::MATRIX>
{
	using value_type = underlying_value_type_t<t_mat>;

	t_mat operator()(const t_mat& mat, const std::function<value_type(value_type)>& fkt) const
	{
		t_mat m;
		m.resize(mat.size1(), mat.size2());

		for(std::size_t i=0; i<mat.size1(); ++i)
			for(std::size_t j=0; j<mat.size2(); ++j)
				m(i,j) = fkt(mat(i,j));

		return m;
	}
};

template<class T, class t_val=underlying_value_type_t<T>>
T apply_fkt(const T& t, const std::function<t_val(t_val)>& fkt)
{
	apply_fkt_impl<T> impl;
	return impl(t, fkt);
}

template<class T, class t_val=underlying_value_type_t<T>>
inline T apply_fkt(const T& t, t_val(*pfkt)(t_val))
{
	std::function<t_val(t_val)> fkt(pfkt);
	return apply_fkt<T, t_val>(t, fkt);
}
// -----------------------------------------------------------------------------


template<class T, LinalgType ty=get_linalg_type<T>::value>
struct get_minmax_impl
{
	using value_type = underlying_value_type_t<T>;

	std::pair<value_type, value_type>
	operator()(T) const
	{
		throw Err("No implementation of get_minmax!");
	}
};

template<class T>
struct get_minmax_impl<T, LinalgType::REAL>
{
	std::pair<T, T>
	operator()(T t) const
	{
		return std::pair<T,T>(t,t);
	}
};

template<class t_vec>
struct get_minmax_impl<t_vec, LinalgType::VECTOR>
{
	using t_val = underlying_value_type_t<t_vec>;

	std::pair<t_val, t_val>
	operator()(const t_vec& vec) const
	{
		t_val tmin = std::numeric_limits<t_val>::max();
		t_val tmax = -tmin;

		for(std::size_t i=0; i<vec.size(); ++i)
		{
			if(vec[i] < tmin) tmin = vec[i];
			if(vec[i] > tmax) tmax = vec[i];
		}

		return std::pair<t_val, t_val>(tmin, tmax);
	}
};

template<class t_mat>
struct get_minmax_impl<t_mat, LinalgType::MATRIX>
{
	using t_val = underlying_value_type_t<t_mat>;

	std::pair<t_val, t_val>
	operator()(const t_mat& mat) const
	{
		t_val tmin = std::numeric_limits<t_val>::max();
		t_val tmax = -tmin;

		for(std::size_t i=0; i<mat.size1(); ++i)
			for(std::size_t j=0; j<mat.size2(); ++j)
			{
				if(mat(i,j) < tmin) tmin = mat(i,j);
				if(mat(i,j) > tmax) tmax = mat(i,j);
			}

		return std::pair<t_val, t_val>(tmin, tmax);
	}
};

template<class T>
std::pair<underlying_value_type_t<T>, underlying_value_type_t<T>>
get_minmax(const T& t)
{
	get_minmax_impl<T> impl;
	return impl(t);
}
// -----------------------------------------------------------------------------


// ! for large matrices use eigenvec_sym from linalg2.h !
template<class t_mat = ublas::matrix<double>,
	class t_vec = ublas::vector<typename t_mat::value_type>,
	typename T = typename t_mat::value_type>
bool eigenvec_sym_simple(const t_mat& mat, std::vector<t_vec>& evecs, std::vector<T>& evals)
{
	if(mat.size1() != mat.size2())
	{
		log_err("Matrix ", mat, " is not square.");
		return false;
	}

//#ifndef NDEBUG
	t_mat matAbs = apply_fkt(mat, std::function<T(T)>((T(*)(T))std::abs));
	T _dEps = get_minmax(matAbs).second / 100.;	// 1% accuracy
	if(!tl::is_symmetric(mat, _dEps)) log_warn("Matrix ", mat, " is not symmetric.");
//#endif

	const std::size_t n = mat.size1();
	t_mat I = ublas::identity_matrix<T>(n);
	t_mat M = mat;

	const T tEps = std::cbrt(std::numeric_limits<T>::epsilon());
	const std::size_t MAX_ITER = 512;
	std::size_t iIter = 0;
	for(iIter=0; iIter<MAX_ITER; ++iIter)
	{
		t_mat Q, R;
		if(!qr_decomp(M, Q, R))
		{
			log_err("QR decomposition failed for matrix ", M);
			return false;
		}

		t_mat Mlast = M;
		M = ublas::prod(R, Q);
		I = ublas::prod(I, Q);


		bool bConverged = 1;
		for(std::size_t iVal=0; iVal<n; ++iVal)
		{
			if(std::abs(M(iVal,iVal)-Mlast(iVal,iVal)) > tEps)
			{
				bConverged = 0;
				break;
			}
		}

		if(bConverged)
			break;
	}

	evals.resize(n);
	evecs.resize(n);

	for(std::size_t iVal=0; iVal<n; ++iVal)
	{
		evals[iVal] = M(iVal, iVal);
		evecs[iVal] = get_column(I, iVal);
	}

	return true;
}

template<typename T=double>
void sort_eigenvecs(std::vector<ublas::vector<T> >& evecs,
	std::vector<T>& evals, bool bOrder=0, T (*pEvalFkt)(T)=0)
{
	if(evecs.size() != evals.size())
		return;

	struct Evec
	{
		ublas::vector<T> vec;
		T val;
	};

	std::vector<Evec> myevecs;
	myevecs.reserve(evecs.size());

	for(std::size_t i=0; i<evecs.size(); ++i)
	{
		Evec ev;
		ev.vec = evecs[i];
		ev.val = evals[i];

		myevecs.push_back(ev);
	}


	std::sort(myevecs.begin(), myevecs.end(),
		[&](const Evec& evec1, const Evec& evec2) -> bool
		{
			bool b;
			if(pEvalFkt)
				b = pEvalFkt(evec1.val) < pEvalFkt(evec2.val);
			else
				b = evec1.val < evec2.val;

			if(bOrder) b = !b;
			return b;
		});


	for(std::size_t i=0; i<evecs.size(); ++i)
	{
		evecs[i] = myevecs[i].vec;
		evals[i] = myevecs[i].val;
	}
}


// --------------------------------------------------------------------------------


template<typename t_vec /*= ublas::vector<double>*/,
	typename T /*= typename t_vec::value_type*/ >
std::vector<t_vec> gram_schmidt(const std::vector<t_vec>& vecs, bool bNorm/*=1*/)
{
	std::vector<t_vec> vecsOut;
	if(vecs.size() == 0)
		return vecsOut;

	vecsOut.resize(vecs.size());
	for(std::size_t i=0; i<vecs.size(); ++i)
	{
		vecsOut[i] = vecs[i];
		for(std::size_t j=0; j<i; ++j)
		{
			T tnum = ublas::inner_prod(vecs[i], vecsOut[j]);
			T tden = ublas::inner_prod(vecsOut[j], vecsOut[j]);

			vecsOut[i] -= tnum/tden * vecsOut[j];
		}
	}

	if(bNorm)
		for(t_vec& vec : vecsOut)
			vec /= ublas::norm_2(vec);

	return vecsOut;
}

template<typename t_vec = ublas::vector<double>,
	typename T = typename t_vec::value_type>
std::vector<t_vec> get_ortho_rhs(const std::vector<t_vec>& vecs)
{
	assert(vecs.size() == 2);

	std::vector<t_vec> vecOrtho = gram_schmidt(vecs, true);
	t_vec vecUp = cross_3(vecOrtho[0], vecOrtho[1]);
	vecOrtho.push_back(vecUp);

	return vecOrtho;
}

}

#endif
