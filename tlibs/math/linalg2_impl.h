/**
 * advanced linalg helpers
 *
 * @author: Tobias Weber
 * @date: 2013-2016
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIBS_LINALG2_IMPL_H__
#define __TLIBS_LINALG2_IMPL_H__

#include "linalg2.h"
#include <memory>

#ifndef NO_LAPACK
extern "C"
{
	#define lapack_complex_float std::complex<float>
	#define lapack_complex_float_real(c) (c.real())
	#define lapack_complex_float_imag(c) (c.imag())
	#define lapack_complex_double std::complex<double>
	#define lapack_complex_double_real(c) (c.real())
	#define lapack_complex_double_imag(c) (c.imag())

	#include <lapacke.h>
}

namespace tl {

// selects the float or double version of a lapack function
template<class T1, class T2, class F1, class F2>
struct select_func
{
	F1* m_f1 = nullptr;
	F2* m_f2 = nullptr;

	select_func(F1* f1, F2* f2) : m_f1(f1), m_f2(f2) {}

	template<class T>
	typename std::enable_if<std::is_same<T, T1>::value, F1*>::type 
		get_func() { return m_f1; }
	template<class T>
	typename std::enable_if<std::is_same<T, T2>::value, F2*>::type 
		get_func() { return m_f2; }
};


template<class T>
bool eigenvec(const ublas::matrix<T>& mat,
	std::vector<ublas::vector<T>>& evecs_real,
	std::vector<ublas::vector<T>>& evecs_imag,
	std::vector<T>& evals_real,
	std::vector<T>& evals_imag)
{
	bool bOk = true;
	select_func<float, double, decltype(LAPACKE_sgeev), decltype(LAPACKE_dgeev)>
		sfunc(LAPACKE_sgeev, LAPACKE_dgeev);
	auto pfunc = sfunc.get_func<T>();

	if(mat.size1() != mat.size2())
		return false;
	if(mat.size1()==0 || mat.size1()==1)
		return false;

	const std::size_t iOrder = mat.size1();
	evecs_real.resize(iOrder); evecs_imag.resize(iOrder);
	evals_real.resize(iOrder); evals_imag.resize(iOrder);
	for(std::size_t i=0; i<iOrder; ++i)
	{
		evecs_real[i].resize(iOrder);
		evecs_imag[i].resize(iOrder);
	}

	std::unique_ptr<T, std::default_delete<T[]>>
		uptrMem(new T[iOrder*iOrder + iOrder*iOrder + iOrder*iOrder]);
	T *pMatrix = uptrMem.get();
	T *pEVs = pMatrix + iOrder*iOrder;

	for(std::size_t i=0; i<iOrder; ++i)
		for(std::size_t j=0; j<iOrder; ++j)
			pMatrix[i*iOrder + j] = mat(i,j);

	int iInfo = (*pfunc)(LAPACK_ROW_MAJOR, 'N', 'V', iOrder,
		pMatrix, iOrder, evals_real.data(), evals_imag.data(),
		nullptr, iOrder, pEVs, iOrder);

	if(iInfo!=0)
	{
		log_err("Could not solve general real eigenproblem", 
			" (lapack error ", iInfo , ").");
		bOk = false;
	}

	for(std::size_t i=0; i<iOrder; ++i)
	{
		bool bIsReal = 0;
		if(float_equal<T>(evals_imag[i], 0.))
			bIsReal = 1;

		if(bIsReal)
		{
			for(std::size_t j=0; j<iOrder; ++j)
			{
				evecs_real[i][j] = pEVs[j*iOrder + i];
				evecs_imag[i][j] = 0.;
			}
		}
		else
		{
			for(std::size_t j=0; j<iOrder; ++j)
			{
				evecs_real[i][j] = pEVs[j*iOrder + i];
				evecs_imag[i][j] = pEVs[j*iOrder + i+1];

				evecs_real[i+1][j] = pEVs[j*iOrder + i];
				evecs_imag[i+1][j] = -pEVs[j*iOrder + i+1];
			}
			++i; // check: (next eigenval) == -(currrent eigenval)
		}

		//evecs_real[i] /= ublas::norm_2(evecs_real[i]);
		//evecs_imag[i] /= ublas::norm_2(evecs_imag[i]);
	}
	return bOk;
}

template<class T>
bool eigenvec_cplx(const ublas::matrix<std::complex<T>>& mat,
	std::vector<ublas::vector<std::complex<T>>>& evecs,
	std::vector<std::complex<T>>& evals)
{
	using t_cplx = std::complex<T>;
	bool bOk = true;
	select_func<float, double, decltype(LAPACKE_cgeev), decltype(LAPACKE_zgeev)>
		sfunc(LAPACKE_cgeev, LAPACKE_zgeev);
	auto pfunc = sfunc.get_func<T>();

	if(mat.size1() != mat.size2())
		return false;
	if(mat.size1()==0 || mat.size1()==1)
		return false;

	const std::size_t iOrder = mat.size1();
	evecs.resize(iOrder);
	evals.resize(iOrder);
	for(std::size_t i=0; i<iOrder; ++i)
		evecs[i].resize(iOrder);

	std::unique_ptr<t_cplx, std::default_delete<t_cplx[]>>
		uptrMem(new t_cplx[iOrder*iOrder + iOrder*iOrder + iOrder]);
	t_cplx *pMatrix = uptrMem.get();

	for(std::size_t i=0; i<iOrder; ++i)
		for(std::size_t j=0; j<iOrder; ++j)
			pMatrix[i*iOrder + j] = mat(i,j);

	t_cplx *pEVs = pMatrix + iOrder*iOrder;
	t_cplx *pEVals = pEVs + iOrder*iOrder;

	int iInfo = (*pfunc)(LAPACK_ROW_MAJOR, 'N', 'V', iOrder,
		pMatrix, iOrder, pEVals,
		nullptr, iOrder, pEVs, iOrder);

	if(iInfo!=0)
	{
		log_err("Could not solve general complex eigenproblem", 
			" (lapack error ", iInfo , ").");
		bOk = false;
	}

	for(std::size_t i=0; i<iOrder; ++i)
	{
		for(std::size_t j=0; j<iOrder; ++j)
			evecs[i][j] = pEVs[j*iOrder + i];
		evals[i] = pEVals[i];
	}

	return bOk;
}


template<class T>
bool eigenvec_sym(const ublas::matrix<T>& mat,
	std::vector<ublas::vector<T>>& evecs,
	std::vector<T>& evals)
{
	bool bOk = true;
	select_func<float, double, decltype(LAPACKE_ssyev), decltype(LAPACKE_dsyev)> 
		sfunc(LAPACKE_ssyev, LAPACKE_dsyev);
	auto pfunc = sfunc.get_func<T>();

	if(mat.size1() != mat.size2())
		return false;
	if(mat.size1()==0 || mat.size1()==1)
		return false;

	const std::size_t iOrder = mat.size1();
	evecs.resize(iOrder);
	evals.resize(iOrder);
	for(std::size_t i=0; i<iOrder; ++i)
		evecs[i].resize(iOrder);

	std::unique_ptr<T, std::default_delete<T[]>> 
		uptrMat(new T[iOrder*iOrder]);
	T *pMatrix = uptrMat.get();

	for(std::size_t i=0; i<iOrder; ++i)
		for(std::size_t j=0; j<iOrder; ++j)
			pMatrix[i*iOrder + j] = mat(i,j);

	int iInfo = (*pfunc)(LAPACK_ROW_MAJOR, 'V', 'U',
		iOrder, pMatrix, iOrder, evals.data());

	if(iInfo!=0)
	{
		log_err("Could not solve symmetric eigenproblem",
			" (lapack error ", iInfo, ").");
		bOk = false;
	}

	for(std::size_t i=0; i<iOrder; ++i)
	{
		for(std::size_t j=0; j<iOrder; ++j)
			evecs[i][j] = pMatrix[j*iOrder + i];
		//evecs[i] /= ublas::norm_2(evecs[i]);
	}

	if(determinant<ublas::matrix<T>>(column_matrix(evecs)) < 0.)
		evecs[0] = -evecs[0];

	return bOk;
}


template<class T>
bool eigenvec_herm(const ublas::matrix<std::complex<T>>& mat,
	std::vector<ublas::vector<std::complex<T>>>& evecs,
	std::vector<T>& evals)
{
	using t_cplx = std::complex<T>;
	bool bOk = true;

	select_func<float, double, decltype(LAPACKE_cheev), decltype(LAPACKE_zheev)> 
		sfunc(LAPACKE_cheev, LAPACKE_zheev);
	auto pfunc = sfunc.get_func<T>();

	if(mat.size1() != mat.size2())
		return false;
	if(mat.size1()==0 || mat.size1()==1)
		return false;

	const std::size_t iOrder = mat.size1();
	evecs.resize(iOrder);
	evals.resize(iOrder);
	for(std::size_t i=0; i<iOrder; ++i)
		evecs[i].resize(iOrder);

	std::unique_ptr<t_cplx, std::default_delete<t_cplx[]>>
		uptrMat(new t_cplx[iOrder*iOrder]);
	t_cplx *pMatrix = uptrMat.get();

	for(std::size_t i=0; i<iOrder; ++i)
		for(std::size_t j=0; j<iOrder; ++j)
			pMatrix[i*iOrder + j] = mat(i,j);

	int iInfo = (*pfunc)(LAPACK_ROW_MAJOR, 'V', 'U',
		iOrder, pMatrix, iOrder, evals.data());

	if(iInfo!=0)
	{
		log_err("Could not solve hermitian eigenproblem",
			" (lapack error ", iInfo, ").");
		bOk = false;
	}

	for(std::size_t i=0; i<iOrder; ++i)
		for(std::size_t j=0; j<iOrder; ++j)
			evecs[i][j] = pMatrix[j*iOrder + i];
	return bOk;
}


template<class T>
bool qr(const ublas::matrix<T>& M,
	ublas::matrix<T>& Q, ublas::matrix<T>& R)
{
	select_func<float, double, decltype(LAPACKE_sgeqrf), decltype(LAPACKE_dgeqrf)> 
		sfunc(LAPACKE_sgeqrf, LAPACKE_dgeqrf);
	auto pfunc = sfunc.get_func<T>();

	const typename ublas::matrix<T>::size_type m = M.size1();
	const typename ublas::matrix<T>::size_type n = M.size2();

	const std::size_t iTauSize = m;//std::min<std::size_t>(m,n);

	std::unique_ptr<T, std::default_delete<T[]>>
		uptrMem(new T[n*m + iTauSize]);
	T *pMem = uptrMem.get();

	T *pMat = pMem;
	T *pTau = pMem + n*m;

	for(std::size_t i=0; i<m; ++i)
		for(std::size_t j=0; j<n; ++j)
			pMat[i*n + j] = M(i,j);

	// see: http://www.math.utah.edu/software/lapack/lapack-d/dgeqrf.html
	int iInfo = (*pfunc)(LAPACK_ROW_MAJOR, m, n, pMat, n, pTau);

	R = ublas::matrix<T>(m,n);
	for(std::size_t i=0; i<m; ++i)
		for(std::size_t j=0; j<n; ++j)
		{
			if(j>=i)
				R(i,j) = pMat[i*n + j];
			else
				R(i,j) = 0.;
		}

	ublas::vector<T> v(iTauSize);

	const ublas::matrix<T> ident = ublas::identity_matrix<T>(iTauSize);
	Q = ident;

	for(std::size_t k=1; k<=iTauSize; ++k)
	{
		T dTau = pTau[k-1];

		for(std::size_t i=1; i<=k-1; ++i)
			v[i-1] = 0.;
		v[k-1] = 1.;

		for(std::size_t i=k+1; i<=iTauSize; ++i)
			v[i-1] = pMat[(i-1)*n + (k-1)];

		ublas::matrix<T> VV = ublas::outer_prod(v, ublas::trans(v));
		ublas::matrix<T> H = ident - dTau*VV;

		Q = ublas::prod(Q, H);
	}

	return (iInfo==0);
}

}

#endif

#endif
