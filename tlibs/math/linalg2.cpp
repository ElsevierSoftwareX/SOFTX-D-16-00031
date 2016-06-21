/*
 * advanced linalg helpers (which depend on lapack/e)
 *
 * @author: tweber
 * @date: 30-apr-2013
 * @license GPLv2 or GPLv3
 */

#include "linalg2.h"
#include "linalg2_impl.h"

#ifndef NO_LAPACK

namespace tl {

template bool qr(const ublas::matrix<double>& M,
	ublas::matrix<double>& Q, ublas::matrix<double>& R);


template bool eigenvec(const ublas::matrix<double>& mat,
	std::vector<ublas::vector<double> >& evecs_real, std::vector<ublas::vector<double>>& evecs_imag,
	std::vector<double>& evals_real, std::vector<double>& evals_imag);

template bool eigenvec_cplx(const ublas::matrix<std::complex<double>>& mat,
	std::vector<ublas::vector<std::complex<double>> >& evecs,
	std::vector<std::complex<double>>& evals);


template bool eigenvec_sym(const ublas::matrix<double>& mat,
	std::vector<ublas::vector<double>>& evecs, std::vector<double>& evals);

template bool eigenvec_herm(const ublas::matrix<std::complex<double>>& mat,
	std::vector<ublas::vector<std::complex<double>>>& evecs,
	std::vector<double>& evals);

}

#endif
