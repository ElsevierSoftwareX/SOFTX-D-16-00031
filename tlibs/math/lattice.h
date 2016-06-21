/**
 * Bravais Lattice Calculations
 * @author tweber
 * @date 13-feb-2014
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIBS_LATTICE_H__
#define __TLIBS_LATTICE_H__

#include "linalg.h"
#include "math.h"
#include "neutrons.hpp"
#include <ostream>

namespace tl {

template<typename T=double>
bool reciprocal(const ublas::matrix<T>& matReal, ublas::matrix<T>& matRecip)
{
	ublas::matrix<T> matInv;
	if(!inverse<ublas::matrix<T>>(ublas::trans(matReal), matInv))
		return false;

	matRecip = T(2)*get_pi<T>()*matInv;
	return true;
}


template<typename T=double>
class Lattice
{
	public:
		using t_vec = ublas::vector<T>;
		using t_mat = ublas::matrix<T>;

	protected:
		t_vec m_vecs[3];

	public:
		Lattice(T a, T b, T c, T alpha, T beta, T gamma);
		Lattice(const t_vec& vec0, const t_vec& vec1, const t_vec& vec2);
		Lattice(const Lattice<T>& lattice);
		Lattice() = default;
		~Lattice() = default;

		bool IsInited() const { return m_vecs[0].size()!=0; }

		// Euler ZXZ rotation
		void RotateEuler(T dPhi, T dTheta, T dPsi);
		// Euler vecRecipZ vecRecipX vecRecipZ rotation
		void RotateEulerRecip(const t_vec& vecRecipX,
			const t_vec& vecRecipY, const t_vec& vecRecipZ,
			T dPhi, T dTheta, T dPsi);

		Lattice GetRecip() const;
		Lattice GetAligned() const;

		t_vec GetPos(T h, T k, T l) const;
		t_vec GetHKL(const t_vec& vec) const;

		T GetAlpha() const;
		T GetBeta() const;
		T GetGamma() const;

		T GetA() const;
		T GetB() const;
		T GetC() const;

		T GetVol() const;

		const t_vec& GetVec(std::size_t i) const { return m_vecs[i]; }
		t_mat GetMetric() const;
};

template<typename T>
Lattice<T>::Lattice(T a, T b, T c, T alpha, T beta, T gamma)
{
	m_vecs[0].resize(3,0);
	m_vecs[1].resize(3,0);
	m_vecs[2].resize(3,0);

	fractional_basis_from_angles(a,b,c, alpha,beta,gamma, m_vecs[0],m_vecs[1],m_vecs[2]);
}

template<typename T>
Lattice<T>::Lattice(const t_vec& vec0, const t_vec& vec1, const t_vec& vec2)
{
	this->m_vecs[0] = vec0;
	this->m_vecs[1] = vec1;
	this->m_vecs[2] = vec2;
}

template<typename T>
Lattice<T>::Lattice(const Lattice<T>& lattice)
{
	this->m_vecs[0] = lattice.m_vecs[0];
	this->m_vecs[1] = lattice.m_vecs[1];
	this->m_vecs[2] = lattice.m_vecs[2];
}

template<typename T>
void Lattice<T>::RotateEuler(T dPhi, T dTheta, T dPsi)
{
	t_mat mat1 = rotation_matrix_3d_z(dPhi);
	t_mat mat2 = rotation_matrix_3d_x(dTheta);
	t_mat mat3 = rotation_matrix_3d_z(dPsi);

	t_mat mat21 = ublas::prod(mat2,mat1);
	t_mat mat = ublas::prod(mat3, mat21);

	for(std::size_t i=0; i<3; ++i)
		m_vecs[i] = ublas::prod(mat, m_vecs[i]);
}

template<typename T>
void Lattice<T>::RotateEulerRecip(const t_vec& vecRecipX,
	const t_vec& vecRecipY, const t_vec& vecRecipZ,
	T dPhi, T dTheta, T dPsi)
{
	// get real vectors
	const std::size_t iDim=3;
	t_mat matReal = column_matrix({vecRecipX, vecRecipY, vecRecipZ});
	if(matReal.size1()!=matReal.size2() || matReal.size1()!=iDim)
		throw Err("Invalid real matrix.");

	t_mat matRecip;
	if(!reciprocal(matReal, matRecip))
		throw Err("Reciprocal matrix could not be calculated.");

	t_vec vecX = get_column(matRecip,0);
	t_vec vecY = get_column(matRecip,1);
	t_vec vecZ = get_column(matRecip,2);

	T dLenX = ublas::norm_2(vecX);
	T dLenY = ublas::norm_2(vecY);
	T dLenZ = ublas::norm_2(vecZ);

	if(float_equal(dLenX, 0.) || float_equal(dLenY, 0.) || float_equal(dLenZ, 0.)
		|| std::isnan(dLenX) || std::isnan(dLenY) || std::isnan(dLenZ))
	{
		throw Err("Invalid reciprocal matrix.");
	}

	vecX /= dLenX;
	vecY /= dLenY;
	vecZ /= dLenZ;

	// rotate around real vectors
	t_mat mat1 = rotation_matrix(vecZ, dPhi);
	t_mat mat2 = rotation_matrix(vecX, dTheta);
	t_mat mat3 = rotation_matrix(vecZ, dPsi);

	t_mat mat21 = ublas::prod(mat2,mat1);
	t_mat mat = ublas::prod(mat3, mat21);

	for(std::size_t i=0; i<3; ++i)
		m_vecs[i] = ublas::prod(mat, m_vecs[i]);
}

template<typename T> T Lattice<T>::GetAlpha() const
{ return std::acos(ublas::inner_prod(m_vecs[1]/GetB(), m_vecs[2]/GetC())); }
template<typename T> T Lattice<T>::GetBeta() const
{ return std::acos(ublas::inner_prod(m_vecs[0]/GetA(), m_vecs[2]/GetC())); }
template<typename T> T Lattice<T>::GetGamma() const
{ return std::acos(ublas::inner_prod(m_vecs[0]/GetA(), m_vecs[1]/GetB())); }

template<typename T> T Lattice<T>::GetA() const { return ublas::norm_2(m_vecs[0]); }
template<typename T> T Lattice<T>::GetB() const { return ublas::norm_2(m_vecs[1]); }
template<typename T> T Lattice<T>::GetC() const { return ublas::norm_2(m_vecs[2]); }

template<typename T>
T Lattice<T>::GetVol() const
{
	return get_volume(column_matrix({m_vecs[0], m_vecs[1], m_vecs[2]}));
}

/*
 (x)   (v0_x v1_x v2_x) (h)
 (y) = (v0_y v1_y v2_y) (k)
 (z)   (v0_z v1_z v2_z) (l)
 */
template<typename T>
typename Lattice<T>::t_vec Lattice<T>::GetPos(T h, T k, T l) const
{
	return h*m_vecs[0] + k*m_vecs[1] + l*m_vecs[2];
}

/*
 (h)   (v0_x v1_x v2_x)^(-1) (x)
 (k) = (v0_y v1_y v2_y)      (y)
 (l)   (v0_z v1_z v2_z)      (z)
 */
template<typename T>
typename Lattice<T>::t_vec Lattice<T>::GetHKL(const t_vec& vec) const
{
	t_mat mat = column_matrix({m_vecs[0], m_vecs[1], m_vecs[2]});

	t_mat matInv;
	if(!inverse(mat, matInv))
		throw Err("Miller indices could not be calculated.");

	return ublas::prod(matInv, vec);
}

template<typename T>
Lattice<T> Lattice<T>::GetRecip() const
{
	const std::size_t iDim=3;
	t_mat matReal = column_matrix({m_vecs[0], m_vecs[1], m_vecs[2]});
	if(matReal.size1()!=matReal.size2() || matReal.size1()!=iDim)
		throw Err("Invalid real lattice matrix.");

	t_mat matRecip;
	if(!reciprocal(matReal, matRecip))
		throw Err("Reciprocal lattice could not be calculated.");

	// warning: first axis does not (necessarily) coincide with assumed first
	//			orientation vector [0,0,1] anymore!
	return Lattice<T>(get_column(matRecip,0), get_column(matRecip,1),
		get_column(matRecip,2));
}

template<typename T>
Lattice<T> Lattice<T>::GetAligned() const
{
	// construct new, correctly oriented reciprocal lattice with first axis along
	// [0,0,1]
	return Lattice<T>(GetA(), GetB(), GetC(), GetAlpha(), GetBeta(), GetGamma());
}


template<typename T>
ublas::matrix<T> Lattice<T>::GetMetric() const
{
	return column_matrix({m_vecs[0], m_vecs[1], m_vecs[2]});
}


// -----------------------------------------------------------------------------

/**
 * B matrix converts rlu to 1/A
 */
template<typename T=double>
ublas::matrix<T> get_B(const Lattice<T>& lattice, bool bIsRealLattice=1)
{
	using t_mat = ublas::matrix<T>;

	t_mat matB;
	if(bIsRealLattice)
		matB = lattice.GetRecip()/*.GetAligned()*/.GetMetric();
	else
		matB = lattice/*.GetAligned()*/.GetMetric();

	return matB;
}


/**
 * U matrix expresses the coordinates in the basis of the scattering plane
 */
template<typename T=double>
ublas::matrix<T> get_U(const ublas::vector<T>& _vec1, const ublas::vector<T>& _vec2,
	const ublas::matrix<T>* pmatB=nullptr)
{
	using t_vec = ublas::vector<T>;
	using t_mat = ublas::matrix<T>;

	t_vec vec1;
	t_vec vec2;

	if(pmatB)
	{
		// in 1/A
		vec1 = ublas::prod(*pmatB, _vec1);
		vec2 = ublas::prod(*pmatB, _vec2);
	}
	else
	{
		// in rlu
		vec1 = _vec1;
		vec2 = _vec2;
	}

	// U: scattering plane coordinate system
	t_mat matU = row_matrix(get_ortho_rhs({vec1, vec2}));
	return matU;
}


/**
 * UB matrix converts rlu to 1/A and expresses it in the scattering plane coords:
 * Q = U*B*hkl
 */
template<typename T=double>
ublas::matrix<T> get_UB(const Lattice<T>& lattice_real,
	const ublas::vector<T>& _vec1, const ublas::vector<T>& _vec2)
{
	using t_mat = ublas::matrix<T>;

	t_mat matB = get_B(lattice_real, 1);		// rlu to 1/A
	t_mat matU = get_U(_vec1, _vec2, &matB);	// scattering in 1/A

	t_mat matUB = ublas::prod(matU, matB);
	return matUB;
}

// -----------------------------------------------------------------------------


// distance for point to be considered inside scattering plane
template<typename T> constexpr T get_plane_dist_tolerance() { return T(1e-6); }

template<typename T=double>
void get_tas_angles(const Lattice<T>& lattice_real,
	const ublas::vector<T>& _vec1, const ublas::vector<T>& _vec2,
	T dKi, T dKf,
	T dh, T dk, T dl,
	bool bSense,
	T *pTheta, T *pTwoTheta,
	ublas::vector<T>* pVecQ = 0)
{
	const T dDelta = get_plane_dist_tolerance<T>();

	using t_vec = ublas::vector<T>;
	using t_mat = ublas::matrix<T>;

	//try
	{
		t_mat matUB = get_UB(lattice_real, _vec1, _vec2);

		t_vec vechkl = make_vec({dh, dk, dl});
		t_vec vecQ = ublas::prod(matUB, vechkl);
		T dQ = ublas::norm_2(vecQ);

		if(pVecQ) *pVecQ = vecQ;

		if(std::fabs(vecQ[2]) > dDelta)
		{
			std::string strErr("Position not in scattering plane.");
			throw Err(strErr);
		}

		*pTwoTheta = get_sample_twotheta(dKi/get_one_angstrom<T>(), dKf/get_one_angstrom<T>(), dQ/get_one_angstrom<T>(), bSense) / get_one_radian<T>();
		T dKiQ = get_angle_ki_Q(dKi/get_one_angstrom<T>(), dKf/get_one_angstrom<T>(), dQ/get_one_angstrom<T>(), /*bSense*/1) / get_one_radian<T>();
		vecQ.resize(2, true);

		// sample rotation = angle between ki and first orientation reflex (plus an arbitrary, but fixed constant)
		T dAngleKiOrient1 = -dKiQ - vec_angle(vecQ);
		*pTheta = dAngleKiOrient1 + get_pi<T>()/T(2);	// a3 convention would be: kiorient1 + pi
		if(!bSense) *pTheta = -*pTheta;
	}
	/*catch(const std::exception& ex)
	{
		log_err(ex.what());
	}*/
}

template<typename T=double>
void get_hkl_from_tas_angles(const Lattice<T>& lattice_real,
	const ublas::vector<T>& _vec1, const ublas::vector<T>& _vec2,
	T dm, T da, T th_m, T th_a, T _th_s, T _tt_s,
	bool bSense_m, bool bSense_a, bool bSense_s,
	T* h, T* k, T* l,
	T* pki=0, T* pkf=0, T* pE=0, T* pQ=0,
	ublas::vector<T>* pVecQ = 0)
{
	using t_vec = ublas::vector<T>;
	using t_mat = ublas::matrix<T>;

	T th_s = _th_s;
	T tt_s = _tt_s;
	if(!bSense_s)
	{
		th_s = -th_s;
		tt_s = -tt_s;
	}

	T ki = get_mono_k(th_m*get_one_radian<T>(), dm*get_one_angstrom<T>(), bSense_m)*get_one_angstrom<T>();
	T kf = get_mono_k(th_a*get_one_radian<T>(), da*get_one_angstrom<T>(), bSense_a)*get_one_angstrom<T>();
	T E = get_energy_transfer(ki/get_one_angstrom<T>(), kf/get_one_angstrom<T>()) / get_one_meV<T>();
	T Q = get_sample_Q(ki/get_one_angstrom<T>(), kf/get_one_angstrom<T>(), tt_s*get_one_radian<T>())*get_one_angstrom<T>();
	T kiQ = get_angle_ki_Q(ki/get_one_angstrom<T>(), kf/get_one_angstrom<T>(), Q/get_one_angstrom<T>(), /*bSense_s*/1) / get_one_radian<T>();

	th_s += get_pi<T>()/T(2);		// theta here
	T Qvec1 = get_pi<T>() - th_s - kiQ;	// a3 convention


	t_mat matUB = get_UB(lattice_real, _vec1, _vec2);
	t_mat matUBinv;
	if(!inverse(matUB, matUBinv))
		throw Err("Cannot invert UB.");

	t_mat rot = rotation_matrix_3d_z(Qvec1);
	t_vec vecQ = ublas::prod(rot, make_vec({Q,0.,0.}));
	t_vec vechkl = ublas::prod(matUBinv, vecQ);

	if(pVecQ) *pVecQ = vecQ;

	if(vechkl.size() != 3)
		throw Err("Cannot determine hkl.");

	*h = vechkl[0];
	*k = vechkl[1];
	*l = vechkl[2];

	if(pki) *pki = ki;
	if(pkf) *pkf = kf;
	if(pE) *pE = E;
	if(pQ) *pQ = Q;
}


template<typename T=double>
std::ostream& operator<<(std::ostream& ostr, const Lattice<T>& lat)
{
	ostr << "a = " << lat.GetA() << ", ";
	ostr << "b = " << lat.GetB() << ", ";
	ostr << "c = " << lat.GetC() << "; ";

	ostr << "alpha = " << tl::r2d(lat.GetAlpha()) << " deg, ";
	ostr << "beta = " << tl::r2d(lat.GetBeta()) << " deg, ";
	ostr << "gamma = " << tl::r2d(lat.GetGamma()) << " deg; ";

	return ostr;
}

}

#endif
