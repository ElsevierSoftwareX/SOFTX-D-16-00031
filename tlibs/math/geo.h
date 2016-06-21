/*
 * linalg and geometry helpers
 *
 * @author: tweber
 * @date: 30-apr-2013
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIBS_GEO_H__
#define __TLIBS_GEO_H__

#include "../helper/flags.h"
#include "../helper/exception.h"
#include "linalg.h"
#include "linalg2.h"
#include "../log/log.h"

#include <iostream>
#include <cmath>
#include <boost/algorithm/minmax_element.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/exception.hpp>


namespace tl {

namespace ublas = boost::numeric::ublas;
namespace math = boost::math;


//------------------------------------------------------------------------------

template<typename T> class Line;

template<typename T> class Plane
{
public:
	using t_vec = ublas::vector<T>;
	using t_mat = ublas::matrix<T>;

protected:
	bool m_bValid = 0;
	t_vec m_vecX0;
	t_vec m_vecDir0, m_vecDir1;
	t_vec m_vecNorm;
	T m_d;

public:
	Plane(const t_vec& vec0,
		const t_vec& dir0, const t_vec& dir1)
		: m_vecX0(vec0), m_vecDir0(dir0), m_vecDir1(dir1)
	{
		m_vecNorm = cross_3(dir0, dir1);
		T tLenNorm = ublas::norm_2(m_vecNorm);
		if(float_equal<T>(tLenNorm, 0.) || tLenNorm!=tLenNorm)
		{
			m_bValid = 0;
			return;
		}
		m_vecNorm /= tLenNorm;

		// Hessian form: vecX0*vecNorm - d = 0
		m_d = ublas::inner_prod(m_vecX0, m_vecNorm);
		m_bValid = 1;
	}

	Plane() = default;
	virtual ~Plane() = default;

	const t_vec& GetX0() const { return m_vecX0; }
	const t_vec& GetDir0() const { return m_vecDir0; }
	const t_vec& GetDir1() const { return m_vecDir1; }
	const t_vec& GetNorm() const { return m_vecNorm; }
	const T& GetD() const { return m_d; }

	T GetDist(const t_vec& vecPt) const
	{
		return ublas::inner_prod(vecPt, m_vecNorm) - m_d;
	}

	T GetAngle(const Plane<T>& plane) const
	{
		return std::acos(GetNorm(), plane.GetNorm());
	}

	// "Lotfusspunkt"
	t_vec GetDroppedPerp(const t_vec& vecP, T *pdDist=0) const
	{
		T dist = GetDist(vecP);
		t_vec vecdropped = vecP - dist*m_vecNorm;

		if(pdDist)
		{
			t_vec vecD = vecP - vecdropped;
			*pdDist = std::sqrt(ublas::inner_prod(vecD, vecD));
		}

		return vecdropped;
	}

	bool IsParallel(const Plane<T>& plane, T eps = std::numeric_limits<T>::epsilon()) const
	{
		return vec_is_collinear<t_vec>(GetNorm(), plane.GetNorm(), eps);
	}

	// http://mathworld.wolfram.com/Plane-PlaneIntersection.html
	bool intersect(const Plane<T>& plane2, Line<T>& lineRet) const
	{
		if(IsParallel(plane2))
			return false;

		const Plane<T>& plane1 = *this;

		// direction vector
		t_vec vecDir = cross_3(plane1.GetNorm(), plane2.GetNorm());

		// find common point in the two planes
		t_mat M = row_matrix({plane1.GetNorm(), plane2.GetNorm()});

		t_vec vecD(2);
		vecD[0] = plane1.GetD();
		vecD[1] = plane2.GetD();

		t_vec vec0(3);
		if(!tl::solve_linear(M, vecD, vec0))
			return 0;

		lineRet = Line<T>(vec0, vecDir);
		return 1;
	}

	bool IsValid() const { return m_bValid; }
};


//------------------------------------------------------------------------------


template<typename T> class Line
{
public:
	using t_vec = ublas::vector<T>;
	using t_mat = ublas::matrix<T>;

protected:
	t_vec m_vecX0;
	t_vec m_vecDir;

public:
	Line() {}
	Line(const t_vec& vec0, const t_vec& dir)
		: m_vecX0(vec0), m_vecDir(dir)
	{}

	virtual ~Line() {}

	t_vec operator()(T t) const
	{
		return m_vecX0 + t*m_vecDir;
	}

	const t_vec& GetX0() const { return m_vecX0; }
	const t_vec& GetDir() const { return m_vecDir; }

	T GetDist(const Line<T>& l1) const
	{
		const Line<T>& l0 = *this;

		t_vec vecNorm = cross_3(l0.GetDir(), l1.GetDir());

		T tnum = std::fabs(ublas::inner_prod(l1.GetX0()-l0.GetX0(), vecNorm));
		T tdenom = ublas::norm_2(vecNorm);

		return tnum/tdenom;
	}

	T GetDist(const t_vec& vecPt) const
	{
		T tnum = ublas::norm_2(cross_3(m_vecDir, vecPt-m_vecX0));
		T tdenom = ublas::norm_2(m_vecDir);

		return tnum / tdenom;
	}

	bool IsParallel(const Line<T>& line, T eps = std::numeric_limits<T>::epsilon()) const
	{
		return vec_is_collinear<t_vec>(GetDir(), line.GetDir(), eps);
	}


	// "Lotfusspunkt"
	t_vec GetDroppedPerp(const t_vec& vecP, T *pdDist=0) const
	{
		T t = ublas::inner_prod(vecP-GetX0(), GetDir()) / ublas::inner_prod(GetDir(), GetDir());
		t_vec vecdropped = operator()(t);

		if(pdDist)
		{
			t_vec vecD = vecP - vecdropped;
			*pdDist = std::sqrt(ublas::inner_prod(vecD, vecD));
		}

		return vecdropped;
	}


	bool GetSide(const t_vec& vecP, T *pdDist=0) const
	{
		const std::size_t N = m_vecDir.size();
		if(N != 2)
		{
			log_err("\"Side of line\" only defined for 2d vectors.");
			return false;
		}

		t_vec vecDropped = GetDroppedPerp(vecP, pdDist);


		t_vec vecNorm(2);
		vecNorm[0] = m_vecDir[1];
		vecNorm[1] = -m_vecDir[0];

		T tDot = ublas::inner_prod(vecP-vecDropped, vecNorm);
		return tDot < T(0);
	}


	// http://mathworld.wolfram.com/Line-PlaneIntersection.html
	bool intersect(const Plane<T>& plane, T& t) const
	{
		const std::size_t N = m_vecDir.size();
		if(N != 3)
		{
			log_err("Line-plane intersection only implemented for 3d vectors.");
			return false;
		}

		const t_vec& posl = this->GetX0();
		const t_vec& dirl = this->GetDir();

		const t_vec& xp0 = plane.GetX0();
		const t_vec xp1 = plane.GetX0() + plane.GetDir0();
		const t_vec xp2 = plane.GetX0() + plane.GetDir1();

		t_mat matDenom(N+1,N+1);
		matDenom(0,0) = 1;		matDenom(0,1) = 1;		matDenom(0,2) = 1;		matDenom(0,3) = 0;
		matDenom(1,0) = xp0[0];	matDenom(1,1) = xp1[0];	matDenom(1,2) = xp2[0];	matDenom(1,3) = dirl[0];
		matDenom(2,0) = xp0[1];	matDenom(2,1) = xp1[1];	matDenom(2,2) = xp2[1];	matDenom(2,3) = dirl[1];
		matDenom(3,0) = xp0[2];	matDenom(3,1) = xp1[2];	matDenom(3,2) = xp2[2];	matDenom(3,3) = dirl[2];

		T denom = determinant(matDenom);
		if(tl::float_equal(denom, 0.))
			return false;

		t_mat matNum(N+1,N+1);
		matNum(0,0) = 1;		matNum(0,1) = 1;		matNum(0,2) = 1;		matNum(0,3) = 1;
		matNum(1,0) = xp0[0];	matNum(1,1) = xp1[0];	matNum(1,2) = xp2[0];	matNum(1,3) = posl[0];
		matNum(2,0) = xp0[1];	matNum(2,1) = xp1[1];	matNum(2,2) = xp2[1];	matNum(2,3) = posl[1];
		matNum(3,0) = xp0[2];	matNum(3,1) = xp1[2];	matNum(3,2) = xp2[2];	matNum(3,3) = posl[2];

		T num = determinant(matNum);

		t = -num / denom;
		return true;
	}

	bool intersect(const Line<T>& line, T& t) const
	{
		const t_vec& pos0 =  this->GetX0();
		const t_vec& pos1 =  line.GetX0();

		const t_vec& dir0 =  this->GetDir();
		const t_vec& dir1 =  line.GetDir();

		const std::size_t N = pos0.size();

		// pos0 + t0*dir0 = pos1 + t1*dir1
		// pos0 - pos1 = t1*dir1 - t0*dir0

		const t_vec pos = pos0-pos1;
		t_mat mat = ublas::identity_matrix<T>(N);

		for(std::size_t i=0; i<N; ++i)
		{
			mat(i, 0) = -dir0[i];
			mat(i, 1) = dir1[i];
		}

		t_mat inv;
		if(!tl::inverse(mat, inv))
		{
			return false;
		}

		t_vec params = ublas::prod(inv, pos);
		t = params[0];

		return true;
	}

	bool GetMiddlePerp(Line<T>& linePerp) const
	{
		const std::size_t N = m_vecDir.size();
		if(N != 2)
		{
			log_err("Perpendicular line only implemented for 2d vectors.");
			return false;
		}

		t_vec vecDir(2);
		vecDir[0] = -m_vecDir[1];
		vecDir[1] = m_vecDir[0];

		t_vec vecPos = this->operator()(0.5);

		linePerp = Line<T>(vecPos, vecDir);
		return true;
	}
};

template<typename T>
std::ostream& operator<<(std::ostream& ostr, const Plane<T>& plane)
{
	ostr << plane.GetX0() << " + s*" << plane.GetDir0()
				<< " + t*" << plane.GetDir1();
	return ostr;
}

template<typename T>
std::ostream& operator<<(std::ostream& ostr, const Line<T>& line)
{
	ostr << line.GetX0() << " + t*" << line.GetDir();
	return ostr;
}


//------------------------------------------------------------------------------

template<class T=double>
class Quadric
{
public:
	using t_vec = ublas::vector<T>;
	using t_mat = ublas::matrix<T>;

protected:
	// general: x^T Q x  +  r x  +  s  =  0
	// here: x^T Q x + s  =  0
	t_mat m_Q = ublas::zero_matrix<T>(3,3);
	//t_vec m_r = ublas::zero_vector<T>(3);
	T m_s = 0;

	t_vec m_vecOffs = ublas::zero_vector<T>(3);

public:
	Quadric() {}
	Quadric(std::size_t iDim)
		: m_Q(ublas::zero_matrix<T>(iDim,iDim))/*, m_r(ublas::zero_vector<T>(iDim))*/
	{}
	Quadric(const t_mat& Q) : m_Q(Q) {}
	Quadric(const t_mat& Q, /*const t_vec& r,*/ T s)
			: m_Q(Q), /*m_r(r),*/ m_s(s) {}
	virtual ~Quadric() {}

	void SetDim(std::size_t iDim) { m_Q.resize(iDim, iDim, 1); }

	const Quadric<T>& operator=(const Quadric<T>& quad)
	{
		this->m_Q = quad.m_Q;
		//this->m_r = quad.m_r;
		this->m_s = quad.m_s;
		this->m_vecOffs = quad.m_vecOffs;

		return *this;
	}

	Quadric<T>& operator=(Quadric<T>&& quad)
	{
		this->m_Q = std::move(quad.m_Q);
		//this->m_r = std::move(quad.m_r);
		this->m_s = std::move(quad.m_s);
		this->m_vecOffs = std::move(quad.m_vecOffs);

		return *this;
	}

	Quadric(const Quadric<T>& quad) { *this = quad; }
	Quadric(Quadric<T>&& quad) { *this = quad; }

	void SetOffset(const t_vec& vec) { m_vecOffs = vec; }
	const t_vec& GetOffset() const { return m_vecOffs; }

	const t_mat& GetQ() const { return m_Q; }
	//const t_vec& GetR() const { return m_r; }
	T GetS() const { return m_s; }

	void SetQ(const t_mat& Q) { m_Q = Q; }
	//void SetR(const t_vec& r) { m_r = r; }
	void SetS(T s) { m_s = s; }

	T operator()(const t_vec& _x) const
	{
		t_vec x = _x-m_vecOffs;

		t_vec vecQ = ublas::prod(m_Q, x);
		T dQ = ublas::inner_prod(x, vecQ);
		//T dR = ublas::inner_prod(m_r, x);

		return dQ /*+ dR*/ + m_s;
	}

	// remove column and row iIdx
	void RemoveElems(std::size_t iIdx)
	{
		m_Q = remove_elems(m_Q, iIdx);
		//m_r = remove_elem(m_r, iIdx);
		m_vecOffs = remove_elem(m_vecOffs, iIdx);
	}

	void transform(const t_mat& S)
	{
		m_Q = tl::transform<t_mat>(m_Q, S, 1);
	}

	// Q = O D O^T
	// O: eigenvecs, D: eigenvals
	bool GetPrincipalAxes(t_mat& matEvecs, std::vector<T>& vecEvals,
		Quadric<T>* pquadPrincipal=nullptr) const
	{
		std::vector<t_vec > evecs;
		if(!eigenvec_sym(m_Q, evecs, vecEvals))
		{
			log_err("Cannot determine eigenvectors.");
			return false;
		}

		sort_eigenvecs<T>(evecs, vecEvals, 1,
			[](T d) -> T { return 1./std::sqrt(d); });

		matEvecs = column_matrix(evecs);

		if(pquadPrincipal)
		{
			pquadPrincipal->SetDim(vecEvals.size());
			pquadPrincipal->SetQ(diag_matrix(vecEvals));
		}
		return true;
	}

	// quad: x^T Q x + s = 0; line: x = x0 + t d
	// (x0 + t d)^T Q (x0 + t d) + s = 0
	// (x0 + t d)^T Q x0 + (x0 + t d)^T Q t d + s = 0
	// (x0^T + t d^T) Q x0 + (x0^T + t d^T) Q t d + s = 0
	// x0^T Q x0 + s  +  (d^T Q x0 + x0^T Q d) t  +  d^T Q d t^2 = 0
	std::vector<T> intersect(const Line<T>& line) const
	{
		const t_mat& Q = GetQ();
		const T& s = m_s;
		const t_vec& d = line.GetDir();
		const t_vec x0 = line.GetX0() - m_vecOffs;;

		// solving at^2 + bt + c = 0 for t
		t_vec vecQd = ublas::prod(Q, d);
		T a = ublas::inner_prod(d, vecQd);

		t_vec vecQx0 = ublas::prod(Q, x0);
		T c = ublas::inner_prod(x0, vecQx0) + s;

		T b = ublas::inner_prod(x0, vecQd);
		b += ublas::inner_prod(d, vecQx0);

		//std::cout << "a=" << a << ", b=" << b << ", c=" << c << std::endl;
		return quadratic_solve(a,b,c);
	}
};

template<class T=double>
std::ostream& operator<<(std::ostream& ostr, const Quadric<T>& quad)
{
	ostr << "Q = " << quad.GetQ() << ", ";
	ostr << "s = " << quad.GetS();
	return ostr;
}


template<class T=double>
class QuadSphere : public Quadric<T>
{
protected:

public:
	QuadSphere() {}
	QuadSphere(std::size_t iDim) : Quadric<T>(iDim) {}

	QuadSphere(T r) : Quadric<T>(3)
	{
		this->m_Q(0,0) =
		this->m_Q(1,1) =
		this->m_Q(2,2) = 1./(r*r);

		this->m_s = -1.;
	}

	QuadSphere(std::size_t iDim, T r) : Quadric<T>(iDim)
	{
		for(std::size_t i=0; i<iDim; ++i)
			this->m_Q(i,i) = 1./(r*r);

		this->m_s = -1.;
	}

	// only valid in principal axis system
	T GetRadius() const { return 1./std::sqrt(this->m_Q(0,0)); }
	T GetVolume() const { return get_ellipsoid_volume(this->m_Q); }

	virtual ~QuadSphere() {}
};


template<class T=double>
class QuadEllipsoid : public Quadric<T>
{
protected:

public:
	QuadEllipsoid() {}
	QuadEllipsoid(std::size_t iDim) : Quadric<T>(iDim) {}

	QuadEllipsoid(T a, T b) : Quadric<T>(2)
	{
		this->m_Q(0,0) = 1./(a*a);
		this->m_Q(1,1) = 1./(b*b);

		this->m_s = -1.;
	}

	QuadEllipsoid(T a, T b, T c) : Quadric<T>(3)
	{
		this->m_Q(0,0) = 1./(a*a);
		this->m_Q(1,1) = 1./(b*b);
		this->m_Q(2,2) = 1./(c*c);

		this->m_s = -1.;
	}

	QuadEllipsoid(T a, T b, T c, T d) : Quadric<T>(4)
	{
		this->m_Q(0,0) = 1./(a*a);
		this->m_Q(1,1) = 1./(b*b);
		this->m_Q(2,2) = 1./(c*c);
		this->m_Q(3,3) = 1./(d*d);

		this->m_s = -1.;
	}

	virtual ~QuadEllipsoid() {}

	// only valid in principal axis system
	T GetRadius(std::size_t i) const { return 1./std::sqrt(std::abs(this->m_Q(i,i))); }
	T GetVolume() const { return get_ellipsoid_volume(this->m_Q); }
};


//------------------------------------------------------------------------------


template<typename T=double>
std::vector<std::size_t> find_zeroes(std::size_t N, const T* pIn)
{
	using t_vec = ublas::vector<T>;
	std::vector<std::size_t> vecIndices;

	for(std::size_t i=0; i<N-1; ++i)
	{
		t_vec zero(2);
		zero[0] = zero[1] = 0.;
		t_vec xdir(2);
		xdir[0] = 1.; xdir[1] = 0.;
		Line<T> xaxis(zero, xdir);

		t_vec pos0(2);
		pos0[0] = 0.; pos0[1] = pIn[i];
		t_vec pos1(2);
		pos1[0] = 1.; pos1[1] = pIn[i+1];
		Line<T> line(pos0, pos1-pos0);

		T param;
		if(!line.intersect(xaxis, param))
		{
			continue;
		}

		t_vec posInters = line(param);
		if(posInters[0]>=0. && posInters[0]<=1.)
			vecIndices.push_back(i);
	}

	return vecIndices;
}

}

#endif
