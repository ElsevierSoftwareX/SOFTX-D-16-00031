/*
 * basic quaternion helpers
 *
 * @author: tweber
 * @date: 2013, 3-dec-2014
 * @license GPLv2 or GPLv3
 */

 #ifndef __TLIBS_QUAT_H__
 #define __TLIBS_QUAT_H__

#include <boost/math/quaternion.hpp>
#include "linalg.h"

namespace tl {

namespace math = boost::math;

// algo from:
// http://www.j3d.org/matrix_faq/matrfaq_latest.html#Q55
template<class mat_type=ublas::matrix<double>,
	class quat_type=math::quaternion<typename mat_type::value_type>>
quat_type rot3_to_quat(const mat_type& rot)
{
	using T = typename quat_type::value_type;
	const T tr = trace(rot);
	T v[3], w;

	if(tr > T(0))								// scalar component is largest
	{
		w    = T(0.5) * std::sqrt(tr+T(1));
		v[0] = (rot(2,1) - rot(1,2)) / (T(4)*w);
		v[1] = (rot(0,2) - rot(2,0)) / (T(4)*w);
		v[2] = (rot(1,0) - rot(0,1)) / (T(4)*w);
	}
	else
	{
		for(std::size_t iComp=0; iComp<3; ++iComp)			// find largest vector component
		{
			const std::size_t iM = iComp;		// major comp.
			const std::size_t im1 = (iComp+1)%3;	// minor comp. 1
			const std::size_t im2 = (iComp+2)%3;	// minor comp. 2

			if(rot(iM,iM)>=rot(im1,im1) && rot(iM,iM)>=rot(im2,im2))
			{
				v[iM]  = T(0.5) * std::sqrt(T(1) + rot(iM,iM) - rot(im1,im1) - rot(im2,im2));
				v[im1] = (rot(im1, iM) + rot(iM, im1)) / (v[iM]*T(4));
				v[im2] = (rot(iM, im2) + rot(im2, iM)) / (v[iM]*T(4));
				w      = (rot(im2,im1) - rot(im1,im2)) / (v[iM]*T(4));

				break;
			}

			if(iComp>=2) throw Err("rot3_to_quat: Invalid condition.");
		}
	}

	quat_type quatRet(w, v[0],v[1],v[2]);
	T norm_eucl = math::abs(quatRet);
	return quatRet/norm_eucl;
}


template<class mat_type=ublas::matrix<double>,
	class quat_type=math::quaternion<typename mat_type::value_type>>
mat_type quat_to_rot3(const quat_type& quat)
{
	const quat_type cquat = math::conj(quat);
	const quat_type i(0,1,0,0), j(0,0,1,0), k(0,0,0,1);

	const quat_type cols[] = {quat*i*cquat,
										quat*j*cquat,
										quat*k*cquat};

	mat_type mat(3,3);
	for(std::size_t icol=0; icol<3; ++icol)
	{
		mat(0, icol) = cols[icol].R_component_2();
		mat(1, icol) = cols[icol].R_component_3();
		mat(2, icol) = cols[icol].R_component_4();
	}

	return mat;
}


template<class quat_type=math::quaternion<double>,
	typename T=typename quat_type::value_type>
std::vector<T> quat_to_euler(const quat_type& quat)
{
	T q[] = {quat.R_component_1(),
		quat.R_component_2(),
		quat.R_component_3(),
		quat.R_component_4()};

	// formulas from:
	// http://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
	T phi = std::atan2(T(2)*(q[0]*q[1] + q[2]*q[3]), T(1)-T(2)*(q[1]*q[1] + q[2]*q[2]));
	T theta = std::asin(T(2)*(q[0]*q[2] - q[3]*q[1]));
	T psi = std::atan2(T(2)*(q[0]*q[3] + q[1]*q[2]), T(1)-T(2)*(q[2]*q[2] + q[3]*q[3]));

	std::vector<T> vec = { phi, theta, psi };
	return vec;
}

template<typename T=double, class... Args>
std::vector<T> rotation_angle(const ublas::matrix<T, Args...>& rot)
{
	std::vector<T> vecResult;

	if(rot.size1()!=rot.size2())
		return vecResult;
	if(rot.size1()<2)
		return vecResult;

	if(rot.size2()==2)
	{
		// rot = ( c -s )
		//       ( s  c )
		T angle = std::atan2(rot(1,0), rot(0,0));
		vecResult.push_back(angle);
	}
	else if(rot.size2()==3)
	{
		math::quaternion<T> quat = rot3_to_quat(rot);
		vecResult = quat_to_euler(quat);
	}

	return vecResult;
}



// see: https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation#Quaternion-derived_rotation_matrix
template<typename T=double>
T rotation_angle(const math::quaternion<T>& quat)
{
	//return 2.*std::asin(math::abs(math::unreal(quat)));
	return T(2)*std::acos(quat.R_component_1());
}

template<class t_vec=ublas::vector<double>>
t_vec rotation_axis(const math::quaternion<typename t_vec::value_type>& quat)
{
	using T = typename t_vec::value_type;

	t_vec vec(3);
	vec[0] = quat.R_component_2();
	vec[1] = quat.R_component_3();
	vec[2] = quat.R_component_4();

	typename t_vec::value_type angle = rotation_angle(quat);
	vec /= std::sin(T(0.5)*angle);

	return vec;
}

template<class quat_type=math::quaternion<double>,
	class vec_type=ublas::vector<typename quat_type::value_type>,
	typename T = typename quat_type::value_type>
quat_type rotation_quat(const vec_type& vec, const T angle)
{
	const T s = std::sin(T(0.5)*angle);
	const T c = std::cos(T(0.5)*angle);
	const T n = std::sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

	const T x = s * vec[0] / n;
	const T y = s * vec[1] / n;
	const T z = s * vec[2] / n;
	const T r = c;

	return quat_type(r, x,y,z);
}

template<class quat_type=math::quaternion<double>, typename T = typename quat_type::value_type>
quat_type rotation_quat_x(typename quat_type::value_type angle)
{ return quat_type(std::cos(T(0.5)*angle), std::sin(T(0.5)*angle), T(0), T(0)); }
template<class quat_type=math::quaternion<double>, typename T = typename quat_type::value_type>
quat_type rotation_quat_y(typename quat_type::value_type angle)
{ return quat_type(std::cos(T(0.5)*angle), T(0), std::sin(T(0.5)*angle), T(0)); }
template<class quat_type=math::quaternion<double>, typename T = typename quat_type::value_type>
quat_type rotation_quat_z(typename quat_type::value_type angle)
{ return quat_type(std::cos(T(0.5)*angle), T(0), T(0), std::sin(T(0.5)*angle)); }




template<class quat_type=math::quaternion<double>, typename T = typename quat_type::value_type>
quat_type stereo_proj(const quat_type& quat)
{ return (T(1)+quat)/(T(1)-quat); }

template<class quat_type=math::quaternion<double>, typename T = typename quat_type::value_type>
quat_type stereo_proj_inv(const quat_type& quat)
{ return (T(1)-quat)/(T(1)+quat); }


template<class QUAT>
struct vec_angle_unsigned_impl<QUAT, LinalgType::QUATERNION>
{
	typename QUAT::value_type operator()(const QUAT& q1, const QUAT& q2) const
	{
		typedef typename QUAT::value_type REAL;
		REAL dot = q1.R_component_1()*q2.R_component_1() +
			q1.R_component_2()*q2.R_component_2() +
			q1.R_component_3()*q2.R_component_3() +
			q1.R_component_4()*q2.R_component_4();

		dot /= math::abs(q1);
		dot /= math::abs(q2);

		return std::acos(dot);
	}
};

}
#endif
