/*
 * dft stuff
 * @author tweber
 * @date 2012, jan-2015
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIBS_DFT_H__
#define __TLIBS_DFT_H__

#include <complex>
#include <vector>
#include <memory>
#include <unordered_map>
#include <boost/functional/hash.hpp>

#include "math.h"
#include "linalg.h"

namespace tl
{
//------------------------------------------------------------------------------
// standard dft
// dft formulas from here:
// http://www.fftw.org/fftw3_doc/The-1d-Discrete-Fourier-Transform-_0028DFT_0029.html#The-1d-Discrete-Fourier-Transform-_0028DFT_0029
template<typename T=double>
std::complex<T> dft_coeff(int k,
	const T *pReal, const T *pImag, std::size_t n,
	bool bInv=0)
{
	std::complex<T> imag(0., 1.);

	std::complex<T> f(0.,0.);
	for(std::size_t j=0; j<n; ++j)
	{
		std::complex<T> t(pReal?pReal[j]:T(0), pImag?pImag[j]:T(0));

		T dv = T(-2)*get_pi<T>()*T(j)*T(k)/T(n);
		if(bInv) dv = -dv;
		f += t * (std::cos(dv) + imag*std::sin(dv));
	}

	return f;
}

template<typename T=double>
void dft_direct(const T *pRealIn, const T *pImagIn,
	T *pRealOut, T *pImagOut, std::size_t n,
	bool bInv=0, bool bNorm=0)
{
	for(std::size_t k=0; k<n; ++k)
	{
		std::complex<T> f = dft_coeff<T>(k, pRealIn, pImagIn, n, bInv);
		pRealOut[k] = f.real();
		pImagOut[k] = f.imag();

		if(bNorm && bInv)
		{
			pRealOut[k] /= n;
			pImagOut[k] /= n;
		}
	}
}

template<typename T=double>
std::vector<std::complex<T>> dft_direct(const std::vector<std::complex<T>>& vecIn,
	bool bInv=0, bool bNorm=0)
{
	const std::size_t n = vecIn.size();
	std::vector<std::complex<T>> vecOut;
	vecOut.reserve(n);

	std::unique_ptr<T[]> in_real(new T[n]);
	std::unique_ptr<T[]> in_imag(new T[n]);
	std::unique_ptr<T[]> out_real(new T[n]);
	std::unique_ptr<T[]> out_imag(new T[n]);

	T* pInR = in_real.get();
	T* pInI = in_imag.get();
	T* pOutR = out_real.get();
	T* pOutI = out_imag.get();

	for(std::size_t i=0; i<n; ++i)
	{
		pInR[i] = vecIn[i].real();
		pInI[i] = vecIn[i].imag();
	}

	dft_direct(pInR, pInI, pOutR, pOutI, n, bInv, bNorm);

	for(std::size_t i=0; i<n; ++i)
		vecOut.push_back(std::complex<T>(pOutR[i], pOutI[i]));

	return vecOut;
}


//------------------------------------------------------------------------------

template<typename T=double>
std::vector<std::complex<T>> arrs_to_cvec(const T* pReal, const T* pImag, std::size_t N)
{
	std::vector<std::complex<T>> vec;
	vec.reserve(N);

	for(std::size_t n=0; n<N; ++n)
		vec.push_back(std::complex<T>(pReal?pReal[n]:0., pImag?pImag[n]:0.));

	return vec;
}

template<typename T=double>
void cvec_to_arrs(const std::vector<std::complex<T>>& vec, T* pReal, T* pImag)
{
	for(std::size_t n=0; n<vec.size(); ++n)
	{
		if(pReal) pReal[n] = vec[n].real();
		if(pImag) pImag[n] = vec[n].imag();
	}
}


//------------------------------------------------------------------------------


// tShift=1: shift 1 sample to the right
// tShift=-11: shift 1 sample to the left
template<typename T=double>
std::vector<std::complex<T>> dft_shift(const std::vector<std::complex<T>>& vecIn, T tShift=1.)
{
	const std::size_t N = vecIn.size();
	std::vector<std::complex<T>> vecOut;
	vecOut.reserve(N);

	for(std::size_t i=0; i<N; ++i)
	{
		const T c = std::cos(T(2)*get_pi<T>()*i*tShift / N);
		const T s = std::sin(T(2)*get_pi<T>()*i*tShift / N);
		std::complex<T> ph(c, -s);

		vecOut.push_back(vecIn[i]*ph);
	}

	return vecOut;
}

template<typename T=double>
std::vector<std::complex<T>> dft_double(const std::vector<std::complex<T>>& vecIn)
{
	const std::size_t N = vecIn.size();
	std::vector<std::complex<T>> vecOut;
	vecOut.reserve(N*2);

	for(std::size_t i=0; i<2*N; ++i)
		vecOut.push_back(vecIn[i%N]);

	return vecOut;
}


//------------------------------------------------------------------------------

template<typename T=std::size_t>
T count_bits(T imax)
{
	T inum = 0;
	for(; imax!=0; imax>>=1) ++inum;
	return inum;
}

template<typename T=std::size_t>
T bit_reverse(T imax, T inum)
{
	if(imax<2) return inum;

	T irev = 0;
	T ibitcnt = count_bits(imax)-2;

	for(T i=1; i<imax; i<<=1)
	{
		if(inum & i)
			irev |= (1 << ibitcnt);
		--ibitcnt;
	}

	return irev;
}

template<typename T=std::size_t>
std::vector<T> bit_reverse_indices(T imax)
{
	std::vector<T> vec;
	vec.reserve(imax);

	for(T i=0; i<imax; ++i)
		vec.push_back(bit_reverse(imax,i));

	return vec;
}

template<typename T=double>
std::complex<T> fft_factor(T N, T k, bool bInv=0)
{
	T ph = 1.;
	if(bInv) ph = -1.;

	T c = std::cos(T(2)*get_pi<T>()*k/N * ph);
	T s = std::sin(T(2)*get_pi<T>()*k/N * ph);
	return std::complex<T>(c, -s);
}

template<typename T=double>
std::vector<std::complex<T>> fft_reorder(const std::vector<std::complex<T>>& vecIn)
{
	std::vector<std::size_t> vecIdx = bit_reverse_indices(vecIn.size());

	std::vector<std::complex<T>> vecInRev;
	vecInRev.reserve(vecIn.size());

	for(std::size_t i=0; i<vecIn.size(); ++i)
		vecInRev.push_back(vecIn[vecIdx[i]]);

	return vecInRev;
}

template<class t_vec>
std::pair<t_vec, t_vec> split_vec(const t_vec& vec)
{
	std::size_t N = vec.size();
	std::pair<t_vec, t_vec> pair;
	pair.first.reserve(N/2);
	pair.second.reserve(N/2);

	for(std::size_t i=0; i<N/2; ++i)
	{
		pair.first.push_back(vec[i]);
		pair.second.push_back(vec[N/2 + i]);
	}

	return pair;
}


template<typename T=double>
std::vector<std::complex<T>> fft_merge(const std::vector<std::complex<T>>& vecIn,
	bool bInv=0)
{
	typedef std::vector<std::complex<T>> t_vec;

	const std::size_t N = vecIn.size();
	const std::size_t N2 = N/2;

	if(N==0 || N==1)
		return vecIn;

	std::pair<t_vec, t_vec> pair = split_vec(vecIn);
	t_vec vec1 = fft_merge(pair.first, bInv);
	t_vec vec2 = fft_merge(pair.second, bInv);

	t_vec vecOut;
	vecOut.resize(N);

	for(std::size_t i=0; i<N2; ++i)
	{
		vecOut[i] = vec1[i] + vec2[i]*fft_factor<T>(N, i, bInv);
		vecOut[N2+i] = vec1[i] + vec2[i]*fft_factor<T>(N, N2+i, bInv);
	}

	return vecOut;
}

template<typename T=double>
std::vector<std::complex<T>> fft_direct(const std::vector<std::complex<T>>& vecIn,
	bool bInv=0, bool bNorm=0)
{
	const std::size_t n = vecIn.size();
	std::vector<std::complex<T>> vecOut;
	vecOut.resize(n);

	vecOut = fft_reorder(vecIn);
	vecOut = fft_merge(vecOut, bInv);

	if(bInv && bNorm)
		for(std::complex<T>& c : vecOut)
			c /= n;

	return vecOut;
}


//==============================================================================


template<class T=double>
class Fourier_base
{
	public:
		Fourier_base() = default;
		virtual ~Fourier_base() = default;

		virtual void trafo(const T* pInR, const T* pInI,
			T* pOutR, T* pOutI, bool bInv=0) = 0;
};

// dft with pre-calculated coefficients
template<class T=double>
class DFT : public Fourier_base<T>
{
	protected:
		ublas::matrix<std::complex<T>> m_matCoeff;
		ublas::matrix<std::complex<T>> m_matCoeffInv;

	protected:
		void InitCoeffMatrices(std::size_t n)
		{
			m_matCoeff.resize(n,n,0);

			for(std::size_t i=0; i<n; ++i)
				for(std::size_t j=0; j<n; ++j)
				{
					const T c = std::cos(T(2)*get_pi<T>()*i*j / n);
					const T s = std::sin(T(2)*get_pi<T>()*i*j / n);
					m_matCoeff(i,j) = std::complex<T>(c, -s);
				}

			m_matCoeffInv = ublas::conj(m_matCoeff);
		}

	public:
		DFT(std::size_t n)
		{
			InitCoeffMatrices(n);
		}

		virtual ~DFT() = default;

		ublas::vector<std::complex<T>> trafo(
			const ublas::vector<std::complex<T>>& vec, bool bInv=0)
		{
			if(!bInv)
				return ublas::prod(m_matCoeff, vec);
			else
				return ublas::prod(m_matCoeffInv, vec);
		}

		virtual void trafo(const T* pInR, const T* pInI,
			T* pOutR, T* pOutI,
			bool bInv=0) override
		{
			const std::size_t iSize = m_matCoeff.size1();

			typedef ublas::vector<std::complex<T>> t_vec;
			t_vec vecIn(iSize), vecOut;

			for(std::size_t i=0; i<iSize; ++i)
			{
				T dR = pInR ? pInR[i] : 0.;
				T dI = pInI ? pInI[i] : 0.;
				vecIn[i] = std::complex<T>(dR, dI);
			}

			vecOut = trafo(vecIn, bInv);

			for(std::size_t i=0; i<iSize; ++i)
			{
				if(pOutR) pOutR[i] = vecOut[i].real();
				if(pOutI) pOutI[i] = vecOut[i].imag();
			}
		}
};


// fft with pre-calculated coefficients
template<class T=double>
class FFT : public Fourier_base<T>
{
	protected:
		struct CoeffKey
		{
			std::size_t N;
			std::size_t i;

			CoeffKey(std::size_t iN, std::size_t ii) : N(iN), i(ii)
			{}
		};

		struct CoeffKeyHash
		{
			std::size_t operator()(const CoeffKey& key) const
			{
				std::size_t iHsh = key.N;
				boost::hash_combine<std::size_t>(iHsh, key.i);
				return iHsh;
			}
		};

		struct CoeffKeyEqual
		{
			bool operator()(const CoeffKey& key0, const CoeffKey& key1) const
			{
				return key0.N==key1.N && key0.i==key1.i;
			}
		};


	protected:
		typedef std::unordered_map<CoeffKey, std::complex<T>,
			CoeffKeyHash, CoeffKeyEqual> t_mapCoeff;
		t_mapCoeff m_mapCoeff, m_mapCoeffInv;
		std::size_t m_n;

	protected:
		void InitCoeffs(std::size_t n)
		{
			for(std::size_t i=1; i<=n; i<<=1)
			for(std::size_t j=0; j<i; ++j)
			{
				m_mapCoeff[CoeffKey(i,j)] = fft_factor<T>(i, j, 0);
				m_mapCoeffInv[CoeffKey(i,j)] = fft_factor<T>(i, j, 1);
			}
		}

		std::vector<std::complex<T>> FFTMerge(const std::vector<std::complex<T>>& vecIn,
			bool bInv=0)
		{
			typedef std::vector<std::complex<T>> t_vec;
			/*const*/ t_mapCoeff *pCoeff = bInv ? &m_mapCoeffInv : &m_mapCoeff;

			const std::size_t N = vecIn.size();
			const std::size_t N2 = N/2;

			t_vec vecOut;
			vecOut.reserve(N);

			if(N==0 || N==1)
				return vecIn;
			else if(N==2)	// N=2 specialisation for efficiency
			{
				for(std::size_t i=0; i<2; ++i)
					vecOut.push_back(vecIn[0] + vecIn[1]*pCoeff->operator[](CoeffKey(2,i)));
				return vecOut;
			}
			else if(N==4)	// N=4 specialisation for efficiency
			{
				// 2-point, upper
				std::complex<T> cU[] = { vecIn[0] + vecIn[1]*pCoeff->operator[](CoeffKey(2,0)),
					vecIn[0] + vecIn[1]*pCoeff->operator[](CoeffKey(2,1))};

				// 2-point, lower
				std::complex<T> cL[] = { vecIn[2] + vecIn[3]*pCoeff->operator[](CoeffKey(2,0)),
					vecIn[2] + vecIn[3]*pCoeff->operator[](CoeffKey(2,1))};

				// merge 2-points
				for(std::size_t i=0; i<4; ++i)
					vecOut.push_back(cU[i%2] + cL[i%2]*pCoeff->operator[](CoeffKey(4,i)));

				return vecOut;
			}
			else if(N==8)	// N=8 specialisation for efficiency
			{
				// 2-point, top
				std::complex<T> cT[] = { vecIn[0] + vecIn[1]*pCoeff->operator[](CoeffKey(2,0)),
					vecIn[0] + vecIn[1]*pCoeff->operator[](CoeffKey(2,1)) };

				// 2-point, uppper
				std::complex<T> cU[] = { vecIn[2] + vecIn[3]*pCoeff->operator[](CoeffKey(2,0)),
					vecIn[2] + vecIn[3]*pCoeff->operator[](CoeffKey(2,1)) };

				// 2-point, lower
				std::complex<T> cL[] = { vecIn[4] + vecIn[5]*pCoeff->operator[](CoeffKey(2,0)),
					vecIn[4] + vecIn[5]*pCoeff->operator[](CoeffKey(2,1)) };

				// 2-point, bottom
				std::complex<T> cB[] = { vecIn[6] + vecIn[7]*pCoeff->operator[](CoeffKey(2,0)),
					vecIn[6] + vecIn[7]*pCoeff->operator[](CoeffKey(2,1)) };

				// 4-points: merge 2-points, top+upper
				t_vec c4U; c4U.reserve(4);
				for(std::size_t i=0; i<4; ++i)
					c4U.push_back(cT[i%2] + cU[i%2]*pCoeff->operator[](CoeffKey(4,i)));

				// 4-points: merge 2-points, lower+bottom
				t_vec c4L; c4L.reserve(4);
				for(std::size_t i=0; i<4; ++i)
					c4L.push_back(cL[i%2] + cB[i%2]*pCoeff->operator[](CoeffKey(4,i)));

				// merge 4-points
				for(std::size_t i=0; i<8; ++i)
					vecOut.push_back(c4U[i%4] + c4L[i%4]*pCoeff->operator[](CoeffKey(8,i)));
				return vecOut;
			}

			std::pair<t_vec, t_vec> pair = split_vec(vecIn);
			t_vec vec1 = FFTMerge(pair.first, bInv);
			t_vec vec2 = FFTMerge(pair.second, bInv);

			for(std::size_t i=0; i<N2; ++i)
				vecOut.push_back(vec1[i] + vec2[i] * pCoeff->operator[](CoeffKey(N,i)));
			for(std::size_t i=0; i<N2; ++i)
				vecOut.push_back(vec1[i] + vec2[i] * pCoeff->operator[](CoeffKey(N, N2+i)));

			return vecOut;
		}


	public:
		FFT(std::size_t n) : m_n(n)
		{
			InitCoeffs(n);
		}

		virtual ~FFT() = default;

		std::vector<std::complex<T>> trafo(const std::vector<std::complex<T>>& vec, bool bInv=0)
		{
			return FFTMerge(fft_reorder(vec), bInv);

			//if(bInv && bNorm)
			//for(std::complex<T>& c : vecOut)
			//	c /= n;
		}

		virtual void trafo(const T* pInR, const T* pInI,
			T* pOutR, T* pOutI, bool bInv=0) override
		{
			typedef std::vector<std::complex<T>> t_vec;
			t_vec vecIn(m_n), vecOut;

			for(std::size_t i=0; i<m_n; ++i)
			{
				T dR = pInR ? pInR[i] : 0.;
				T dI = pInI ? pInI[i] : 0.;
				vecIn[i] = std::complex<T>(dR, dI);
			}

			vecOut = trafo(vecIn, bInv);

			for(std::size_t i=0; i<m_n; ++i)
			{
				if(pOutR) pOutR[i] = vecOut[i].real();
				if(pOutI) pOutI[i] = vecOut[i].imag();
			}
		}
};

}

#endif
