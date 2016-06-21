/**
 * Wrapper to choose between FFT and FFTw
 *
 * @author tweber
 * @date jan-2015
 * @license GPLv2 or GPLv3
 */

#ifndef __FOURIER_CHOOSER__
#define __FOURIER_CHOOSER__

#ifdef USE_FFTW
	#include "fftw.h"
#else
	#include "dft.h"
#endif


namespace tl {

template<class t_real>
class Fourier_gen
{
protected:
	std::size_t m_iSize;

#ifdef USE_FFTW
	FFTw m_dft;
#else
	FFT<t_real> m_dft;
#endif

public:
	Fourier_gen(std::size_t iSize) : m_iSize(iSize), m_dft(m_iSize)
	{}
	virtual ~Fourier_gen() = default;

	void fft(const t_real *pRealIn, const t_real *pImagIn,
		t_real *pRealOut, t_real *pImagOut)
	{ m_dft.trafo(pRealIn, pImagIn, pRealOut, pImagOut, 0); }

	void ifft(const t_real* pRealIn, const t_real *pImagIn,
		t_real *pRealOut, t_real *pImagOut)
	{ m_dft.trafo(pRealIn, pImagIn, pRealOut, pImagOut, 1); }
};

using Fourier = Fourier_gen<double>;

}

#endif
