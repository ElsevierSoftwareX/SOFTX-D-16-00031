/*
 * fftw interface
 *
 * @author Tobias Weber <tweber@frm2.tum.de>
 * @date August 2012
 * @license GPLv2 or GPLv3
 */

#ifndef __FOURIER_FFTW__
#define __FOURIER_FFTW__

#include "dft.h"

namespace tl {


class FFTw : public Fourier_base<double>
{
	protected:
		unsigned int m_iSize;

		void *m_pIn, *m_pOut;
		void *m_pPlan, *m_pPlan_inv;


	public:
		FFTw(unsigned int iSize);
		virtual ~FFTw();

		virtual void trafo(const double* pInR, const double *pInI,
			double *pOutR, double *pOutI, bool bInv=0) override;
};

}

#endif
