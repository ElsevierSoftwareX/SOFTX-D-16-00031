/*
 * fftw interface
 *
 * @author Tobias Weber <tweber@frm2.tum.de>
 * @date August 2012
 * @license GPLv2 or GPLv3
 */

#include "fftw.h"
#include "math.h"
#include "../log/log.h"

#include <iostream>
#include <cstring>
#include <memory>

#include <fftw3.h>

namespace tl {

FFTw::FFTw(unsigned int iSize) : m_iSize(iSize)
{
	m_pIn = fftw_malloc(m_iSize * sizeof(fftw_complex));
	m_pOut = fftw_malloc(m_iSize * sizeof(fftw_complex));

	m_pPlan = malloc(sizeof(fftw_plan));
	m_pPlan_inv = malloc(sizeof(fftw_plan));

	fftw_plan* pPlan = (fftw_plan*) m_pPlan;
	fftw_plan* pPlan_inv = (fftw_plan*) m_pPlan_inv;

	*pPlan = fftw_plan_dft_1d(iSize,
				(fftw_complex*)m_pIn, (fftw_complex*)m_pOut,
				FFTW_FORWARD, FFTW_MEASURE);

	*pPlan_inv = fftw_plan_dft_1d(iSize,
				(fftw_complex*)m_pIn, (fftw_complex*)m_pOut,
				FFTW_BACKWARD, FFTW_MEASURE);

	if(!*pPlan || !*pPlan_inv)
		log_err("Fourier: Could not create plan.");
}

FFTw::~FFTw()
{
	fftw_destroy_plan(*(fftw_plan*)m_pPlan);
	fftw_destroy_plan(*(fftw_plan*)m_pPlan_inv);

	free(m_pPlan);
	free(m_pPlan_inv);

	fftw_free(m_pIn);
	fftw_free(m_pOut);
}


void FFTw::trafo(const double* pRealIn, const double *pImagIn,
		double *pRealOut, double *pImagOut,
		bool bInv)
{
	fftw_plan* pPlan = bInv ? (fftw_plan*) m_pPlan_inv : (fftw_plan*) m_pPlan;

	fftw_complex* pIn = (fftw_complex*)m_pIn;
	fftw_complex* pOut = (fftw_complex*)m_pOut;

	for(unsigned int i=0; i<m_iSize; ++i)
	{
		pIn[i][0] = pRealIn ? pRealIn[i] : 0.;
		pIn[i][1] = pImagIn ? pImagIn[i] : 0.;

		pOut[i][0] = 0.;
		pOut[i][1] = 0.;
	}

	fftw_execute(*pPlan);

	for(unsigned int i=0; i<m_iSize; ++i)
	{
		if(pRealOut) pRealOut[i] = pOut[i][0];
		if(pImagOut) pImagOut[i] = pOut[i][1];
	}
}

}
