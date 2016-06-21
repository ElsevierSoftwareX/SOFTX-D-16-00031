/*
 * simple numerical integration
 * @author tweber
 * @date june-2015
 * @license GPLv2 or GPLv3
 */

#ifndef __NUMINT_H__
#define __NUMINT_H__

#include <functional>
#include <vector>

namespace tl {


template<class R=double, class A=double>
R numint_trap(const std::function<R(A)>& fkt,
	A x0, A x1)
{
	return 0.5*(x1-x0) * (fkt(x0) + fkt(x1));
}

template<class R=double, class A=double>
R numint_trapN(const std::function<R(A)>& fkt,
	A x0, A x1, unsigned int N)
{
	const A xstep = (x1-x0)/A(N);

	R xsum = fkt(x0) + fkt(x1);
	for(unsigned int i=1; i<N; ++i)
		xsum += 2.*fkt(x0 + i*xstep);

	xsum *= 0.5*xstep;
	return xsum;
}


template<class R=double, class A=double>
R numint_rect(const std::function<R(A)>& fkt,
	A x0, A x1, unsigned int N)
{
	const A xstep = (x1-x0)/A(N);

	R xsum = R(0);
	for(unsigned int i=0; i<N; ++i)
		xsum += fkt(x0 + i*xstep);

	xsum *= xstep;
	return xsum;
}


template<class R=double, class A=double>
R numint_simp(const std::function<R(A)>& fkt,
	A x0, A x1)
{
	return (fkt(x0) + 4.*fkt(0.5*(x0+x1)) + fkt(x1)) * (x1-x0)/6.;
}

template<class R=double, class A=double>
R numint_simpN(const std::function<R(A)>& fkt,
	A x0, A x1, unsigned int N)
{
	const A xstep = (x1-x0)/A(N);
	R xsum = fkt(x0) + fkt(x1);

	for(unsigned int i=1; i<=N/2; ++i)
	{
		xsum += 2.*fkt(x0 + 2.*i*xstep);
		xsum += 4.*fkt(x0 + (2.*i - 1.)*xstep);
	}
	xsum -= 2.*fkt(x0 + 2.*N/2*xstep);

	xsum *= xstep/3.;
	return xsum;
}


// --------------------------------------------------------------------------------

template<class R=double, class A=double>
R convolute(const std::function<R(A)>& fkt0, const std::function<R(A)>& fkt1, 
	A x, A x0, A x1, unsigned int N)
{
	// convolution of fkt0 and fkt1...
	std::function<R(A,A)> fkt = [&fkt0, &fkt1](A t, A tau) -> R
	{
		return fkt0(tau) * fkt1(t-tau);
	};

	// ... at fixed arg x
	std::function<R(A)> fktbnd = std::bind1st(fkt, x);

	return numint_simpN(fktbnd, x0, x1, N);
}



template<class cont_type = std::vector<double>>
cont_type convolute_discrete(const cont_type& f, const cont_type& g)
{
	const std::size_t M = f.size();
	const std::size_t N = g.size();

	cont_type conv;
	conv.reserve(M+N-1);

	for(std::size_t n=0; n<M+N-1; ++n)
	{
		typename cont_type::value_type val = 0.;

		for(std::size_t m=0; m<M; ++m)
			if(n>=m && n-m<N)
				val += f[m] * g[n-m];

		conv.push_back(val);
	}

	return conv;
}


}
#endif
