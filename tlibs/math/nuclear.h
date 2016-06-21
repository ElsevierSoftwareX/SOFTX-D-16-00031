/**
 * nuclear physics formulas
 * @author Tobias Weber
 * @date 2016
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIBS_NUCLEAR__
#define __TLIBS_NUCLEAR__


#include "math.h"
#include "units.h"
#include "../helper/exception.h"
#include <cmath>


namespace tl {

// see: https://de.wikipedia.org/wiki/Mittleres_logarithmisches_Energiedekrement
template<class Y=double>
Y mean_log_E_loss(Y A)
{
	Y a = Y((1-A)*(1-A)) / Y((1+A)*(1+A));
	return Y(1) + a/(Y(1)-a) * std::log(a);
}

template<class Sys, class Y=double>
Y mean_collisions(Y A, const t_energy<Sys,Y>& E_from, const t_energy<Sys,Y>& E_to)
{
	Y dLogLoss = mean_log_E_loss<Y>(A);
	Y dLogE = std::log(Y(E_from/E_to));

	return dLogE / dLogLoss;
}

}

#endif
