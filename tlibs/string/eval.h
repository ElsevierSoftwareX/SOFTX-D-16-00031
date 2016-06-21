/**
 * minimalistic expression evaluator
 * @author tweber
 * @date apr-2016
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIBS_EVAL_H__
#define __TLIBS_EVAL_H__

#include <string>
#include <utility>

namespace tl
{
	template<class t_str=std::string, class t_val=double>
	std::pair<bool, t_val> eval_expr(const t_str& str) noexcept;
}

#ifdef TLIBS_INC_HDR_IMPLS
	#include "eval_impl.h"
#endif

#endif
