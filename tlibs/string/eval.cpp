/**
 * minimalistic expression evaluator
 * @author tweber
 * @date apr-2016
 * @license GPLv2 or GPLv3
 */

#include "eval.h"
#include "eval_impl.h"

namespace tl
{
	template std::pair<bool, double> eval_expr<std::string, double>(const std::string& str) noexcept;
	//template std::pair<bool, float> eval_expr<std::string, float>(const std::string& str) noexcept;
	//template std::pair<bool, int> eval_expr<std::string, int>(const std::string& str) noexcept;
}
