/**
 * minimalistic expression evaluator
 * @author tweber
 * @date apr-2016
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIBS_EVAL_IMPL_H__
#define __TLIBS_EVAL_IMPL_H__

#define BOOST_SPIRIT_USE_PHOENIX_V3
#include <boost/spirit/include/qi.hpp>
#include <boost/phoenix/phoenix.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <string>
#include <utility>
#include <unordered_map>
#include <type_traits>
#include "string.h"
#include "../log/log.h"
#include "../math/math.h"
#include "../math/units.h"

namespace tl
{
	// real functions with one parameter
	template<class t_str, class t_val,
		typename std::enable_if<std::is_floating_point<t_val>::value>::type* =nullptr>
	t_val call_func1(const t_str& strName, t_val t)
	{
		//std::cout << "calling " << strName << " with arg " << t << std::endl;
		static const std::unordered_map</*t_str*/std::string, t_val(*)(t_val)> s_funcs =
		{
			{ "sin", std::sin }, { "cos", std::cos }, { "tan", std::tan },
			{ "asin", std::asin }, { "acos", std::acos }, { "atan", std::atan },
			{ "sinh", std::sinh }, { "cosh", std::cosh }, { "tanh", std::tanh },
			{ "asinh", std::asinh }, { "acosh", std::acosh }, { "atanh", std::atanh },

			{ "sqrt", std::sqrt }, { "cbrt", std::cbrt },
			{ "exp", std::exp },
			{ "log", std::log }, { "log2", std::log2 }, { "log10", std::log10 },

			{ "erf", std::erf }, { "erfc", std::erfc }, { "erf_inv", boost::math::erf_inv },

			{ "round", std::round }, { "ceil", std::ceil }, { "floor", std::floor },
			{ "abs", std::abs },
		};

		return s_funcs.at(wstr_to_str(strName))(t);
	}

	// real functions with two parameters
	template<class t_str, class t_val,
		typename std::enable_if<std::is_floating_point<t_val>::value>::type* =nullptr>
	t_val call_func2(const t_str& strName, t_val t1, t_val t2)
	{
		static const std::unordered_map</*t_str*/std::string, t_val(*)(t_val, t_val)> s_funcs =
		{
			{ "pow", std::pow }, { "atan2", std::atan2 },
			{ "mod", std::fmod },
		};

		return s_funcs.at(wstr_to_str(strName))(t1, t2);
	}

	// real constants
	template<class t_str, class t_val,
		typename std::enable_if<std::is_floating_point<t_val>::value>::type* =nullptr>
	t_val get_const(const t_str& strName)
	{
		static const tl::t_energy_si<t_val> meV = get_one_meV<t_val>();
		//static const tl::t_time_si<t_val> picosec = tl::get_one_picosecond<t_val>();
		static const tl::t_time_si<t_val> sec = tl::get_one_second<t_val>();
		static const tl::t_temperature_si<t_val> kelvin = tl::get_one_kelvin<t_val>();
		static const tl::t_action_si<t_val> hbar = tl::get_hbar<t_val>();
		static const tl::t_energy_per_temperature_si<t_val> kB = tl::get_kB<t_val>();

		static const std::unordered_map</*t_str*/std::string, t_val> s_consts =
		{
			{ "pi", get_pi<t_val>() },
			{ "hbar",  t_val(hbar/meV/sec) },	// hbar in [meV s]
			{ "kB",  t_val(kB/meV*kelvin) },	// kB in [meV / K]
		};

		return s_consts.at(wstr_to_str(strName));
	}


	// alternative: int functions with one parameter
	template<class t_str, class t_val,
		typename std::enable_if<std::is_integral<t_val>::value>::type* =nullptr>
	t_val call_func1(const t_str& strName, t_val t)
	{
		static const std::unordered_map</*t_str*/std::string, t_val(*)(t_val)> s_funcs =
		{
			{ "abs", std::abs },
		};

		return s_funcs.at(wstr_to_str(strName))(t);
	}

	// alternative: int functions with two parameters
	template<class t_str, class t_val,
		typename std::enable_if<std::is_integral<t_val>::value>::type* =nullptr>
	t_val call_func2(const t_str& strName, t_val t1, t_val t2)
	{
		static const std::unordered_map</*t_str*/std::string, std::function<t_val(t_val, t_val)>> s_funcs =
		{
			{ "pow", [t1, t2](t_val t1, t_val t2) -> t_val { return t_val(std::pow(t1, t2)); } },
			{ "mod", [t1, t2](t_val t1, t_val t2) -> t_val { return t1%t2; } },
		};

		return s_funcs.at(wstr_to_str(strName))(t1, t2);
	}

	// alternative: int constants
	template<class t_str, class t_val,
		typename std::enable_if<std::is_integral<t_val>::value>::type* =nullptr>
	t_val get_const(const t_str& strName)
	{
		static const std::unordered_map</*t_str*/std::string, t_val> s_consts =
		{
		};

		return s_consts.at(wstr_to_str(strName));
	}


	namespace sp = boost::spirit;
	namespace qi = boost::spirit::qi;
	namespace asc = boost::spirit::ascii;
	namespace ph = boost::phoenix;

	template<class t_str, class t_val, class t_skip=asc::space_type>
	class ExprGrammar : public qi::grammar<
		typename t_str::const_iterator, t_val(), t_skip>
	{
		protected:
			using t_ch = typename t_str::value_type;
			using t_iter = typename t_str::const_iterator;
			using t_valparser = typename std::conditional<
				std::is_floating_point<t_val>::value,
				qi::real_parser<t_val>, qi::int_parser<t_val>>::type;

			qi::rule<t_iter, t_val(), t_skip> m_expr, m_term;
			qi::rule<t_iter, t_val(), t_skip> m_val, m_baseval, m_const;
			qi::rule<t_iter, t_val(), t_skip> m_func1, m_func2;
			qi::rule<t_iter, t_val(), t_skip> m_pm, m_pm_opt, m_p, m_m;
			qi::rule<t_iter, t_str(), t_skip> m_ident;

			void SetErrorHandling()
			{
				m_expr.name("expr");
				m_term.name("term");
				m_val.name("val");
				m_baseval.name("baseval");
				m_const.name("const");
				m_func1.name("func_1arg");
				m_func2.name("func_2args");
				m_pm.name("plusminus");
				m_pm_opt.name("plusminus_opt");
				m_p.name("plus");
				m_m.name("minus");
				m_ident.name("ident");

				qi::on_error<qi::fail>(m_expr,
					ph::bind([](t_iter beg, t_iter err, t_iter end, const sp::info& infoErr)
					{
						std::string strBeg(beg, err);
						std::string strEnd(err, end);
						std::string strRem;
						if(strEnd.length())
							strRem = "remaining \"" + strEnd + "\", ";

						log_err("Parsed \"", strBeg, "\", ", strRem,
							"expected token ", infoErr, ".");
					}, qi::labels::_1, qi::labels::_3, qi::labels::_2, qi::labels::_4));
			}

		public:
			ExprGrammar() : ExprGrammar::base_type(m_expr, "expr")
			{
				// + or -
				m_expr = ((m_pm_opt > m_term) [ qi::_val = qi::_1*qi::_2  ]
					> *((m_p|m_m) > m_term) [ qi::_val += qi::_1*qi::_2 ]);

				m_pm_opt = (m_p | m_m) [ qi::_val = qi::_1 ]
					| qi::eps [ qi::_val = t_val(1) ];
				m_p = qi::char_(t_ch('+')) [ qi::_val = t_val(1) ];
				m_m = qi::char_(t_ch('-')) [ qi::_val = t_val(-1) ];

				// * or /
				m_term = (m_val [ qi::_val = qi::_1 ]
					> *((t_ch('*') > m_val) [ qi::_val *= qi::_1 ]
						| (t_ch('/') > m_val) [ qi::_val /= qi::_1 ]))
					| m_val [ qi::_val = qi::_1 ];

				// pow
				m_val = m_baseval [ qi::_val = qi::_1 ]
					> *((t_ch('^') > m_baseval)
					[ qi::_val = ph::bind([](t_val val1, t_val val2) -> t_val
					{ return std::pow(val1, val2); }, qi::_val, qi::_1)]);

				m_baseval = t_valparser() | m_func2 | m_func1 | m_const
					| (t_ch('(') > m_expr > t_ch(')'));

				// lazy evaluation of constants via phoenix bind
				m_const = m_ident [ qi::_val = ph::bind([](const t_str& strName) -> t_val
					{ return get_const<t_str, t_val>(strName); }, qi::_1) ];

				// lazy evaluation of functions via phoenix bind
				m_func2 = (m_ident >> t_ch('(') >> m_expr >> t_ch(',') >> m_expr >> t_ch(')'))
					[ qi::_val = ph::bind([](const t_str& strName, t_val val1, t_val val2) -> t_val
					{ return call_func2<t_str, t_val>(strName, val1, val2); },
					qi::_1, qi::_2, qi::_3) ];
				m_func1 = ((m_ident >> t_ch('(') >> m_expr >> t_ch(')')))
					[ qi::_val = ph::bind([](const t_str& strName, t_val val) -> t_val
					{ return call_func1<t_str, t_val>(strName, val); },
					qi::_1, qi::_2) ];

				m_ident = qi::lexeme[qi::char_("A-Za-z_") > *qi::char_("0-9A-Za-z_")];

				SetErrorHandling();
			}

			~ExprGrammar() {}
	};

	template<class t_str/*=std::string*/, class t_val/*=double*/>
	std::pair<bool, t_val> eval_expr(const t_str& str) noexcept
	{
		if(trimmed(str).length() == 0)
			return std::make_pair(true, t_val(0));

		try
		{
			using t_iter = typename t_str::const_iterator;
			t_iter beg = str.begin(), end = str.end();
			t_val valRes(0);

			ExprGrammar<t_str, t_val> gram;
			bool bOk = qi::phrase_parse(beg, end, gram, asc::space, valRes);
			if(beg != end)
			{
				bOk = 0;
			}
			return std::make_pair(bOk, valRes);
		}
		catch(const std::exception& ex)
		{
			log_err("Parsing failed with error: ", ex.what(), ".");
			return std::make_pair(false, t_val(0));
		}
	}
}

#endif
