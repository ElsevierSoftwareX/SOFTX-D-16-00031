/**
 * Debug helpers
 * @author tweber
 * @date dec-2015
 * @license GPLv2 or GPLv3
 */

#ifndef __TL_DEBUG_H__
#define __TL_DEBUG_H__

#include <string>
#include <boost/version.hpp>


// ----------------------------------------------------------------------------

#if BOOST_VERSION >= 105700

#include <boost/type_index.hpp>

namespace tl{

template<typename T>
std::string get_typename(bool bFull=1)
{
	boost::typeindex::type_index idx;

	if(bFull)
		idx = boost::typeindex::type_id_with_cvr<T>();
	else
		idx = boost::typeindex::type_id<T>();

	return idx.pretty_name();
}
}

#else

#include <typeinfo>

namespace tl{

template<typename T>
std::string get_typename(bool bFull=1)
{
	return std::string(typeid(T).name());
}
}

#endif

// ----------------------------------------------------------------------------


namespace tl{


#ifndef NDEBUG
	extern void log_backtrace();
#else
	static inline void log_backtrace() {}
#endif

}
#endif
