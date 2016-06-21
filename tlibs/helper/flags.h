/*
 * Compiler- and system-specific stuff
 * @author tweber
 * @date 2013
 * @license GPLv2 or GPLv3
 */

#ifndef __COMPILER_FLAGS_H__
#define __COMPILER_FLAGS_H__

#ifdef __CYGWIN__
        #undef __STRICT_ANSI__
#endif

namespace tl {

// normal popen is not thread-safe on all systems
void *my_popen(const char* pcCmd, const char* pcType="w");
int my_pclose(void*);

}


#ifdef __CYGWIN__
	#include "string.h"

	// missing functions
	namespace std
	{
		template<typename T, class t_str=std::string>
		t_str to_string(T t)
		{
			return tl::var_to_str<T, t_str>(t);
		}

		template<class t_str=std::string>
		double stod(const t_str& str)
		{
			return tl::str_to_var<double, t_str>(str);
		}
	}
#endif

#endif
