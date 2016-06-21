/*
 * Chrono helpers
 * @author tweber
 * @date aug-2015
 * @license GPLv2 or GPLv3
 */

#ifndef __CHRONO_HELPERS_H__
#define __CHRONO_HELPERS_H__

#include <chrono>


namespace tl {

// like std::chrono::seconds/minutes/hours, but with variable type
template<typename T = long >
using t_dur_secs = std::chrono::duration<T, std::ratio<1, 1>>;
template<typename T = long >
using t_dur_mins = std::chrono::duration<T, std::ratio<60, 1>>;
template<typename T = long >
using t_dur_hours = std::chrono::duration<T, std::ratio<60*60, 1>>;

template<typename T = long >
using t_dur_days = std::chrono::duration<T, std::ratio<60*60*24, 1>>;

template<typename T = long >
using t_dur_weeks = std::chrono::duration<T, std::ratio<60*60*24*7, 1>>;


// duration since epoch
template<typename t_dur = std::chrono::seconds>
t_dur epoch_dur()
{
	namespace ch = std::chrono;
	return ch::duration_cast<t_dur>(ch::system_clock::now().time_since_epoch());
}

// seconds since epoch
template<typename T=double>
T epoch()
{
	return epoch_dur<t_dur_secs<T>>().count();
}

}
#endif
