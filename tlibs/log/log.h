/**
 * Simple logger
 * @author tweber
 * @date 12-sep-2014
 * @license GPLv2 or GPLv3
 */

#ifndef __LOGGER_H__
#define __LOGGER_H__

#include <vector>
#include <unordered_map>
#include <iostream>
#include <sstream>
#include <thread>
#include <mutex>
#include <utility>

#include "../helper/array.h"


namespace tl {

enum class LogColor
{
	NONE,
	RED, BLUE, GREEN,
	YELLOW, PURPLE, CYAN,
	WHITE, BLACK
};

class Log
{
private:
	int m_iDepth = 0;

protected:
	static std::recursive_mutex s_mtx;

	// pair of ostream and colour flag
	using t_pairOstr = std::pair<std::ostream*, bool>;
	std::vector<t_pairOstr> m_vecOstrs;
	using t_mapthreadOstrs = std::unordered_map<std::thread::id, std::vector<t_pairOstr>>;
	t_mapthreadOstrs m_mapOstrsTh;

	std::string m_strInfo;
	LogColor m_col = LogColor::NONE;

	bool m_bEnabled = 1;
	bool m_bShowDate = 1;

	bool m_bShowThread = 0;
	unsigned int m_iNumThreads = 0;

	using t_threadmap = std::unordered_map<std::thread::id, std::string>;
	t_threadmap m_threadmap;

protected:
	static std::string get_timestamp();
	static std::string get_thread_id();
	static std::string get_color(LogColor col, bool bBold=0);

	std::vector<t_pairOstr>& GetThreadOstrs();

	void begin_log();
	void end_log();

	void inc_depth();
	void dec_depth();

public:
	Log();
	Log(const std::string& strInfo, LogColor col, std::ostream* = nullptr);
	virtual ~Log();

	void AddOstr(std::ostream* pOstr, bool bCol=1, bool bThreadLocal=0);
	void RemoveOstr(std::ostream* pOstr);

#if __cplusplus > 201402L	// C++1z
	template<typename ...t_args>
	void operator()(t_args&&... args)
	{
		if(!m_bEnabled) return;
		begin_log();

		std::vector<t_pairOstr>& vecOstrsTh = GetThreadOstrs();
		std::vector<t_pairOstr> vecOstrs = arrayunion({m_vecOstrs, vecOstrsTh});

		for(t_pairOstr& pair : vecOstrs)
		{
			if(pair.first)
				((*pair.first) << ... << std::forward<t_args>(args));
		}
		end_log();
	}
#else
	template<typename t_arg>
	void operator()(t_arg&& arg)
	{
		if(!m_bEnabled) return;

		inc_depth();
		std::vector<t_pairOstr>& vecOstrsTh = GetThreadOstrs();
		std::vector<t_pairOstr> vecOstrs = arrayunion({m_vecOstrs, vecOstrsTh});

		for(t_pairOstr& pair : vecOstrs)
		{
			if(pair.first)
				(*pair.first) << std::forward<t_arg>(arg);
		}
		dec_depth();
	}

	template<typename t_arg, typename... t_args>
	void operator()(t_arg&& arg, t_args&&... args)
	{
		if(!m_bEnabled) return;

		inc_depth();
		(*this)(std::forward<t_arg>(arg));
		(*this)(std::forward<t_args>(args)...);
		dec_depth();
	}
#endif

	void SetEnabled(bool bEnab) { m_bEnabled = bEnab; }
	void SetShowDate(bool bDate) { m_bShowDate = bDate; }
	void SetShowThread(bool bThread) { m_bShowThread = bThread; }
};


extern Log log_info, log_warn, log_err, log_crit, log_debug;

}
#endif
