/**
 * Simple logger
 * @author tweber
 * @date 12-sep-2014
 * @license GPLv2 or GPLv3
 */

#include "log.h"
#include <iomanip>
#include <cstring>
#include <chrono>
#include <boost/date_time/c_time.hpp>


namespace tl {
std::recursive_mutex Log::s_mtx;

std::string Log::get_timestamp()
{
	namespace ch = std::chrono;
	using ch::system_clock;
	using boost::date_time::c_time;

	system_clock::time_point now = system_clock::now();

	// milliseconds
	ch::milliseconds msecs = ch::duration_cast<ch::milliseconds>(now.time_since_epoch());
	msecs -= ch::duration_cast<ch::milliseconds>(ch::duration_cast<ch::seconds>(msecs));
	std::ostringstream ostrmsecs;
	ostrmsecs << std::setw(3) << std::setfill('0') << msecs.count();

	// time and date
	std::time_t tm = system_clock::to_time_t(now);
	std::tm tmNow;
	c_time::localtime(&tm, &tmNow);
	// /*std::*/localtime_r(&tm, &tmNow);

	std::string::value_type cTime[64];
	std::strftime(cTime, sizeof cTime, "%Y-%b-%d %H:%M:%S.", &tmNow);
	std::strcat(cTime, ostrmsecs.str().c_str());
	return std::string(cTime);
}

std::string Log::get_thread_id()
{
	std::ostringstream ostr;
	ostr << std::hex << std::this_thread::get_id();
	return ostr.str();
}

std::string Log::get_color(LogColor col, bool bBold)
{
	switch(col)
	{
		case LogColor::RED: return bBold ? "\033[1;31m" : "\033[0;31m";
		case LogColor::YELLOW: return bBold ? "\033[1;33m" : "\033[0;33m";
		case LogColor::BLUE: return bBold ? "\033[1;34m" : "\033[0;34m";
		case LogColor::GREEN: return bBold ? "\033[1;32m" : "\033[0;32m";
		case LogColor::PURPLE: return bBold ? "\033[1;35m" : "\033[0;35m";
		case LogColor::CYAN: return bBold ? "\033[1;36m" : "\033[0;36m";
		case LogColor::WHITE: return bBold ? "\033[1;37m" : "\033[0;37m";
		case LogColor::BLACK: return bBold ? "\033[1;30m" : "\033[0;30m";
		case LogColor::NONE: 
		default: return "\033[0m";
	}
}

void Log::begin_log()
{
	s_mtx.lock();

	std::vector<t_pairOstr>& vecOstrsTh = GetThreadOstrs();
	std::vector<t_pairOstr> vecOstrs = arrayunion({m_vecOstrs, vecOstrsTh});

	for(t_pairOstr &pairOstr : vecOstrs)
	{
		std::ostream *pOstr = pairOstr.first;
		bool bCol = pairOstr.second;

		if(!pOstr)
			continue;
		if(bCol)
			(*pOstr) << get_color(m_col, 1);
		if(m_bShowDate)
			(*pOstr) << get_timestamp() << ", ";
		if(m_bShowThread)
		{
			using t_mapKey = typename t_threadmap::value_type::first_type;

			t_mapKey idThread = std::this_thread::get_id();
			typename t_threadmap::const_iterator iterMap = m_threadmap.find(idThread);
			if(iterMap == m_threadmap.end())
			{
				++m_iNumThreads;
				std::ostringstream ostrThread;
				ostrThread << "Thread " << m_iNumThreads;

				iterMap = m_threadmap.insert({idThread, ostrThread.str()}).first;
			}

			(*pOstr) << iterMap->second << ", ";
			//(*pOstr) << get_thread_id() << ", ";
		}
		(*pOstr) << m_strInfo << ": ";
		if(bCol)
			(*pOstr) << get_color(m_col, 0);
	}
}

void Log::end_log()
{
        std::vector<t_pairOstr>& vecOstrsTh = GetThreadOstrs();
        std::vector<t_pairOstr> vecOstrs = arrayunion({m_vecOstrs, vecOstrsTh});

	for(t_pairOstr& pairOstr : vecOstrs)
	{
		std::ostream *pOstr = pairOstr.first;
		bool bCol = pairOstr.second;

		if(!pOstr)
			continue;
		if(bCol)
			(*pOstr) << get_color(LogColor::NONE);
		(*pOstr) << std::endl;
	}
	s_mtx.unlock();
}

void Log::inc_depth()
{
	std::lock_guard<decltype(s_mtx)> _lck(s_mtx);

	if(m_iDepth++ == 0)
		begin_log();
}

void Log::dec_depth()
{
	std::lock_guard<decltype(s_mtx)> _lck(s_mtx);
	if(--m_iDepth <= 0)
	{
		m_iDepth = 0;
		end_log();
	}
}

Log::Log() : m_vecOstrs{{&std::cerr, 1}}
{}

Log::Log(const std::string& strInfo, LogColor col, std::ostream* pOstr)
	: m_vecOstrs{{pOstr ? pOstr : &std::cerr, 1}},
	  m_strInfo(strInfo), m_col(col)
{}

Log::~Log()
{}

std::vector<Log::t_pairOstr>& Log::GetThreadOstrs()
{
	return m_mapOstrsTh[std::this_thread::get_id()];
}

void Log::AddOstr(std::ostream* pOstr, bool bCol, bool bThreadLocal)
{
	std::lock_guard<decltype(s_mtx)> _lck(s_mtx);
	if(bThreadLocal)
	{
		std::vector<t_pairOstr>& vecOstrsTh = GetThreadOstrs();
		vecOstrsTh.push_back({pOstr, bCol});
	}
	else
	{
		m_vecOstrs.push_back({pOstr, bCol});
	}
}

void Log::RemoveOstr(std::ostream* pOstr)
{
	std::lock_guard<decltype(s_mtx)> _lck(s_mtx);
	using t_iter = std::vector<t_pairOstr>::iterator;

	auto fktDel = [pOstr](const t_iter::value_type& pairOstr) -> bool
	{
		return pairOstr.first == pOstr;
	};

	t_iter iterNewEnd = std::remove_if(m_vecOstrs.begin(), m_vecOstrs.end(), fktDel);
	if(iterNewEnd != m_vecOstrs.end())
		m_vecOstrs.resize(iterNewEnd-m_vecOstrs.begin());

	for(t_mapthreadOstrs::value_type& pairTh : m_mapOstrsTh)
	{
		t_iter iterNewEndTh = std::remove_if(pairTh.second.begin(), pairTh.second.end(), fktDel);
		if(iterNewEndTh != pairTh.second.end())
			pairTh.second.resize(iterNewEndTh - pairTh.second.begin());
	}
}


Log log_info("INFO", LogColor::WHITE, &std::cerr),
	log_warn("WARNING", LogColor::YELLOW, &std::cerr), 
	log_err("ERROR", LogColor::RED, &std::cerr),
	log_crit("CRITICAL", LogColor::PURPLE, &std::cerr),
	log_debug("DEBUG", LogColor::CYAN, &std::cerr);
}
