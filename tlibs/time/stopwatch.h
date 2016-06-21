/**
 * simple stopwatch for debugging
 * @author tweber
 * @date jan-2015
 */

#ifndef __STOPWATCH_H__
#define __STOPWATCH_H__

#include <chrono>
#include <ctime>

namespace tl {

template<class T = double>
class Stopwatch
{
	public:
		typedef std::chrono::system_clock::time_point t_tp_sys;
		typedef std::chrono::steady_clock::time_point t_tp_st;
		typedef std::chrono::duration<T> t_dur;
		typedef std::chrono::system_clock::duration t_dur_sys;

	protected:
		t_tp_sys m_timeStart;
		t_tp_st m_timeStart_st, m_timeStop_st;

		t_dur m_dur;
		t_dur_sys m_dur_sys;

		T m_dDur = T(0);

	public:
		Stopwatch() = default;
		virtual ~Stopwatch() = default;

		void start()
		{
			m_timeStart = std::chrono::system_clock::now();
			m_timeStart_st = std::chrono::steady_clock::now();
		}

		void stop()
		{
                	m_timeStop_st = std::chrono::steady_clock::now();

			m_dur = std::chrono::duration_cast<t_dur>(m_timeStop_st-m_timeStart_st);
			m_dur_sys = std::chrono::duration_cast<t_dur_sys>(m_dur);

			m_dDur = T(t_dur::period::num)/T(t_dur::period::den) * T(m_dur.count());
		}

		T GetDur() const
		{
			return m_dDur;
		}

		static std::string to_str(const t_tp_sys& t)
		{
			std::time_t tStart = std::chrono::system_clock::to_time_t(t);
			std::tm tmStart = *std::localtime(&tStart);

			char cTime[256];
			std::strftime(cTime, sizeof cTime, "%a %Y-%b-%d %H:%M:%S %Z", &tmStart);
			return std::string(cTime);
		}

		std::string GetStartTimeStr() const { return to_str(m_timeStart); }
		std::string GetStopTimeStr() const { return to_str(m_timeStart+m_dur_sys); }
};

}
#endif
