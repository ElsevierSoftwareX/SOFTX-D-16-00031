/*
 * Exception
 * @author tweber
 * @date 04-mar-2014
 * @license GPLv2 or GPLv3
 */

#ifndef __MY_EXCEPT_H__
#define __MY_EXCEPT_H__

#include <exception>
#include <string>

namespace tl {

class Err : public std::exception
{
	protected:
		std::string m_strErr;

	public:
		Err(const std::string& strErr, bool bErrStr=0) noexcept
			: m_strErr((bErrStr? "Exception: " : "") + strErr)
		{}

		Err(const char* pcErr) noexcept : Err(std::string(pcErr))
		{}

		virtual ~Err() noexcept
		{}

		virtual const char* what() const noexcept override
		{
			return m_strErr.c_str();
		}
};

}
#endif
