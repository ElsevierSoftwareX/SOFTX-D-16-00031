/*
 * file helper
 * @author tweber
 * @date 07-mar-2013
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIB_TMPFILE__
#define __TLIB_TMPFILE__

#include <string>


namespace tl {


class TmpFile
{
private:
	static const std::size_t s_iRndLen = 6;

protected:
	std::string m_strFile;
	std::string m_strPrefix;
	int m_iHandle;

public:
	TmpFile();
	virtual ~TmpFile();

	bool open();
	void close();
	const std::string& GetFileName() const;

	void SetPrefix(const char* pcStr);

	static int mkstemp(std::string& strFile);
};

}

#endif
