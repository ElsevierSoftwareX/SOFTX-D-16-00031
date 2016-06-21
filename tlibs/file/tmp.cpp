/*
 * file helper
 * @author tweber
 * @date 07-mar-2013
 * @license GPLv2 or GPLv3
 */

#include "tmp.h"
#include "file.h"
#include "../math/rand.h"
#include "../string/string.h"

#include <boost/filesystem.hpp>
#include <unistd.h>
#include <fcntl.h>
#include <cstdio>


namespace tl {

TmpFile::TmpFile() : m_strPrefix("tlibs_tmp"), m_iHandle(-1) {}
TmpFile::~TmpFile() { close(); }

// cygwin does not seem to have a ::mkstemp...
int TmpFile::mkstemp(std::string& strFile)
{
	static const std::string strChars = "abcdefghijklmnopqrstuvwxyz1234567890_";
	static const std::size_t iLenChars = strChars.length();

	std::string strRnd;
	strRnd.reserve(s_iRndLen);
	for(unsigned int iRnd=0; iRnd<s_iRndLen; ++iRnd)
		strRnd.push_back(strChars[simple_rand(iLenChars)]);

	if(!find_and_replace(strFile, std::string(s_iRndLen, 'X'), strRnd))
		return -1;

	int iFile = ::open(strFile.c_str(), O_RDWR | O_CREAT | O_EXCL, 0600);
	return iFile;
}

bool TmpFile::open()
{
	namespace fs = boost::filesystem;

	fs::path pathTmpDir = fs::temp_directory_path();
	pathTmpDir /= m_strPrefix + "_" + std::string(s_iRndLen, 'X');

	m_strFile = pathTmpDir.string();
	m_iHandle = mkstemp(m_strFile);
	if(m_iHandle == -1)
		return false;

	return true;
}

void TmpFile::close()
{
	if(m_iHandle != -1)
		::close(m_iHandle);

	if(m_strFile != "")
		remove(m_strFile.c_str());
}

const std::string& TmpFile::GetFileName() const
{
	return m_strFile;
}

void TmpFile::SetPrefix(const char* pcStr)
{
	m_strPrefix = pcStr;
}

}
