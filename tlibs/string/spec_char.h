/*
 * Special chars
 * @author tweber
 * @date 09-mar-14
 * @license GPLv2 or GPLv3
 */

#ifndef __SPEC_CHAR_H__
#define __SPEC_CHAR_H__

#include <string>
#include <unordered_map>

namespace tl {

struct SpecChar
{
	std::string strUTF8;
	std::wstring strUTF16;

	SpecChar() {}
	SpecChar(const char* pcUTF8, const wchar_t* pcUTF16)
			: strUTF8(pcUTF8), strUTF16(pcUTF16)
	{}
};

typedef std::unordered_map<std::string, SpecChar> t_mapSpecChars;


extern void init_spec_chars();
extern void deinit_spec_chars();
extern const std::string& get_spec_char_utf8(const std::string& strChar);
extern const std::wstring& get_spec_char_utf16(const std::string& strChar);

extern const t_mapSpecChars& get_spec_chars();

}

#endif
