/*
 * Compression
 * @author tweber
 * @date 15-aug-2013
 * @license GPLv2 or GPLv3
 */

// TODO: Fix: The (de)compressors don't work with wchar_t!
#ifndef __TLIBS_COMP_H__
#define __TLIBS_COMP_H__

#include <memory>
#include <type_traits>
#include <vector>
#include <fstream>
#include <iostream>

#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/copy.hpp>

#include "../string/string.h"
#include "../log/log.h"


namespace tl {

namespace ios = boost::iostreams;


enum class Compressor
{
	GZ, BZ2, Z, XZ,
	AUTO, INVALID
};


// -----------------------------------------------------------------------------


template<class t_char=unsigned char>
Compressor comp_from_magic(const t_char* pcMagic, std::size_t iLen)
{
	if(iLen<2)
		return Compressor::INVALID;
	if(pcMagic[0]==0x1f && pcMagic[1]==0x8b)
		return Compressor::GZ;
	if(pcMagic[0]==0x78 && (pcMagic[1]==0xda ||
		pcMagic[1]==0x9c || pcMagic[1]==0x01))
		return Compressor::Z;


	if(iLen<3)
		return Compressor::INVALID;
	if(pcMagic[0]=='B' && pcMagic[1]=='Z' && pcMagic[2]=='h')
		return Compressor::BZ2;


	if(iLen<6)
		return Compressor::INVALID;
	if(pcMagic[0]==0xfd && pcMagic[1]=='7' && pcMagic[2]=='z' &&
		pcMagic[3]=='X' && pcMagic[4]=='Z' && pcMagic[5]==0x00)
		return Compressor::XZ;


	return Compressor::INVALID;
}


template<class t_str=std::string>
Compressor comp_from_ext(const t_str& strExt)
{
	if(strExt == "gz") return Compressor::GZ;
	else if(strExt == "bz2") return Compressor::BZ2;
	else if(strExt == "xz") return Compressor::XZ;
	else if(strExt == "z") return Compressor::Z;
	return Compressor::INVALID;
}


// -----------------------------------------------------------------------------



template<class t_char=char>
bool decomp_stream_to_stream(std::basic_istream<t_char>& istr,
	std::basic_ostream<t_char>& ostr, Compressor comp=Compressor::AUTO)
{
	typedef typename std::make_unsigned<t_char>::type t_uchar;

	if(comp == Compressor::AUTO)
	{
		t_uchar pcMagic[] = {0,0,0,0,0,0};
		istr.read((t_char*)pcMagic, 6);
		istr.seekg(0, std::ios::beg);

		comp = comp_from_magic(pcMagic, 6);
	}

	using filtering_streambuf = typename std::conditional<
		std::is_same<t_char, char>::value,
		ios::filtering_streambuf<ios::input>, ios::filtering_wstreambuf<ios::input>>::type;
	filtering_streambuf bufInput;

	using t_alloc = std::allocator<t_char>;
	if(comp == Compressor::GZ)
		bufInput.push(ios::basic_gzip_decompressor<t_alloc>());
	else if(comp == Compressor::BZ2)
		bufInput.push(ios::basic_bzip2_decompressor<t_alloc>());
	else if(comp == Compressor::Z)
		bufInput.push(ios::basic_zlib_decompressor<t_alloc>());
	else if(comp == Compressor::XZ)
	{
		log_err("XZ decompression not yet supported.");
		return false;
	}
	else
	{
		log_err("Unknown decompression selected.");
		return false;
	}

	bufInput.push(istr);
	ios::copy(bufInput, ostr);

	return true;
}

template<class t_char=char>
bool comp_stream_to_stream(std::basic_istream<t_char>& istr,
	std::basic_ostream<t_char>& ostr, Compressor comp=Compressor::GZ)
{
	if(comp == Compressor::AUTO)
		comp = Compressor::GZ;

	using filtering_streambuf = typename std::conditional<
		std::is_same<t_char, char>::value,
		ios::filtering_streambuf<ios::input>, ios::filtering_wstreambuf<ios::input>>::type;
	filtering_streambuf bufInput;

	using t_alloc = std::allocator<t_char>;
	if(comp == Compressor::GZ)
		bufInput.push(ios::basic_gzip_compressor<t_alloc>(/*ios::gzip_params(9)*/));
	else if(comp == Compressor::BZ2)
		bufInput.push(ios::basic_bzip2_compressor<t_alloc>());
	else if(comp == Compressor::Z)
		bufInput.push(ios::basic_zlib_compressor<t_alloc>(/*ios::zlib_params(9)*/));
	else if(comp == Compressor::XZ)
	{
		log_err("XZ compression not yet supported.");
		return false;
	}
	else
	{
		log_err("Unknown compression selected.");
		return false;
	}

	bufInput.push(istr);
	ios::copy(bufInput, ostr);

	return true;
}



// -----------------------------------------------------------------------------



template<class t_char=char>
inline bool __comp_file_to_file(const char* pcFileIn, const char* pcFileOut, Compressor comp, bool bDecomp=0)
{
	std::basic_ifstream<t_char> ifstr(pcFileIn, std::ios_base::binary);
	if(!ifstr.is_open())
	{
		log_err("Cannot open \"", pcFileIn, "\".");
		return false;
	}


	std::basic_ofstream<t_char> ofstr(pcFileOut, std::ios_base::binary);
	if(!ofstr.is_open())
	{
		log_err("Cannot open \"", pcFileOut, "\".");
		return false;
	}


	if(comp == Compressor::AUTO)
	{
		std::string strExt;
		if(bDecomp)
			strExt = get_fileext(std::string(pcFileIn));
		else
			strExt = get_fileext(std::string(pcFileOut));

		if(strExt == "gz")
			comp = Compressor::GZ;
		else if(strExt == "bz2")
			comp = Compressor::BZ2;
	}

	if(bDecomp)
		return decomp_stream_to_stream<t_char>(ifstr, ofstr, comp);
	else
		return comp_stream_to_stream<t_char>(ifstr, ofstr, comp);
}

template<class t_char=char>
bool comp_file_to_file(const char* pcFileIn, const char* pcFileOut, Compressor comp=Compressor::AUTO)
{
	return __comp_file_to_file<t_char>(pcFileIn, pcFileOut, comp, 0);
}

template<class t_char=char>
bool decomp_file_to_file(const char* pcFileIn, const char* pcFileOut, Compressor comp=Compressor::AUTO)
{
	return __comp_file_to_file<t_char>(pcFileIn, pcFileOut, comp, 1);
}


//--------------------------------------------------------------------------------


template<class T> struct _comp_alloc
{
	using value_type = T;
	std::shared_ptr<bool> m_pbKeepMem;

	_comp_alloc()
	{
		m_pbKeepMem = std::make_shared<bool>(0);
	}

	_comp_alloc(const _comp_alloc<T>& alloc)
	{
		m_pbKeepMem = alloc.m_pbKeepMem;
	}

	T* allocate(std::size_t iNum)
	{
		T* pT = new T[iNum];
		return pT;
	}
	void deallocate(T* pT, std::size_t)
	{
		if(*m_pbKeepMem) return;
		if(pT) delete[] pT;
	}
};

template<class T> using _t_comp_arr = std::vector<T, _comp_alloc<T>>;


template<class t_char=char>
inline bool __comp_mem_to_mem(const t_char* pcIn, std::size_t iLenIn,
	t_char*& pcOut, std::size_t& iLenOut, Compressor comp, bool bDecomp=0)
{
	ios::stream<ios::basic_array_source<t_char>> istr(pcIn, iLenIn);
	_comp_alloc<t_char> alloc;
	_t_comp_arr<t_char> vecOut(alloc);
	vecOut.reserve(iLenIn);

	using filtering_ostream = typename std::conditional<
		std::is_same<t_char,char>::value,
		ios::filtering_ostream, ios::filtering_wostream>::type;
	filtering_ostream arrOut(ios::back_inserter(vecOut));

	bool bOk=0;
	if(bDecomp)
		bOk = decomp_stream_to_stream<t_char>(istr, arrOut, comp);
	else
		bOk = comp_stream_to_stream<t_char>(istr, arrOut, comp);

	*alloc.m_pbKeepMem = 1;
	iLenOut = vecOut.size();
	pcOut = vecOut.data();

	return bOk;
}


/*
 * example:
 *	char pc[] = "123456\nABCDEF\n\n";
 *	char *pcTst;
 *	std::size_t iLenOut=0;
 *	comp_mem_to_mem<char>(pc, strlen(pc), pcTst, iLenOut, Compressor::BZ2);
 *	std::ofstream ofstr("tst.txt.bz2");
 *	ofstr.write((char*)pcTst, iLenOut);
 *	ofstr.close();
 *	delete[] pcTst;
 */
template<class t_char=char>
bool comp_mem_to_mem(const t_char* pcIn, std::size_t iLenIn, 
	t_char*& pcOut, std::size_t& iLenOut, Compressor comp=Compressor::GZ)
{
	return __comp_mem_to_mem<t_char>(pcIn, iLenIn, pcOut, iLenOut, comp, 0);
}

template<class t_char=char>
bool decomp_mem_to_mem(const t_char* pcIn, std::size_t iLenIn, 
	t_char*& pcOut, std::size_t& iLenOut, Compressor comp=Compressor::AUTO)
{
	return __comp_mem_to_mem<t_char>(pcIn, iLenIn, pcOut, iLenOut, comp, 1);
}



//--------------------------------------------------------------------------------


template<class t_char=char>
inline bool __comp_mem_to_stream(const void* pvIn, std::size_t iLenIn,
	std::basic_ostream<t_char>& ostr, Compressor comp, bool bDecomp=0)
{
	ios::stream<ios::basic_array_source<t_char>> istr((t_char*)pvIn, iLenIn);

	bool bOk=0;
	if(bDecomp)
		bOk = decomp_stream_to_stream<t_char>(istr, ostr, comp);
	else
		bOk = comp_stream_to_stream<t_char>(istr, ostr, comp);

	return bOk;
}

template<class t_char=char>
bool comp_mem_to_stream(const void* pvIn, std::size_t iLenIn, 
	std::basic_ostream<t_char>& ostr, Compressor comp=Compressor::GZ)
{
	return __comp_mem_to_stream<t_char>(pvIn, iLenIn, ostr, comp, 0);
}

template<class t_char=char>
bool decomp_mem_to_stream(const void* pvIn, std::size_t iLenIn, 
	std::ostream& ostr, Compressor comp=Compressor::AUTO)
{
	return __comp_mem_to_stream<t_char>(pvIn, iLenIn, ostr, comp, 1);
}



//--------------------------------------------------------------------------------


template<class t_char=char>
inline bool __comp_mem_to_mem_fix(const void* pvIn, std::size_t iLenIn,
	void* pvOut, std::size_t iLenOut, Compressor comp, bool bDecomp=0)
{
	t_char *pcIn = (t_char*)pvIn;
	t_char *pcOut = (t_char*)pvOut;
	ios::stream<ios::basic_array_source<t_char>> istr(pcIn, iLenIn);
	ios::stream<ios::basic_array_sink<t_char>> ostr(pcOut, iLenOut);

	bool bOk = 0;
	if(bDecomp)
		bOk = decomp_stream_to_stream<t_char>(istr, ostr, comp);
	else
		bOk = comp_stream_to_stream<t_char>(istr, ostr, comp);

	return bOk;
}

template<class t_char=char>
bool comp_mem_to_mem_fix(const void* pvIn, std::size_t iLenIn, 
	void* pvOut, std::size_t iLenOut, Compressor comp=Compressor::GZ)
{
	return __comp_mem_to_mem_fix<t_char>(pvIn, iLenIn, pvOut, iLenOut, comp, 0);
}

template<class t_char=char>
bool decomp_mem_to_mem_fix(const void* pvIn, std::size_t iLenIn, 
	void* pvOut, std::size_t iLenOut, Compressor comp=Compressor::AUTO)
{
	return __comp_mem_to_mem_fix<t_char>(pvIn, iLenIn, pvOut, iLenOut, comp, 1);
}


//--------------------------------------------------------------------------------


template<class t_char=char>
std::shared_ptr<std::basic_istream<t_char>>
create_autodecomp_istream(std::basic_istream<t_char>& istr)
{
	typedef typename std::make_unsigned<t_char>::type t_uchar;

	t_uchar pcMagic[] = {0,0,0,0,0,0};
	std::streampos pos = istr.tellg();
	istr.read((t_char*)pcMagic, 6);
	istr.seekg(pos, std::ios::beg);
	Compressor comp = comp_from_magic(pcMagic, 6);

	using filtering_istream = typename std::conditional<
		std::is_same<t_char,char>::value,
		ios::filtering_istream, ios::filtering_wistream>::type;
	std::shared_ptr<filtering_istream> ptrIstr
		= std::make_shared<filtering_istream>();

	using t_alloc = std::allocator<t_char>;
	if(comp == Compressor::GZ)
		ptrIstr->push(ios::basic_gzip_decompressor<t_alloc>());
	else if(comp == Compressor::BZ2)
		ptrIstr->push(ios::basic_bzip2_decompressor<t_alloc>());
	else if(comp == Compressor::Z)
		ptrIstr->push(ios::basic_zlib_decompressor<t_alloc>());
	else if(comp == Compressor::XZ)
	{
		log_err("XZ decompression not yet supported.");
		return nullptr;
	}

	ptrIstr->push(istr);
	return ptrIstr;
}

template<class t_char=char>
std::shared_ptr<std::basic_ostream<t_char>>
create_comp_ostream(std::basic_ostream<t_char>& ostr,
	Compressor comp = Compressor::INVALID)
{
	using filtering_ostream = typename std::conditional<
		std::is_same<t_char,char>::value,
		ios::filtering_ostream, ios::filtering_wostream>::type;
	std::shared_ptr<filtering_ostream> ptrOstr
		= std::make_shared<filtering_ostream>();

	using t_alloc = std::allocator<t_char>;
	if(comp == Compressor::GZ)
		ptrOstr->push(ios::basic_gzip_compressor<t_alloc>());
	else if(comp == Compressor::BZ2)
		ptrOstr->push(ios::basic_bzip2_compressor<t_alloc>());
	else if(comp == Compressor::Z)
		ptrOstr->push(ios::basic_zlib_compressor<t_alloc>());
	else if(comp == Compressor::XZ)
		log_err("XZ compression not yet supported.");

	ptrOstr->push(ostr);
	return ptrOstr;
}


}
#endif
