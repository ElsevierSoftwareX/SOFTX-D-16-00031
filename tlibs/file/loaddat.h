/*
 * raw data loader
 * @author tweber
 * @date mar-2016
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIBS_LOADDAT__
#define __TLIBS_LOADDAT__

#include <string>
#include <fstream>
#include <vector>
#include <unordered_map>

#include "../log/log.h"
#include "../string/string.h"
#if !defined NO_IOSTR
	#include "../file/comp.h"
#endif


namespace tl {

template<class t_real=double, class t_char=char>
class DatFile
{
	public:
		using t_str = std::basic_string<t_char>;
		using t_col = std::vector<t_real>;
		using t_dat = std::vector<t_col>;
		using t_map = std::unordered_map<t_str, t_str>;

	protected:
		bool m_bOk = 0;

		t_char m_chComm = '#';			// comment
		t_str m_strSeps = {'=', ':'};		// key-value separator
		t_str m_strDatSep = {' ', '\t'};	// data separators

		t_dat m_vecCols;
		std::size_t m_iCurLine = 0;

		t_map m_mapHdr;
		std::vector<t_str> m_vecRawHdr;

	protected:
		void ReadHeaderLine(const t_str& _strLine)
		{
			std::size_t iBeg = _strLine.find_first_not_of(m_chComm);
			if(iBeg == t_str::npos) return;

			t_str strLine(_strLine, iBeg, t_str::npos);
			trim(strLine);
			if(strLine.length() == 0)
				return;

			bool bInsert = 1;
			std::pair<t_str, t_str> pair =
				split_first<t_str>(strLine, t_str({m_strSeps}), 1);
			if(pair.first.length()==0 && pair.second.length()==0)
				bInsert=0;

			if(bInsert)
				m_mapHdr.insert(std::move(pair));
			m_vecRawHdr.emplace_back(std::move(strLine));
		}

		void ReadDataLine(const t_str& strLine)
		{
			std::vector<t_real> vecToks;
			get_tokens<t_real, t_str>(strLine, m_strDatSep, vecToks);
			if(vecToks.size() == 0)
				return;

			if(vecToks.size() > m_vecCols.size())
			{
				if(m_vecCols.size() == 0)	// first line with data
				{
					m_vecCols.resize(vecToks.size());
				}
				else				// add zero columns
				{
					std::size_t iRest = vecToks.size()-m_vecCols.size();
					for(std::size_t iCol=0; iCol<iRest; ++iCol)
					{
						t_col vecNew(m_vecCols[0].size(), t_real(0));
						m_vecCols.emplace_back(std::move(vecNew));
					}

					log_warn("Data file loader: Too many elements in line ", m_iCurLine, ".");
				}
			}
			if(m_vecCols.size() < vecToks.size())
				log_warn("Data file loader: Too few elements in line ", m_iCurLine, ".");

			for(std::size_t iCol=0; iCol<vecToks.size(); ++iCol)
				m_vecCols[iCol].push_back(vecToks[iCol]);
			// fill with 0 if less tokens than columns
			for(std::size_t iCol=vecToks.size(); iCol<m_vecCols.size(); ++iCol)
				m_vecCols[iCol].push_back(t_real(0));
		}

	public:
		void Unload()
		{
			m_bOk = 0;
			m_iCurLine = 0;
			m_vecRawHdr.clear();
			m_mapHdr.clear();
			m_vecCols.clear();
		}

		bool Save(std::basic_ostream<t_char>& ostr)
		{
			if(!m_bOk) return 0;

			for(const typename t_map::value_type& val : m_mapHdr)
				ostr << m_chComm << " " <<
					val.first << " " << m_strSeps[0] << " " << val.second << "\n";

			std::size_t iRows = GetRowCount();
			std::size_t iCols = GetColumnCount();

			for(std::size_t iRow=0; iRow<iRows; ++iRow)
			{
				for(std::size_t iCol=0; iCol<iCols; ++iCol)
					ostr << std::setw(16) << m_vecCols[iCol][iRow] << m_strDatSep[0];
				ostr << "\n";
			}
			return 1;
		}

		bool Save(const t_str& strFile)
		{
			std::ofstream ofstr(strFile);
			if(!ofstr) return 0;
			return Save(ofstr);
		}

		bool Load(std::basic_istream<t_char>& istr)
		{
			Unload();
			while(!istr.eof())
			{
				t_str strLine;
				std::getline(istr, strLine);
				++m_iCurLine;
				trim<t_str>(strLine);
				if(strLine.length() == 0)
					continue;

				if(strLine[0] == m_chComm)
					ReadHeaderLine(strLine);
				else
					ReadDataLine(strLine);
			}

			m_bOk = 1;
			return m_bOk;
		}

		bool Load(const t_str& strFile)
		{
			m_bOk = 0;
			std::basic_ifstream<t_char> ifstr(wstr_to_str(strFile));
			if(!ifstr)
				return false;
#if !defined NO_IOSTR
			std::shared_ptr<std::basic_istream<t_char>> ptrIstr =
				create_autodecomp_istream(ifstr);
			if(!ptrIstr)
				return false;
			std::basic_istream<t_char> *pIstr = ptrIstr.get();
#else
			std::basic_istream<t_char> *pIstr = &ifstr;
#endif
			return Load(*pIstr);
		}

		bool IsOk() const { return m_bOk; }

		void SetCommentChar(t_char ch) { m_chComm = ch; }
		void SetSeparatorChars(const t_str& str) { m_strSeps = str; }
		void SetDataSeparators(const t_str& str) { m_strDatSep = str; };

		const t_dat& GetData() const { return m_vecCols; }
		const t_col& GetColumn(std::size_t iCol) const { return m_vecCols[iCol]; }
		t_col& GetColumn(std::size_t iCol) { return m_vecCols[iCol]; }
		std::size_t GetColumnCount() const { return m_vecCols.size(); }
		std::size_t GetRowCount() const { return m_vecCols[0].size(); }

		const std::vector<t_str>& GetRawHeader() const { return m_vecRawHdr; }
		const t_map& GetHeader() const { return m_mapHdr; }

	public:
		DatFile() = default;
		virtual ~DatFile() = default;

		DatFile(const t_str& strFile)
		{
			Load(strFile);
		}
};

}
#endif
