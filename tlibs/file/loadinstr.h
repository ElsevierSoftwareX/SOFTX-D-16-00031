/*
 * Loads instrument-specific data files
 * @author tweber
 * @date feb-2015
 * @license GPLv2 or GPLv3
 */

#ifndef __LOADINSTR_H__
#define __LOADINSTR_H__

#include <unordered_map>
#include <map>
#include <vector>
#include <array>
#include <iostream>
#include "../string/string.h"
#include "loaddat.h"

namespace tl{

// interface for instrument-specific data files
template<class _t_real = double>
class FileInstrBase
{
	public:
		using t_real = _t_real;

		typedef std::unordered_map<std::string, std::string> t_mapParams;
		typedef std::vector<std::string> t_vecColNames;
		typedef std::vector<t_real> t_vecVals;
		typedef std::vector<t_vecVals> t_vecDat;

	protected:
		std::array<t_real, 5> GetScanHKLKiKf(const char* pcH, const char* pcK,
			const char* pcL, const char* pcE, std::size_t i) const;

	public:
		FileInstrBase() = default;
		virtual ~FileInstrBase() = default;

		virtual bool Load(const char* pcFile) = 0;

		virtual std::array<t_real, 3> GetSampleLattice() const = 0;
		virtual std::array<t_real, 3> GetSampleAngles() const = 0;
		virtual std::array<t_real, 2> GetMonoAnaD() const = 0;

		virtual std::array<bool, 3> GetScatterSenses() const = 0;
		virtual std::array<t_real, 3> GetScatterPlane0() const = 0;
		virtual std::array<t_real, 3> GetScatterPlane1() const = 0;

		virtual t_real GetKFix() const = 0;
		virtual bool IsKiFixed() const = 0;

		virtual const t_vecVals& GetCol(const std::string& strName) const = 0;
		virtual t_vecVals& GetCol(const std::string& strName) = 0;

		virtual std::size_t GetScanCount() const = 0;
		virtual std::array<t_real, 5> GetScanHKLKiKf(std::size_t i) const = 0;
		virtual bool MergeWith(const FileInstrBase<t_real>* pDat) = 0;

		virtual std::string GetTitle() const = 0;
		virtual std::string GetUser() const = 0;
		virtual std::string GetLocalContact() const = 0;
		virtual std::string GetScanNumber() const = 0;
		virtual std::string GetSampleName() const = 0;
		virtual std::string GetSpacegroup() const = 0;
		virtual std::string GetTimestamp() const = 0;

		virtual const t_vecDat& GetData() const = 0;
		virtual const t_vecColNames& GetColNames() const = 0;
		virtual const t_mapParams& GetAllParams() const = 0;

		virtual std::vector<std::string> GetScannedVars() const = 0;
		virtual std::string GetCountVar() const = 0;
		virtual std::string GetMonVar() const = 0;

		virtual std::string GetScanCommand() const = 0;

	public:
		virtual bool MatchColumn(const std::string& strRegex,
			std::string& strColName, bool bSortByCounts=0, bool bFilterEmpty=1) const;

		static FileInstrBase<t_real>* LoadInstr(const char* pcFile);
};


// psi & ill files
template<class _t_real = double>
class FilePsi : public FileInstrBase<_t_real>
{
	public:
		using t_real = _t_real;
                using t_mapParams = typename FileInstrBase<t_real>::t_mapParams;
                using t_vecColNames = typename FileInstrBase<t_real>::t_vecColNames;
                using t_vecVals = typename FileInstrBase<t_real>::t_vecVals;
                using t_vecDat = typename FileInstrBase<t_real>::t_vecDat;

		// internal parameters in m_mapParams
		typedef std::map<std::string, t_real> t_mapIParams;

	protected:
		t_mapParams m_mapParams;
		t_mapIParams m_mapParameters, m_mapZeros, m_mapVariables, m_mapPosHkl, m_mapScanSteps;
		t_vecColNames m_vecColNames;
		t_vecDat m_vecData;

	protected:
		void ReadData(std::istream& istr);
		void GetInternalParams(const std::string& strAll, t_mapIParams& mapPara);

	public:
		FilePsi() = default;
		virtual ~FilePsi() = default;

		virtual bool Load(const char* pcFile) override;

		void PrintParams(std::ostream& ostr) const;
		const t_mapParams& GetParams() const { return m_mapParams; }

		const std::string& GetColName(std::size_t iCol) const { return m_vecColNames[iCol]; }
		std::size_t GetColCount() const { return m_vecColNames.size(); }

		const t_vecVals& GetCol(std::size_t iCol) const { return m_vecData[iCol]; }
		virtual const t_vecVals& GetCol(const std::string& strName) const override;
		virtual t_vecVals& GetCol(const std::string& strName) override;

	public:
		virtual std::array<t_real, 3> GetSampleLattice() const override;
		virtual std::array<t_real, 3> GetSampleAngles() const override;
		virtual std::array<t_real, 2> GetMonoAnaD() const override;

		virtual std::array<bool, 3> GetScatterSenses() const override;
		virtual std::array<t_real, 3> GetScatterPlane0() const override;
		virtual std::array<t_real, 3> GetScatterPlane1() const override;

		virtual t_real GetKFix() const override;
		virtual bool IsKiFixed() const override;

		std::array<t_real, 4> GetPosHKLE() const;	// zero pos.
		std::array<t_real, 4> GetDeltaHKLE() const;	// scan steps

		virtual std::size_t GetScanCount() const override;
		virtual std::array<t_real, 5> GetScanHKLKiKf(std::size_t i) const override;
		virtual bool MergeWith(const FileInstrBase<t_real>* pDat) override;

		virtual std::string GetTitle() const override;
		virtual std::string GetUser() const override;
		virtual std::string GetLocalContact() const override;
		virtual std::string GetScanNumber() const override;
		virtual std::string GetSampleName() const override;
		virtual std::string GetSpacegroup() const override;
		virtual std::string GetTimestamp() const override;

		virtual const t_vecDat& GetData() const override { return m_vecData; }
		virtual const t_vecColNames& GetColNames() const override { return m_vecColNames; }
		virtual const t_mapParams& GetAllParams() const override { return m_mapParams; }

		virtual std::vector<std::string> GetScannedVars() const override;
		virtual std::string GetCountVar() const override;
		virtual std::string GetMonVar() const override;

		virtual std::string GetScanCommand() const override;
};


// frm/nicos files
template<class _t_real = double>
class FileFrm : public FileInstrBase<_t_real>
{
	public:
		using t_real = _t_real;
                using t_mapParams = typename FileInstrBase<t_real>::t_mapParams;
                using t_vecColNames = typename FileInstrBase<t_real>::t_vecColNames;
                using t_vecVals = typename FileInstrBase<t_real>::t_vecVals;
                using t_vecDat = typename FileInstrBase<t_real>::t_vecDat;

	protected:
		t_mapParams m_mapParams;
		t_vecColNames m_vecQuantities, m_vecUnits;
		t_vecDat m_vecData;

	public:
		FileFrm() = default;
		virtual ~FileFrm() = default;

	protected:
		void ReadHeader(std::istream& istr);
		void ReadData(std::istream& istr);

	public:
		virtual bool Load(const char* pcFile) override;

		virtual std::array<t_real, 3> GetSampleLattice() const override;
		virtual std::array<t_real, 3> GetSampleAngles() const override;
		virtual std::array<t_real, 2> GetMonoAnaD() const override;

		virtual std::array<bool, 3> GetScatterSenses() const override;
		virtual std::array<t_real, 3> GetScatterPlane0() const override;
		virtual std::array<t_real, 3> GetScatterPlane1() const override;

		virtual t_real GetKFix() const override;
		virtual bool IsKiFixed() const override;

		virtual std::size_t GetScanCount() const override;
		virtual std::array<t_real, 5> GetScanHKLKiKf(std::size_t i) const override;
		virtual bool MergeWith(const FileInstrBase<t_real>* pDat) override;

		virtual const t_vecVals& GetCol(const std::string& strName) const override;
		virtual t_vecVals& GetCol(const std::string& strName) override;

		virtual std::string GetTitle() const override;
		virtual std::string GetUser() const override;
		virtual std::string GetLocalContact() const override;
		virtual std::string GetScanNumber() const override;
		virtual std::string GetSampleName() const override;
		virtual std::string GetSpacegroup() const override;
		virtual std::string GetTimestamp() const override;

		virtual const t_vecDat& GetData() const override { return m_vecData; }
		virtual const t_vecColNames& GetColNames() const override { return m_vecQuantities; }
		virtual const t_mapParams& GetAllParams() const override { return m_mapParams; }

		virtual std::vector<std::string> GetScannedVars() const override;
		virtual std::string GetCountVar() const override;
		virtual std::string GetMonVar() const override;

		virtual std::string GetScanCommand() const override;
};


// macs files
template<class _t_real = double>
class FileMacs : public FileInstrBase<_t_real>
{
	public:
		using t_real = _t_real;
                using t_mapParams = typename FileInstrBase<t_real>::t_mapParams;
                using t_vecColNames = typename FileInstrBase<t_real>::t_vecColNames;
                using t_vecVals = typename FileInstrBase<t_real>::t_vecVals;
                using t_vecDat = typename FileInstrBase<t_real>::t_vecDat;

	protected:
		t_mapParams m_mapParams;
		t_vecColNames m_vecQuantities;
		t_vecDat m_vecData;

	public:
		FileMacs() = default;
		virtual ~FileMacs() = default;

	protected:
		void ReadHeader(std::istream& istr);
		void ReadData(std::istream& istr);

	public:
		virtual bool Load(const char* pcFile) override;

		virtual std::array<t_real, 3> GetSampleLattice() const override;
		virtual std::array<t_real, 3> GetSampleAngles() const override;
		virtual std::array<t_real, 2> GetMonoAnaD() const override;

		virtual std::array<bool, 3> GetScatterSenses() const override;
		virtual std::array<t_real, 3> GetScatterPlane0() const override;
		virtual std::array<t_real, 3> GetScatterPlane1() const override;

		virtual t_real GetKFix() const override;
		virtual bool IsKiFixed() const override;

		virtual std::size_t GetScanCount() const override;
		virtual std::array<t_real, 5> GetScanHKLKiKf(std::size_t i) const override;
		virtual bool MergeWith(const FileInstrBase<t_real>* pDat) override;

		virtual const t_vecVals& GetCol(const std::string& strName) const override;
		virtual t_vecVals& GetCol(const std::string& strName) override;

		virtual std::string GetTitle() const override;
		virtual std::string GetUser() const override;
		virtual std::string GetLocalContact() const override;
		virtual std::string GetScanNumber() const override;
		virtual std::string GetSampleName() const override;
		virtual std::string GetSpacegroup() const override;
		virtual std::string GetTimestamp() const override;

		virtual const t_vecDat& GetData() const override { return m_vecData; }
		virtual const t_vecColNames& GetColNames() const override { return m_vecQuantities; }
		virtual const t_mapParams& GetAllParams() const override { return m_mapParams; }

		virtual std::vector<std::string> GetScannedVars() const override;
		virtual std::string GetCountVar() const override;
		virtual std::string GetMonVar() const override;

		virtual std::string GetScanCommand() const override;
};


// trisp files
template<class _t_real = double>
class FileTrisp : public FileInstrBase<_t_real>
{
	public:
		using t_real = _t_real;
                using t_mapParams = typename FileInstrBase<t_real>::t_mapParams;
                using t_vecColNames = typename FileInstrBase<t_real>::t_vecColNames;
                using t_vecVals = typename FileInstrBase<t_real>::t_vecVals;
                using t_vecDat = typename FileInstrBase<t_real>::t_vecDat;

	protected:
		t_mapParams m_mapParams;
		t_vecColNames m_vecQuantities;
		t_vecDat m_vecData;

	public:
		FileTrisp() = default;
		virtual ~FileTrisp() = default;

	protected:
		void ReadHeader(std::istream& istr);
		void ReadData(std::istream& istr);

	public:
		virtual bool Load(const char* pcFile) override;

		virtual std::array<t_real, 3> GetSampleLattice() const override;
		virtual std::array<t_real, 3> GetSampleAngles() const override;
		virtual std::array<t_real, 2> GetMonoAnaD() const override;

		virtual std::array<bool, 3> GetScatterSenses() const override;
		virtual std::array<t_real, 3> GetScatterPlane0() const override;
		virtual std::array<t_real, 3> GetScatterPlane1() const override;

		virtual t_real GetKFix() const override;
		virtual bool IsKiFixed() const override;

		virtual std::size_t GetScanCount() const override;
		virtual std::array<t_real, 5> GetScanHKLKiKf(std::size_t i) const override;
		virtual bool MergeWith(const FileInstrBase<t_real>* pDat) override;

		virtual const t_vecVals& GetCol(const std::string& strName) const override;
		virtual t_vecVals& GetCol(const std::string& strName) override;

		virtual std::string GetTitle() const override;
		virtual std::string GetUser() const override;
		virtual std::string GetLocalContact() const override;
		virtual std::string GetScanNumber() const override;
		virtual std::string GetSampleName() const override;
		virtual std::string GetSpacegroup() const override;
		virtual std::string GetTimestamp() const override;

		virtual const t_vecDat& GetData() const override { return m_vecData; }
		virtual const t_vecColNames& GetColNames() const override { return m_vecQuantities; }
		virtual const t_mapParams& GetAllParams() const override { return m_mapParams; }

		virtual std::vector<std::string> GetScannedVars() const override;
		virtual std::string GetCountVar() const override;
		virtual std::string GetMonVar() const override;

		virtual std::string GetScanCommand() const override;
};


// raw data files
template<class _t_real = double>
class FileRaw : public FileInstrBase<_t_real>
{
	public:
		using t_real = _t_real;
                using t_mapParams = typename FileInstrBase<t_real>::t_mapParams;
                using t_vecColNames = typename FileInstrBase<t_real>::t_vecColNames;
                using t_vecVals = typename FileInstrBase<t_real>::t_vecVals;
                using t_vecDat = typename FileInstrBase<t_real>::t_vecDat;

	protected:
		DatFile<t_real, char> m_dat;
		t_vecColNames m_vecCols;

	public:
		FileRaw() = default;
		virtual ~FileRaw() = default;

	public:
		virtual bool Load(const char* pcFile) override;

		virtual std::array<t_real, 3> GetSampleLattice() const override;
		virtual std::array<t_real, 3> GetSampleAngles() const override;
		virtual std::array<t_real, 2> GetMonoAnaD() const override;

		virtual std::array<bool, 3> GetScatterSenses() const override;
		virtual std::array<t_real, 3> GetScatterPlane0() const override;
		virtual std::array<t_real, 3> GetScatterPlane1() const override;

		virtual t_real GetKFix() const override;
		virtual bool IsKiFixed() const override;

		virtual std::size_t GetScanCount() const override;
		virtual std::array<t_real, 5> GetScanHKLKiKf(std::size_t i) const override;
		virtual bool MergeWith(const FileInstrBase<t_real>* pDat) override;

		virtual const t_vecVals& GetCol(const std::string& strName) const override;
		virtual t_vecVals& GetCol(const std::string& strName) override;

		virtual std::string GetTitle() const override;
		virtual std::string GetUser() const override;
		virtual std::string GetLocalContact() const override;
		virtual std::string GetScanNumber() const override;
		virtual std::string GetSampleName() const override;
		virtual std::string GetSpacegroup() const override;
		virtual std::string GetTimestamp() const override;

		virtual const t_vecDat& GetData() const override;
		virtual const t_vecColNames& GetColNames() const override;
		virtual const t_mapParams& GetAllParams() const override;

		virtual std::vector<std::string> GetScannedVars() const override;
		virtual std::string GetCountVar() const override;
		virtual std::string GetMonVar() const override;

		virtual std::string GetScanCommand() const override;
};



// for legacy compatibility
using FileInstr = FileInstrBase<>;

}

#ifdef TLIBS_INC_HDR_IMPLS
	#include "loadinstr_impl.h"
#endif

#endif
