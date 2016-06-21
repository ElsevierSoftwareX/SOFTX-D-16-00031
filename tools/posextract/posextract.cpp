/*
 * Extracts (hkl), E positions from scan files for monteconvo
 * @author tweber
 * @date 25-jul-2015
 * @license GPLv2
 */

#include "tlibs/file/loadinstr.h"
#include "tlibs/log/log.h"
#include "tlibs/math/neutrons.hpp"
#include <memory>
#include <fstream>


static inline void usage(const char* pcProg)
{
	tl::log_err("Usage:\n", "\t", pcProg, " <scan file> <out file>");
}


static inline void extract_pos(const char* pcIn, const char* pcOut)
{
	std::shared_ptr<tl::FileInstr> ptrInstr(tl::FileInstr::LoadInstr(pcIn));
	tl::FileInstr *pInstr = ptrInstr.get();

	if(!pInstr)
	{
		tl::log_err("Cannot load scan file \"", pcIn, "\".");
		return;
	}


	std::ofstream ofstr(pcOut);
	if(!ofstr.is_open())
	{
		tl::log_err("Cannot load output file \"", pcOut, "\".");
		return;
	}

	ofstr.precision(16);
	ofstr << "#\n";
	ofstr << "# num_neutrons: 1000\n";
	ofstr << "# algo: 1\n";
	ofstr << "#\n";
	ofstr << "# fixed_ki: " << pInstr->IsKiFixed() << "\n";
	ofstr << "# kfix: " << pInstr->GetKFix() << "\n";
	ofstr << "#\n";
	ofstr << "# h k l E data follows\n";
	ofstr << "#\n";

	tl::log_info("Extracting ", pInstr->GetScanCount(), " scan positions.");
	for(std::size_t iScan=0; iScan<pInstr->GetScanCount(); ++iScan)
	{
		std::array<double,5> arrPos = pInstr->GetScanHKLKiKf(iScan);
		double dh = arrPos[0];
		double dk = arrPos[1];
		double dl = arrPos[2];
		double dki = arrPos[3];
		double dkf = arrPos[4];
		double dE = (tl::k2E(dki/tl::angstrom) - tl::k2E(dkf/tl::angstrom))/tl::meV;

		ofstr << std::left << std::setw(20) << dh << " "
			<< std::left << std::setw(20) << dk << " " 
			<< std::left << std::setw(20) << dl << " " 
			<< std::left << std::setw(20) << dE << std::endl;
	}

	tl::log_info("Done.");
	ofstr.close();
}


int main(int argc, char** argv)
{
	if(argc < 3)
	{
		usage(argv[0]);
		return -1;
	}

	extract_pos(argv[1], argv[2]);
	return 0;
}
