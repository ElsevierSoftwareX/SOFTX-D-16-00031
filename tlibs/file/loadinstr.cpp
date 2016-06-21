/*
 * Loads instrument-specific data files
 * @author tweber
 * @date feb-2015
 * @license GPLv2 or GPLv3
 */

#include "loadinstr.h"
#include "loadinstr_impl.h"

namespace tl
{
	template FileInstrBase<double>* FileInstrBase<double>::LoadInstr(const char* pcFile);

	template class FilePsi<double>;
	template class FileFrm<double>;
	template class FileMacs<double>;
	template class FileTrisp<double>;
	template class FileRaw<double>;


/*	template FileInstrBase<float>* FileInstrBase<float>::LoadInstr(const char* pcFile);

	template class FilePsi<float>;
	template class FileFrm<float>;
	template class FileMacs<float>;
	template class FileTrisp<float>;
	template class FileRaw<float>;
*/
}
