/*
 * Simple Script
 * @author tweber
 * @date 2013
 * @license GPLv2 or GPLv3
 */

#include "script_helper.h"
#include "log/log.h"
#include <fstream>
#include <iostream>

t_char* load_file(const char* pcFile)
{
	t_ifstream ifstr(pcFile/*, std::ios::binary*/);
	//ifstr.imbue(std::locale(""));

	if(!ifstr.is_open())
	{
		tl::log_err("Cannot open \"", pcFile, "\".");
		return 0;
	}

	ifstr.seekg(0, std::ios::end);
	std::size_t iFileLen = ifstr.tellg();
	ifstr.seekg(0, std::ios::beg);

	t_char* pcInput = new t_char[iFileLen+1];
	ifstr.read(pcInput, iFileLen);
	pcInput[iFileLen] = 0;

	if(ifstr.fail())
		tl::log_err("Failed loading file \"", pcFile, "\".");

	return pcInput;
}
