/*
 * CL helpers
 * @author tweber
 * @date 20-jan-2015
 * @copyright GPLv2 or GPLv3
 */

#ifndef __CL_WRAP_H__
#define __CL_WRAP_H__

#include <vector>
#include <unordered_map>
#include <algorithm>
#include <iostream>
#include <type_traits>
#include <CL/cl.hpp>

namespace tl {


template<typename t_real=double>
bool get_best_cl_dev(cl::Platform& platRet, cl::Device& devRet, 
	cl::Context* pctxRet=nullptr)
{
	std::vector<cl::Platform> vecPlat;
	cl::Platform::get(&vecPlat);

	if(vecPlat.size() == 0)
		return false;


	struct _Dev
	{
		cl::Platform* pPlat = nullptr;
		cl::Device dev;

		cl_device_type devtype;
	};

	std::vector<_Dev> vecAllDevs;

	for(cl::Platform& plat : vecPlat)
	{
		std::vector<cl::Device> vecDevs;
		plat.getDevices(CL_DEVICE_TYPE_ALL, &vecDevs);

		for(cl::Device& dev : vecDevs)
		{
			_Dev _dev;
			_dev.pPlat = &plat;
			_dev.dev = dev;
			_dev.devtype = dev.getInfo<CL_DEVICE_TYPE>();

			std::string strExtensions = dev.getInfo<CL_DEVICE_EXTENSIONS>();

			// needs double type support?
			if(std::is_same<t_real, double>::value)
			{
				bool bHasDouble = (strExtensions.find("cl_khr_fp64") != std::string::npos);
				if(!bHasDouble) continue;
			}

			vecAllDevs.push_back(_dev);
		}
	}

	if(vecAllDevs.size() == 0)
		return false;

	std::sort(vecAllDevs.begin(), vecAllDevs.end(),
		[](const _Dev& dev1, const _Dev& dev2) -> bool
		{
			int (*get_device_score)(cl_device_type ty) = [](cl_device_type ty) -> int
			{
				int iScore = 0;

				if(ty & CL_DEVICE_TYPE_GPU)
					iScore += 1000;
				if(ty & CL_DEVICE_TYPE_ACCELERATOR)
					iScore += 100;
				if(ty & CL_DEVICE_TYPE_CPU)
					iScore += 10;

				return iScore;
			};

			int iScore1 = get_device_score(dev1.devtype);
			int iScore2 = get_device_score(dev2.devtype);

			return iScore1 > iScore2;
		});

	platRet = *vecAllDevs[0].pPlat;
	devRet = vecAllDevs[0].dev;
	if(pctxRet) *pctxRet = cl::Context({devRet});

	return true;
}



template<typename t_real=double>
bool build_cl_src(cl::Context& ctx, const std::string& strSrc,
	cl::Program& prog, std::unordered_map<std::string, cl::Kernel>& mapKerns)
{
	using t_map = std::remove_reference<decltype(mapKerns)>::type;

	cl::Program::Sources vecSrc;
	vecSrc.emplace_back(cl::Program::Sources::value_type(strSrc.c_str(), strSrc.length()+1));

	prog = cl::Program(ctx, vecSrc);
	std::string strBuildOpts;
	if(std::is_same<t_real, double>::value)
		strBuildOpts = "-DUSE_DOUBLE";

	if(prog.build(/*{dev}, */strBuildOpts.c_str()) != CL_SUCCESS)
	{
		std::vector<cl::Device> vecDev;
		ctx.getInfo(CL_CONTEXT_DEVICES, &vecDev);

		for(cl::Device& dev : vecDev)
		{
			std::cerr << prog.getBuildInfo<CL_PROGRAM_BUILD_LOG>(dev) << std::endl;
		}
		return 0;
	}


	std::vector<cl::Kernel> vecKerns;
	if(prog.createKernels(&vecKerns) != CL_SUCCESS)
		return false;

	for(cl::Kernel& kern : vecKerns)
	{
		std::string strName = kern.getInfo<CL_KERNEL_FUNCTION_NAME>();
		mapKerns.insert(t_map::value_type(strName, std::move(kern)));
	}

	return 1;
}



template<typename t_real=double>
const std::string& get_cl_typedefs()
{
	if(std::is_same<t_real, double>::value)
	{
		static const std::string strTypedefs =
		R"RAWSTR(

		#pragma OPENCL EXTENSION cl_khr_fp64: enable
		typedef double t_real;
		typedef double2 t_real2;
		typedef double4 t_real4;

		)RAWSTR";

		return strTypedefs;
	}
	else
	{
		static const std::string strTypedefs =
		R"RAWSTR(

		typedef float t_real;
		typedef float2 t_real2;
		typedef float4 t_real4;

		)RAWSTR";

		return strTypedefs;
	}
}



// ----------------------------------------------------------------------------
// buffers

template<class T, template<class...> class t_cont = std::vector>
inline bool create_cl_readbuf(cl::Context& ctx, const t_cont<T>& cont, cl::Buffer& buf)
{
	cl_int iErr;
	buf = cl::Buffer(ctx, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
		cont.size() * sizeof(typename t_cont<T>::value_type),
		(void*)cont.data(), &iErr);

	return iErr == CL_SUCCESS;
}

template<class T>
inline bool create_cl_writebuf(cl::Context& ctx, std::size_t iNum, cl::Buffer& buf)
{
	cl_int iErr;
	buf = cl::Buffer(ctx, CL_MEM_WRITE_ONLY,
		iNum * sizeof(T),
		nullptr, &iErr);

	return iErr == CL_SUCCESS;
}


template<class T, template<class...> class t_cont = std::vector>
inline bool create_cl_readimg(cl::Context& ctx, const t_cont<T>& cont, cl::Image2D& img,
	std::size_t iW, std::size_t iH,
	cl::ImageFormat fmt = cl::ImageFormat(CL_RGBA/*CL_RGBx*/, CL_UNSIGNED_INT8))
{
	cl_int iErr;
	img = cl::Image2D(ctx, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
		fmt, iW, iH, 0, (void*)cont.data(), &iErr);

	return iErr == CL_SUCCESS;
}

static inline bool create_cl_writeimg(cl::Context& ctx, cl::Image2D& img,
	std::size_t iW, std::size_t iH,
	cl::ImageFormat fmt = cl::ImageFormat(CL_RGBA, CL_UNSIGNED_INT8))
{
	cl_int iErr;
	img = cl::Image2D(ctx, CL_MEM_WRITE_ONLY,
		fmt, iW, iH, 0, nullptr, &iErr);

	return iErr == CL_SUCCESS;
}


}

#endif
