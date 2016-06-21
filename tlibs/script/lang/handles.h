/*
 * Handles
 * @author tweber
 * @date 09-oct-2013
 * @license GPLv2 or GPLv3
 */

#ifndef __SCRIPT_HANDLES__
#define __SCRIPT_HANDLES__

#include "types.h"
#include "symbol.h"

#include <vector>
#include <thread>
#include <future>
#include <mutex>
#include <stdio.h>


enum HandleType : unsigned int
{
	HANDLE_FILE = (1<<0),

	HANDLE_THREAD = (1<<1),
	HANDLE_TASK = (1<<2),
	HANDLE_MUTEX = (1<<3),


	HANDLE_ASYNCS = HANDLE_THREAD|HANDLE_TASK,
};




class Handle
{
public:
	virtual ~Handle() {}

	virtual HandleType GetType() const = 0;
};


class HandleFile : public Handle
{
protected:
	FILE *m_pFile;

public:
	HandleFile(FILE *pFile);
	virtual ~HandleFile();

	virtual HandleType GetType() const { return HANDLE_FILE; }
};


template<class HANDLE, HandleType HANDLE_TYPE>
class GenericHandle : public Handle
{
protected:
	HANDLE* m_pHandle;

public:
	GenericHandle(HANDLE* pHandle) : m_pHandle(pHandle)
	{}

	virtual ~GenericHandle()
	{
		if(m_pHandle)
		{
			delete m_pHandle;
			m_pHandle = 0;
		}
	}


	virtual HandleType GetType() const { return HANDLE_TYPE; }
	HANDLE* GetInternalHandle() { return m_pHandle; }
};


using HandleThread = GenericHandle<std::thread, HANDLE_THREAD>;
using HandleMutex = GenericHandle<std::mutex, HANDLE_MUTEX>;

class HandleTask : public GenericHandle<std::future<Symbol*>, HANDLE_TASK>
{
	protected:
		bool m_bIsThread = 0;

	public:
		HandleTask(std::future<Symbol*>* pTask, bool bIsThread=0)
				: GenericHandle(pTask), m_bIsThread(bIsThread)
		{}
		
		bool IsThread() const { return m_bIsThread; }
};



class HandleManager
{
protected:
	std::vector<Handle*> m_vecHandles;

public:
	HandleManager();
	virtual ~HandleManager();

	Handle* GetHandle(t_int iIdx);
	t_int AddHandle(Handle* pHandle);
	void CloseHandle(t_int iIdx);
	
	unsigned int CountHandles(HandleType) const;
	unsigned int CountAllThreads() const;
};


#endif
