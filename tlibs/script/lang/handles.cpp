/*
 * Handles
 * @author tweber
 * @date 09-oct-2013
 * @license GPLv2 or GPLv3
 */

#include "handles.h"

HandleFile::HandleFile(FILE *pFile) : m_pFile(pFile)
{}

HandleFile::~HandleFile()
{
	if(m_pFile)
	{
		fclose(m_pFile);
		m_pFile = 0;
	}
}




HandleManager::HandleManager()
{
	m_vecHandles.reserve(64);
}

HandleManager::~HandleManager()
{
	for(Handle*& pHandle : m_vecHandles)
	{
		if(pHandle)
		{
			delete pHandle;
			pHandle = 0;
		}
	}

	m_vecHandles.clear();
}

Handle* HandleManager::GetHandle(t_int iIdx)
{
	if(iIdx >= int(m_vecHandles.size()) || iIdx<0)
		return 0;

	return m_vecHandles[iIdx];
}

t_int HandleManager::AddHandle(Handle* pHandle)
{
	m_vecHandles.push_back(pHandle);
	return m_vecHandles.size()-1;
}

void HandleManager::CloseHandle(t_int iIdx)
{
	if(iIdx >= int(m_vecHandles.size()) || iIdx<0)
		return;

	if(m_vecHandles[iIdx])
	{
		delete m_vecHandles[iIdx];
		m_vecHandles[iIdx] = 0;
	}
}

unsigned int HandleManager::CountHandles(HandleType ty) const
{
	unsigned int uiCnt = 0;
	
	for(const Handle* pHandle : m_vecHandles)
	{
		if(!pHandle) continue;

		if(pHandle->GetType() & ty)
			++uiCnt;
	}
	
	return uiCnt;
}

unsigned int HandleManager::CountAllThreads() const
{
	unsigned int uiCnt = 0;
	
	for(const Handle* pHandle : m_vecHandles)
	{
		if(!pHandle) continue;

		if(pHandle->GetType() & HANDLE_THREAD)
			++uiCnt;
		else if(pHandle->GetType() & HANDLE_TASK)
		{
			if(((HandleTask*)pHandle)->IsThread())
				++uiCnt;
		}
	}
	
	return uiCnt;
}
