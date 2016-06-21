/*
 * TcpClient
 * @author tweber
 * @date aug-2014
 * @license GPLv2 or GPLv3
 *
 * @desc based on Boost's example chat client (c) 2003-2015 by C. M. Kohlhoff:
 * @desc http://www.boost.org/doc/libs/1_61_0/doc/html/boost_asio/example/cpp11/chat/chat_client.cpp
 */

#ifndef __TCP_CLIENT__
#define __TCP_CLIENT__

#include <boost/asio.hpp>
#include <boost/signals2.hpp>
#include <boost/lockfree/queue.hpp>

#include <string>
#include <list>
#include <thread>


namespace tl {

namespace sys = boost::system;
namespace asio = boost::asio;
namespace ip = boost::asio::ip;
namespace sig = boost::signals2;
namespace lf = boost::lockfree;


template<class t_ch=char, class t_str=std::basic_string<t_ch>>
class TcpTxtClient_gen
{
public:
	static bool get_cmd_tokens(const t_str& str, const t_str& strDelim,
		std::vector<t_str>& vecStr, t_str& strRemainder);

protected:
	t_str m_strHost, m_strService;

	asio::io_service *m_pservice = nullptr;
	ip::tcp::socket *m_psock = nullptr;
	std::thread* m_pthread = nullptr;

	t_str m_strCmdDelim = "\n";
	lf::queue<const t_str*, lf::fixed_sized<false>> m_listWriteBuffer;
	t_str m_strReadBuffer;

	typedef sig::signal<void(const t_str&)> t_sigRecv;
	typedef sig::signal<void(const t_str&, const t_str&)> t_sigDisconn;
	typedef sig::signal<void(const t_str&, const t_str&)> t_sigConn;

	t_sigRecv m_sigRecv;
	t_sigDisconn m_sigDisconn;
	t_sigConn m_sigConn;

public:
	TcpTxtClient_gen();
	virtual ~TcpTxtClient_gen();
	void set_delim(const t_str& strDelim) { m_strCmdDelim = strDelim; }

	void add_receiver(const typename t_sigRecv::slot_type& conn);
	void add_disconnect(const typename t_sigDisconn::slot_type& disconn);
	void add_connect(const typename t_sigConn::slot_type& conn);

	bool connect(const t_str& strHost, const t_str& strService);
	virtual void disconnect(bool bAlwaysSendSignal=0);
	bool is_connected();

	void write(const t_str& str);
	void wait();

protected:
	void flush_write();
	void read_loop();
};


template<class t_ch=char, class t_str=std::basic_string<t_ch>>
class TcpTxtServer_gen : public TcpTxtClient_gen<t_ch, t_str>
{
protected:
	ip::tcp::endpoint* m_pendpoint = nullptr;
	ip::tcp::acceptor* m_pacceptor = nullptr;

protected:
	typedef sig::signal<void(unsigned short iPort)> t_sigServerStart;

	t_sigServerStart m_sigServerStart;

public:
	TcpTxtServer_gen();
	virtual ~TcpTxtServer_gen();

	virtual void disconnect(bool bAlwaysSendSignal=0) override;
	bool start_server(unsigned short iPort);

	void add_server_start(const typename t_sigServerStart::slot_type& start);
};


// for legacy compatibility
using TcpClient = TcpTxtClient_gen<>;
using TcpServer = TcpTxtServer_gen<>;

}

#ifdef TLIBS_INC_HDR_IMPLS
	#include "tcp_impl.h"
#endif
#endif
