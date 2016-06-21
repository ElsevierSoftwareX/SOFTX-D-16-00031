/*
 * TcpClient
 * @author tweber
 * @date aug-2014
 * @license GPLv2 or GPLv3
 *
 * @desc based on Boost's example chat client (c) 2003-2015 by C. M. Kohlhoff:
 * @desc http://www.boost.org/doc/libs/1_61_0/doc/html/boost_asio/example/cpp11/chat/chat_client.cpp
 */

#include "tcp.h"
#include "tcp_impl.h"

template class tl::TcpTxtClient_gen<>;
template class tl::TcpTxtServer_gen<>;
