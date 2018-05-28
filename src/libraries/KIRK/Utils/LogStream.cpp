/*
 * RT
 * LogStream.cpp
 *
 * @author: Hendrik Schwanekamp
 * @mail:   hendrik.schwanekamp@gmx.net
 *
 * Implements the LogStream class, which provides stream style input for captain KIRKS log.(see Log.h)
 *
 *
 */

// includes
//--------------------
#include "LogStream.h"
//--------------------

// namespace
//--------------------
namespace KIRK {
namespace LOG {
//--------------------

// function definitions of the LogStream class
//-------------------------------------------------------------------
LogStream::LogStream(const LogStream &ls) : logger(ls.logger), sFileline(ls.sFileline), lvl(ls.lvl)
{
	mFilepath = ls.mFilepath;
    str(ls.str());
}

LogStream::LogStream(Log &logger, Level lvl, const char * const filepath, std::string line) : logger(logger), sFileline(line), lvl(lvl)
{
	mFilepath = filepath;
}

LogStream::~LogStream()
{
    logger.log(lvl, mFilepath, sFileline, "%", str());
}


}}