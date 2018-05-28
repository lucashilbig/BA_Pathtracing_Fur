/*
 * RT
 * LogStream.h
 *
 * @author: Hendrik Schwanekamp
 * @mail:   hendrik.schwanekamp@gmx.net
 *
 * Implements the LogStream class, which provides stream style input for captain KIRKS log.(see Log.h)
 *
 *
 */

#ifndef RT_LOGSTREAM_H
#define RT_LOGSTREAM_H

// includes
//--------------------
#include <sstream>
#include <stdexcept>
#include "Log.h"
//--------------------

// namespace
//--------------------
namespace KIRK {
namespace LOG {
//--------------------

//-------------------------------------------------------------------
/**
 * class LogStream
 *
 * usage:
 * The constructor is usually called from the log class. Than you can log using <<. After the ";" the Logstream is
 * destroyed. It writes its message to the log in its destructor.
 */
class LogStream : public std::ostringstream
{
public:
    /**
     * @brief Copy Constructor
     * @param ls Log Stream to copy from
     */
    LogStream(const LogStream& ls);
    /**
     * @brief Constructor
     * @param logger the logger the completed message is written to
     * @param lvl the log level this message is logged with
     * @param sFilepos file name and line number
     */
    LogStream(Log &logger, Level lvl, const char * const filepath, std::string line);
    /**
     * @brief destructor writes internal string buffer to the log
     */
    ~LogStream();

private:
    Level lvl;
    Log &logger;
	const char * mFilepath;
	std::string sFileline;
};

}}
#endif //RT_LOGSTREAM_H
