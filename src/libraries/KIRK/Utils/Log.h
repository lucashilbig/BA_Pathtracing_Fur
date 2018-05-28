#ifndef __LOG_H__
#define __LOG_H__
//** @file */
#include <stdio.h>
#include <mutex>
#include <fstream>
#include <string>
#include <sstream>
#include <ostream>
#include <chrono>
#include <iostream>

// Check if platform is windows, set LOG_PLATFORM_WINDOWS to 1
#ifdef _WIN32
#define LOG_PLATFORM_WINDOWS 1
#endif

// forward declare LogStream
namespace KIRK {
namespace LOG {
class LogStream;
}}

// Super helpful and easy to use macros
/** LOG_DEBUG Makro */
#define LOG_DEBUG(...) KIRK::LOG::Log::getInstance().log(KIRK::LOG::Level::LVL_DEBUG, KIRK::LOG::getRelativePath(__FILE__), LOG_STR(__LINE__), __VA_ARGS__)
/** LOG_INFO Makro */
#define LOG_INFO(...) KIRK::LOG::Log::getInstance().log(KIRK::LOG::Level::LVL_INFO, KIRK::LOG::getRelativePath(__FILE__), LOG_STR(__LINE__), __VA_ARGS__)
/** LOG_WARN Makro */
#define LOG_WARN(...) KIRK::LOG::Log::getInstance().log(KIRK::LOG::Level::LVL_WARN, KIRK::LOG::getRelativePath(__FILE__), LOG_STR(__LINE__),  __VA_ARGS__)
/** LOG_ERROR Makro */
#define LOG_ERROR(...) KIRK::LOG::Log::getInstance().log(KIRK::LOG::Level::LVL_ERROR, KIRK::LOG::getRelativePath(__FILE__), LOG_STR(__LINE__), __VA_ARGS__)
/** LOG_INPUT Makro */
#define LOG_INPUT(...) KIRK::LOG::Log::getInstance().log(KIRK::LOG::Level::LVL_INPUT, KIRK::LOG::getRelativePath(__FILE__), LOG_STR(__LINE__), __VA_ARGS__)


// stream style log makros
#define LOG_DEBUGs() if(KIRK::LOG::Level::LVL_DEBUG < KIRK::LOG::Log::getInstance().getLogLevel()) ; else KIRK::LOG::LogStream(KIRK::LOG::Log::getInstance(), KIRK::LOG::Level::LVL_DEBUG, KIRK::LOG::getRelativePath(__FILE__), LOG_STR(__LINE__))
#define LOG_INFOs() if(KIRK::LOG::Level::LVL_INFO < KIRK::LOG::Log::getInstance().getLogLevel()) ; else KIRK::LOG::LogStream(KIRK::LOG::Log::getInstance(), KIRK::LOG::Level::LVL_INFO, KIRK::LOG::getRelativePath(__FILE__), LOG_STR(__LINE__))
#define LOG_WARNs() if(KIRK::LOG::Level::LVL_WARN < KIRK::LOG::Log::getInstance().getLogLevel()) ; else KIRK::LOG::LogStream(KIRK::LOG::Log::getInstance(), KIRK::LOG::Level::LVL_WARN, KIRK::LOG::getRelativePath(__FILE__), LOG_STR(__LINE__))
#define LOG_ERRORs() if(KIRK::LOG::Level::LVL_ERROR < KIRK::LOG::Log::getInstance().getLogLevel()) ; else KIRK::LOG::LogStream(KIRK::LOG::Log::getInstance(), KIRK::LOG::Level::LVL_ERROR, KIRK::LOG::getRelativePath(__FILE__), LOG_STR(__LINE__))
#define LOG_INPUTs() if(KIRK::LOG::Level::LVL_INPUT < KIRK::LOG::Log::getInstance().getLogLevel()) ; else KIRK::LOG::LogStream(KIRK::LOG::Log::getInstance(), KIRK::LOG::Level::LVL_INPUT, KIRK::LOG::getRelativePath(__FILE__), LOG_STR(__LINE__))

// make Debug log output not to be calculated in release builds
#ifdef NDEBUG
#undef LOG_DEBUG
#define LOG_DEBUG(...) ;
#undef LOG_DEBUGs
#define LOG_DEBUGs() if(false) ; else KIRK::LOG::LogStream(KIRK::LOG::Log::getInstance(), KIRK::LOG::Level::LVL_DEBUG, KIRK::LOG::getRelativePath(__FILE__), LOG_STR(__LINE__))
#endif

#define LOG_STR1(x) #x
#define LOG_STR(x) LOG_STR1(x)

namespace KIRK {
namespace LOG {

constexpr const char * const strEnd(const char * const str)
{
	return *str ? strEnd(str + 1) : str;
}
constexpr const char * const processPath(const char * const start, const char * const end, std::size_t path_level = 0)
{
	return (start < end && *end != '/' && *end != '\\') ?
		processPath(start, end - 1, path_level) :
		((start < end && path_level <= 1) ?
		 processPath(start, end - 1, path_level + 1) :
		 (end + 1)
		 )
		;
}
/*
* @brief constexpr for compiletime path relativization
* @param path Absolute Path.
*/
constexpr const char * const getRelativePath(const char * const path)
{
	return processPath(path, strEnd(path));
}

/**
* @brief Enum defining hierarchic log levels
*/
enum class Level
{
	LVL_ALL,
	LVL_DEBUG,
	LVL_INFO,
	LVL_WARN,
	LVL_ERROR,
	LVL_INPUT,
	LVL_NONE
};

/**
* C(a)pt(ain Kirk)s Log
*
* References:
* http://en.cppreference.com/w/cpp/language/parameter_pack
* https://bitbucket.org/waterreaction/logger
* http://eu.wiley.com/WileyCDA/WileyTitle/productCd-0470932449.html
*/
class Log
{
public:
	/**
	* @brief Returns a reference to the singleton Log object
	* @return Reference to singleton object.
	*/
	static Log& getInstance();

	/**
	* @brief Variadic Log function. Can be used manually though use of macros (LOG_INFO, ...) is recommended.
	* @param log_level The Level the message will be logged with.
	* @param filepath The source file from where the message is logged. Auto filled in when using Macros. Can be empty.
	* @param line The line from where the message is logged. Auto filled in when using Macros. Can be empty.
	* @param in_formatmessage The actual message to be logged. Use % as placeholder for any type that overloads ostream << operator. Use %% to escape %.
	* @param args Replaces % placeholders in in_formatmessage based on order of appearance.
	*/
	template <typename... Args>
	void log(Level log_level, const char * const filepath, std::string line, std::string in_formatmessage, Args... args);

	/**
	* @brief Enable or Disable logging to file
	* @param log_to_file true = enable, false = disable
	*/
	void setLogToFile(bool log_to_file);
	/**
	* @brief Enable or Disable logging to console
	* @param log_to_console true = enable, false = disable
	*/
	void setLogToConsole(bool log_to_console);
	/**
	* @brief Defines the minimum log level to be logged. e.g. Set to Level::WARN to only log warnings and errors.
	* @param log_level Minimum level to be logged.
	*/
	void setLogLevel(Level log_level);
    /**
     * @brief getter for the current log level
     * @return the current log level
     */
    Level getLogLevel() { return m_log_level;}
	/**
	* @brief Use to define a prefix shown in front of every log entry.
	* @param log_prefix String used as prefix.
	*/
	void setLogPrefix(std::string log_prefix);
	/**
	* @brief Use to define a prefix for the logs filename.
	* @param file_prefix String used as filename prefix.
	* @attention Be careful with symbols that cannot be in a filename.
	*/
	void setFilePrefix(std::string file_prefix);
	/**
	* @brief Enable or Disable colored console output
	* @param use_color true = enable, false = disable
	*/
	void setUseColor(bool use_color);
	/**
	* @brief Enable or Disable printing the source of a log call in console. Has no effect on log file output.
	* @param log_source true = enable, false = disable
	*/
	void setLogSource(bool log_source);

protected:
	static Log* m_instance; //!< Holds the instance of the Log for singleton
	std::string m_filename; //!< The full path with name to the logs output file
	std::ofstream m_outstream; //!< The outstream used for writing to output file

	/**
	* @brief Embedded helper class that makes sure the single instance is deleted on shutdown.
	*/
	friend class Cleanup;
	class Cleanup
	{
	public:
		~Cleanup();
	};

	/**
	* @brief Enum for all possible console colors.
	*/
	enum class LogColor
	{
		BLACK,
		DARK_GREY,
		DARK_RED,
		DARK_GREEN,
		DARK_YELLOW,
		DARK_BLUE,
		DARK_MAGENTA,
		DARK_CYAN,
		GREY,
		RED,
		GREEN,
		YELLOW,
		BLUE,
		MAGENTA,
		CYAN,
		WHITE
	};

private:
	Log();
	virtual ~Log();
	Log(const Log&); //Prevent copying
	Log& operator=(const Log&); //Prevent assignment

	/**
	* @brief Generates a filename containing file prefix, time and and Debug Level.
	*/
	void initFilename();
	/**
	* @brief Opens on ofstream if file logging is enabled and closes ofstream if one is open and file logging is disabled.
	*/
	void initOutstream();
	/**
	* @brief Converts Level Enum to string.
	* @param level
	* @return The enums value name as string.
	*/
	static std::string levelToString(Level level);

	/**
	* @brief Recursive variadic helper function for formatting a string by replacing % with args.
	* @param format String to be formatted
	* @param value Current argument to be inserted.
	* @param args Remaing arguments.
	*/
	template<typename T, typename... Args>
	std::string logHelper(std::string format, T value, Args... args);
	/**
	* @brief If no args are passed or no args remain this will be called.
	* @param format String to be formatte
	* @return Returns format param unchanged.
	*/
	std::string logHelper(std::string format);

#if LOG_PLATFORM_WINDOWS
	/**
	* @brief Converts LogColor into WORD color Attr used for Windows console coloring.
	* @param log_color LogColor enum
	* @return WORD color attribute for windows console coloring.
	*/
	unsigned short logColorToColorAttr(LogColor log_color);
#endif
	/**
	* @brief Converts LogColor into ANSI code color number.
	* @param log_color LogColor enum
	* @return const char* ansi code color number
	*/
	const char* logColorToColorAnsi(LogColor log_color);
	/**
	* @brief Converts Level to a predefined LogColor
	* @param level
	* @return LogColor based on level
	*/
	static LogColor levelToLogColor(Level level);
	/**
	* @brief Writes given string into cout or cerr (based on log_level) in the given color.
	* @param log_level The log_level is used to check if write in cout or cerr
	* @param log_color The color the string should be printed in.
	* @param string The string to write.
	*/
    void coloredCoutCerr(LogColor log_color, std::string string);

	static std::mutex m_mutex; //!< Mutex object used to protect data for thread safe use
	std::string m_log_prefix; //!< String prefix for log messages
	std::string m_file_prefix; //!< String prefic for file name
	bool m_compile_debug = false; //!< Set to true if compiled in debug mode
	bool m_log_to_file = true; //!< Use to toggle logging to file
	bool m_log_to_console = true; //!< Use to toggle logging to console
	Level m_log_level = Level::LVL_ALL; //!< Minimum log level to be logged
	bool m_use_color = true; //!< Use to toggle colored console output
	bool m_log_source = true; //!< Use to toggle logs source (filename:line) logging
};


template<typename ...Args>
inline void Log::log(Level log_level, const char * const filepath, std::string line, std::string in_formatmessage, Args ...args)
{
	//Do not log if this log entrys log_level is below m_log_level.
	if (log_level < m_log_level) return;

    // lock for logging from multiple threads
    std::lock_guard<std::mutex> guard(m_mutex);

	std::string log_message_pre;
	std::string log_message;
	std::string log_source = std::string(filepath) + ":" + line;
	//Append log_prefix
	log_message_pre.append(m_log_prefix);
	//Append log_level followed by tab
	log_message_pre.append("[" + levelToString(log_level) + "]\t");
	log_message.append(logHelper(in_formatmessage, args...));

	if (!log_source.empty()) log_source = "(" + log_source + ")\t";

	if (m_log_to_console) {
        coloredCoutCerr(levelToLogColor(log_level), log_message_pre);
		if (m_log_source)
            coloredCoutCerr(LogColor::GREY, log_source);
        std::cout << log_message << std::endl;
	}

	if (m_log_to_file) {
		if (!m_outstream) initOutstream(); // Init outstream if it hasn't been initiated yet.
		if (m_outstream.good()) {
			m_outstream << log_message_pre << log_source << log_message << std::endl;
		}
		else {
			std::cerr << "Error while writing to Log outstream" << std::endl;
		}
	}
}

template<typename T, typename ...Args>
inline std::string Log::logHelper(std::string format, T value, Args ...args)
{
	std::string log_message;
	for (size_t i = 0; i < format.size(); ++i){
		// if current char is a %
		if (format[i] == '%') {
			//check if next character is % too, if so we skip it.
			//using double %% means we escape % so we are still able to print % in log
			if (i < (format.size() - 1) && format[i + 1] == '%') {
				log_message.push_back('%');
				++i;
				continue;
			}
			//read current arg in stringstream and append to log_message
			std::stringstream stringstream;
			stringstream << value;
			log_message.append(stringstream.str());

			//recursive call loghelper with remaining format substring and remaining args
			if (i < (format.size() - 1)) {
				log_message.append(logHelper(format.substr(i + 1), args...));
			}

			break;
		}
		// otherwise insert current char
		else {
			log_message.push_back(format[i]);
		}
	}
	return log_message;
}

inline std::string Log::logHelper(std::string format)
{
	return format;
}
}
}

// include forward declared classes
//--------------------
#include "LogStream.h"
//--------------------

#endif
