#include "Log.h"

#include <sstream>
#include <mutex>
#include <stdexcept>

// include windows h here because it defines some wired stuff
#if LOG_PLATFORM_WINDOWS
#define NOMINMAX
#include <windows.h>
#endif

namespace KIRK {
namespace LOG {

Log* Log::m_instance = nullptr;
std::mutex Log::m_mutex;

Log& Log::getInstance()
{
	static Cleanup cleanup;
	std::lock_guard<std::mutex> guard(m_mutex);
	if (m_instance == nullptr)
		m_instance = new Log();
	return *m_instance;
}
Log::Cleanup::~Cleanup()
{
	std::lock_guard<std::mutex> guard(Log::m_mutex);
	delete Log::m_instance;
	Log::m_instance = nullptr;
}

Log::Log()
{
	// Try to disable console output if it's not compiled in Debug Config.
	// Should work for all OS.
#ifdef NDEBUG
	m_log_to_console = true;
	m_log_to_file = true;
	m_log_source = false;
	m_log_level = Level::LVL_INFO;
#else
	m_log_to_console = true;
	m_log_to_file = true;
	m_log_source = true;
	m_log_level = Level::LVL_ALL;
	m_compile_debug = true;
#endif

	initFilename();
}

Log::~Log()
{
	m_outstream.close();
}

void Log::initOutstream()
{
	if (m_log_to_file) {
		if (m_outstream) m_outstream.close();
		m_outstream.open(m_filename, std::ios_base::app);
		if (!m_outstream.good()) {
			throw std::runtime_error("Unable to open Log file output stream.");
		}
	}
	else {
		if (m_outstream) m_outstream.close();
	}
}

void Log::initFilename()
{
	time_t t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
	tm *timeinfo = localtime(&t);

	char time[100];
	strftime(time, 100, "%Y-%m-%d_%H-%M-%S", timeinfo);

	std::stringstream stringstream;
	stringstream << LOGS_PATH << '/' << m_file_prefix << time << (m_compile_debug ? "_Debug" : "") << ".log";

	m_filename = stringstream.str();
}

std::string Log::levelToString(Level level)
{
	switch (level) {
	case Level::LVL_ALL:
		return "ALL";
	case Level::LVL_DEBUG:
		return "DEBUG";
	case Level::LVL_INFO:
		return "INFO";
	case Level::LVL_WARN:
		return "WARN";
	case Level::LVL_ERROR:
		return "ERROR";
	case Level::LVL_INPUT:
		return "INPUT";
	case Level::LVL_NONE:
		return "NONE";
	default:
		return "NONE";
	}
}

Log::LogColor Log::levelToLogColor(Level level)
{
	switch (level) {
	case Level::LVL_DEBUG:
		return Log::LogColor::MAGENTA;
	case Level::LVL_INFO:
		return Log::LogColor::GREEN;
	case Level::LVL_WARN:
		return Log::LogColor::YELLOW;
	case Level::LVL_ERROR:
		return Log::LogColor::RED;
	case Level::LVL_INPUT:
		return Log::LogColor::BLUE;
	default:
		return Log::LogColor::WHITE;
	}
}

#if LOG_PLATFORM_WINDOWS
unsigned short Log::logColorToColorAttr(Log::LogColor log_color)
{
	switch (log_color) {
	case LogColor::BLACK:
		return 0;
	case LogColor::DARK_BLUE:
		return FOREGROUND_BLUE;
	case LogColor::DARK_CYAN:
		return FOREGROUND_GREEN | FOREGROUND_BLUE;
	case LogColor::DARK_GREEN:
		return FOREGROUND_GREEN;
	case LogColor::DARK_GREY:
		return FOREGROUND_BLUE | FOREGROUND_GREEN | FOREGROUND_RED;
	case LogColor::DARK_MAGENTA:
		return FOREGROUND_BLUE | FOREGROUND_RED;
	case LogColor::DARK_RED:
		return FOREGROUND_RED;
	case LogColor::DARK_YELLOW:
		return FOREGROUND_GREEN | FOREGROUND_RED;
	case LogColor::BLUE:
		return FOREGROUND_INTENSITY | FOREGROUND_BLUE;
	case LogColor::CYAN:
		return FOREGROUND_INTENSITY | FOREGROUND_GREEN | FOREGROUND_BLUE;
	case LogColor::GREEN:
		return FOREGROUND_INTENSITY | FOREGROUND_GREEN;
	case LogColor::GREY:
		return FOREGROUND_INTENSITY;
	case LogColor::MAGENTA:
		return FOREGROUND_INTENSITY | FOREGROUND_BLUE | FOREGROUND_RED;
	case LogColor::RED:
		return FOREGROUND_INTENSITY | FOREGROUND_RED;
	case LogColor::YELLOW:
		return FOREGROUND_INTENSITY | FOREGROUND_GREEN | FOREGROUND_RED;
	case LogColor::WHITE:
		return FOREGROUND_INTENSITY | FOREGROUND_BLUE | FOREGROUND_GREEN | FOREGROUND_RED;
	default:
		return 0;
	}
}
#endif

const char* Log::logColorToColorAnsi(Log::LogColor log_color)
{
	switch (log_color) {
	case LogColor::BLACK:
		return "30";
	case LogColor::DARK_BLUE:
		return "34";
	case LogColor::DARK_CYAN:
		return "36";
	case LogColor::DARK_GREEN:
		return "32";
	case LogColor::DARK_GREY:
		return "30;1";
	case LogColor::DARK_MAGENTA:
		return "35";
	case LogColor::DARK_RED:
		return "31";
	case LogColor::DARK_YELLOW:
		return "33";
	case LogColor::BLUE:
		return "34;1";
	case LogColor::CYAN:
		return "36;1";
	case LogColor::GREEN:
		return "32;1";
	case LogColor::GREY:
		return "37";
	case LogColor::MAGENTA:
		return "35;1";
	case LogColor::RED:
		return "31;1";
	case LogColor::YELLOW:
		return "33;1";
	case LogColor::WHITE:
		return "37;1";
	default:
		return "37;1";
	}
}

void Log::coloredCoutCerr(LogColor log_color, std::string string)
{
	if (!m_use_color)
    {
        std::cout << string;
		return;
	}
#if LOG_PLATFORM_WINDOWS
	HANDLE h_std = GetStdHandle(STD_OUTPUT_HANDLE);

	CONSOLE_SCREEN_BUFFER_INFO cs_buffer_info;
	GetConsoleScreenBufferInfo(h_std, &cs_buffer_info);
	WORD old_color_attrs = cs_buffer_info.wAttributes;

	SetConsoleTextAttribute(h_std, logColorToColorAttr(log_color));
	std::cout << string;
	SetConsoleTextAttribute(h_std, old_color_attrs);
#else
		std::cout << "\033[0;" << logColorToColorAnsi(log_color) << 'm' << string << "\033[m";
#endif
}

void Log::setLogToFile(bool log_to_file)
{
	std::lock_guard<std::mutex> guard(m_mutex);
	m_log_to_file = log_to_file;
}

void Log::setLogToConsole(bool log_to_console)
{
	std::lock_guard<std::mutex> guard(m_mutex);
	m_log_to_console = log_to_console;
}

void Log::setLogLevel(Level log_level)
{
	std::lock_guard<std::mutex> guard(m_mutex);
	m_log_level = log_level;
}

void Log::setLogPrefix(std::string log_prefix)
{
	std::lock_guard<std::mutex> guard(m_mutex);
	m_log_prefix = log_prefix;
}
void Log::setFilePrefix(std::string file_prefix)
{
	std::lock_guard<std::mutex> guard(m_mutex);
	m_file_prefix = file_prefix;
	initFilename();
	if (m_outstream) initOutstream(); //If outstream already exists init it again
}
void Log::setUseColor(bool use_color)
{
	std::lock_guard<std::mutex> guard(m_mutex);
	m_use_color = use_color;
}
void Log::setLogSource(bool log_source)
{
	std::lock_guard<std::mutex> guard(m_mutex);
	m_log_source = log_source;
}
}
}
