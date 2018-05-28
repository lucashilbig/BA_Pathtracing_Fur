#ifndef ARG_PARSER_H
#define ARG_PARSER_H

#include <map>
#include <string>

#include "KIRK/Utils/Log.h"

namespace KIRK 
{
	/** General namespace for executable argument parsing. */
	namespace ArgParser 
	{
		/**
		This class is a special map<char, std::string> wrapper which introduces 
		default values and automatic value parsing to default types.
		*/
		class ArgMap
		{
		public:
			/** Empty default constructor. */
			ArgMap() {}

			/**
			Adds an argument key-value pair to the map. You can fetch it by ArgMap::get(char, ...).
			@param key The mapped key and also the argument key.
			@param value The argument value
			*/
			void put(char key, std::string value) 
			{
				m_args.insert(std::pair<char, std::string>(key, value));
			}

			/**
			Get the value of an argument by it's key char (e.g. -s "value" -> fetch "value" by calling get('s')). You can also give a default value.
			The returned object type, as well as the default value is determined at compile time, so you can call <code>int i=get('k', 100)</code> as well as
			<code>std::string s=get('s', std::string("string"))</code>.
			@param key The argument to fetch the value of.
			@param default_value The value which will be returned if the argument key does not exist.
			@return The converted argument value for the given key or default_value, if the key does not exist.
			*/
			template<typename BasicType>
            BasicType get(char key, BasicType default_value = BasicType()) const
			{
				try
				{
					std::string val_str = m_args.at(key);
					std::istringstream convert(val_str);

					//Use a stringstream to convert from string to whatever type is given.
					//Obviously that does not work for that many types, but the most
					//important ones like int, std::string, double and float should work.
					BasicType result;
					convert >> result;
					return result;
				}
				catch (std::out_of_range e)
				{
					return default_value;
				}
			}

		private:
			std::map<char, std::string> m_args;		//!< The map which maps the argument keys to the argument value strings.
		};

		/**
		Creates an ArgMap from the given args. Just pass the ones from <code>int main(int argc, char *argv[])</code>.<br/>
		The argument count has to be an odd number (that means that the count of parseable segments is even as the executable argument will be discarded).
		@param argc The size of argv or the argument count.
		@param argv All passed arguments including the executable call.
		@return An ArgMap containing all parsed arguments as char keys and string values.
		*/
		ArgMap toMap(int argc, char *argv[]);
	};
}

#endif // !ARG_PARSER_H
