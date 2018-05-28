#include "ArgParser.h"

KIRK::ArgParser::ArgMap KIRK::ArgParser::toMap(int argc, char *argv[])
{
	KIRK::ArgParser::ArgMap arg_map;

	if (argc % 2 != 1) {
		// args are always in pairs of -<key> <value> starting with one executable call.
		// So there has to be an odd number of args.
		LOG_ERROR("Argument count has to be odd. Is %.", argc);
		return arg_map;
	}

	for (int ai = 1; ai < argc; ai += 2)
	{
		char* key = argv[ai];
		char* value = argv[ai + 1];

		// Now:
		// key[0] should always be a dash.
		if (key[0] != '-') {
			LOG_ERROR("Argument key not valid. Has to begin with a dash, followed by a key character.");
			continue;
		}

		// key[1] is the according key char.
		arg_map.put(key[1], std::string(value));
	}

	return arg_map;
}