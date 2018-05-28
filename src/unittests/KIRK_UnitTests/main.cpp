#include <gtest/gtest.h>
#include "KIRK/Utils/Log.h"

/**
	Lists all available tests and asks user for input, which tests should be run.
*/
int main(int argc, char **argv)
{
	// Init and run googletest just to list all existing Test Cases.
	// No Testcases will be run with gtest_list_tests flag true.
	::testing::InitGoogleTest(&argc, argv);
	::testing::FLAGS_gtest_list_tests = true;
	RUN_ALL_TESTS();

	KIRK::LOG::Log::getInstance().setFilePrefix("UnitTests_");
	KIRK::LOG::Log::getInstance().setLogPrefix("[----------] ");

	// Init GoogleTest again
	::testing::InitGoogleTest(&argc, argv);
	// Let user decide which tests should be run
	std::cout << "Enter a filter to choose which tests should be run." << std::endl;
	std::cout << "	e.g. 'CPUDatastructureTest.*' to run all CPUDatastructureTest test cases." << std::endl;
	std::cout << "	or '*' to run all testcases.'" << std::endl;
	std::string filter;
	std::cin >> filter;
	std::cout << "Set Flag --gtest_filter=" << filter << std::endl;
	::testing::FLAGS_gtest_filter = filter;
	::testing::FLAGS_gtest_list_tests = false;

	//Run All Testcases matching the given filter.
	return RUN_ALL_TESTS();
}
