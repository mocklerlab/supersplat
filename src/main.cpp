/*
 * main.cpp
 *
 *  Created on: Apr 21, 2011
 *      Author: Douglas W Bryant Jr
 */

#include <iostream>
#include <string>
#include "Splat.h"
#include "Stack.h"

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "2.0"
#endif
#ifndef CONTACT
#define CONTACT "Doug Bryant <bryantjr@eecs.oregonstate.edu>"
#endif

static bool usage() {
	std::cout << std::endl;
	std::cout << "Program: supersplat (splice junction finder from RNA-seq data)" << std::endl;
	std::cout << "Version: " << PACKAGE_VERSION << std::endl;
	std::cout << "Contact: " << CONTACT << std::endl;
	std::cout << std::endl;
	std::cout << "Usage: supersplat [command] <options>" << std::endl << std::endl;
	std::cout << "Command: splat    find splice junctions from RNA-seq data" << std::endl;
	std::cout << "         stack    determine trusted splice junctions from splat output" << std::endl;
	std::cout << std::endl;
	return true;
}

int main(int argc, char** argv) {
	if(argc < 2) return usage();
	if(strcmp(argv[1], "splat") == 0) Splat* splat = new Splat(argc-1, argv+1);
	else if(strcmp(argv[1], "stack") == 0) Stack* stack = new Stack(argc-1, argv+1);
	else {
		std::cout << "[main] unrecognized command " << argv[1] << std::endl;
		return 1;
	}
	return 0;
}
