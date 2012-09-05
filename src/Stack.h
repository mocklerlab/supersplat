/*
 * Stack.h
 *
 *  Created on: Apr 22, 2011
 *      Author: Douglas W Bryant Jr
 */

#ifndef STACK_H_
#define STACK_H_

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include "StackParameters.h"
#include "Stacker.h"

class Stack {
public:
	Stack(int argc, char** argv);
	virtual ~Stack();
private:
	void outputUsage() ;
	bool run();
	void printResults(std::vector<std::string>*) ;
	StackParameters* m_dParams;
};


#endif /* STACK_H_ */
