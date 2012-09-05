/*
 * Splat.h
 *
 *  Created on: Apr 21, 2011
 *      Author: Douglas W Bryant Jr
 */

#ifndef SPLAT_H_
#define SPLAT_H_

#include <fstream>
#include <iostream>
#include <stdio.h>
#include "stdafx.h"
#include "Reference.h"
#include "Reads.h"
#include "Scheduler.h"

class Splat {
public:
	Splat(int argc, char** argv) ;
	virtual ~Splat();
private:
	void outputUsage() ;

};

#endif /* SPLAT_H_ */
