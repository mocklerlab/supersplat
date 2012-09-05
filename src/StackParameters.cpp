/*
 * StackParameters.cpp
 *
 *  Created on: Apr 26, 2011
 *      Author: Douglas W Bryant Jr
 */

#include "StackParameters.h"

StackParameters::StackParameters() {
	m_sReferenceFileName = new std::string("");
	m_sSplatOutputFileName = new std::string("");
	m_sReferenceNameToAnalyze = new std::string("");
	m_sStackOutputFileName = new std::string("");
	m_iMinReadCopyNumber = new int(-1);
	m_iMaxReadCopyNumber = new int(-1);
	m_iMinNumDiffSequences = new int(-1);
}

StackParameters::~StackParameters() {
	delete m_sReferenceFileName;
	delete m_sSplatOutputFileName;
	delete m_sReferenceNameToAnalyze;
	delete m_sStackOutputFileName;
	delete m_iMinReadCopyNumber;
	delete m_iMaxReadCopyNumber;
	delete m_iMinNumDiffSequences;
}

bool StackParameters::verifyParameters() {
	bool l_bParametersGood = true;
	if(*m_sReferenceFileName == "" || *m_sSplatOutputFileName == "" || *m_iMinReadCopyNumber == -1
			|| *m_iMinNumDiffSequences == -1) l_bParametersGood = false;
	if(l_bParametersGood && *m_sStackOutputFileName == "") *m_sStackOutputFileName = *m_sReferenceFileName + ".stack";
	return l_bParametersGood;
}
