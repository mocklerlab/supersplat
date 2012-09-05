/*
 * StackParameters.h
 *
 *  Created on: Apr 26, 2011
 *      Author: Douglas W Bryant Jr
 */

#ifndef STACKPARAMETERS_H_
#define STACKPARAMETERS_H_

#include <string>

class StackParameters {
public:
	StackParameters();
	virtual ~StackParameters();
	void setReferenceFileName(std::string name) { *m_sReferenceFileName = name; }
	std::string* getReferenceFileName() { return m_sReferenceFileName; }
	void setSplatOutputFileName(std::string name) { *m_sSplatOutputFileName = name; }
	std::string* getSplatOutputFileName() { return m_sSplatOutputFileName; }
	void setReferenceNameToAnalyze(std::string name) { *m_sReferenceNameToAnalyze = name; }
	std::string* getReferenceNameToAnalyze() { return m_sReferenceNameToAnalyze; }
	void setStackOutputFileName(std::string name) { *m_sStackOutputFileName = name; }
	std::string* getStackOutputFileName() { return m_sStackOutputFileName; }
	void setMinReadCopyNumber(int num) { *m_iMinReadCopyNumber = num; }
	int* getMinReadCopyNumber() { return m_iMinReadCopyNumber; }
	void setMaxReadCopyNumber(int num) { *m_iMaxReadCopyNumber = num; }
	int* getMaxReadCopyNumber() { return m_iMaxReadCopyNumber; }
	void setMinNumDiffSequences(int num) { *m_iMinNumDiffSequences = num; }
	int* getMinNumDiffSequences() { return m_iMinNumDiffSequences; }
	bool verifyParameters() ;
private:
	std::string* 	m_sReferenceFileName;
	std::string* 	m_sSplatOutputFileName;
	std::string*	m_sReferenceNameToAnalyze;
	std::string*	m_sStackOutputFileName;
	int*			m_iMinReadCopyNumber;
	int*			m_iMaxReadCopyNumber;
	int*			m_iMinNumDiffSequences;
};

#endif /* STACKPARAMETERS_H_ */
