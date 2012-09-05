/*
 * Stacker.h
 *
 *  Created on: Apr 26, 2011
 *      Author: Douglas W Bryant Jr
 */

#ifndef STACKER_H_
#define STACKER_H_

#include <vector>
#include <sstream>
#include <fstream>
#include <string>
#include <algorithm>
#include "stdafx.h"
#include "StackParameters.h"

struct SplatResultInfo {
	std::string ReferenceName;
	int DonorSite;
	int AcceptorSite;
	std::string FlankingDinucleotides;
	int LeftChunkLength;
	int RightChunkLength;
	int IntronLength;
	std::string Sequence;
	int ReadCopyNum;
	float Score;

	bool operator<(const SplatResultInfo& a) const {
		return LeftChunkLength < a.LeftChunkLength;
	}
};

struct StackerReferenceInfo {
	std::string* ReferenceHeader;
	// This ugly thing goes: DonarPosition -> (AcceptorPosition -> (List of splat results that hit those two positions))
	__gnu_cxx::hash_map<int, __gnu_cxx::hash_map<int, std::vector<SplatResultInfo> >* >* DonarToAcceptorHash;
	int* NumReferencePositions;
};

class Stacker {
public:
	Stacker(StackParameters*);
	virtual ~Stacker();
	bool initReferenceInfos(std::string) ;
	bool stackSplatOutput(std::string) ;
	std::vector<std::string>* getStackedResults() ;
private:
	bool PassCuttoffs(std::vector<SplatResultInfo>*) ;
	std::string formatOutput(std::vector<SplatResultInfo>) ;
	__gnu_cxx::hash_map<const std::string, StackerReferenceInfo>* 	m_hStackerReferenceInfos;
	std::string* 													m_sReferenceNameToAnalyze;
	int* 															m_iMinReadCopyNumber;
	int* 															m_iMaxReadCopyNumber;
	int* 															m_iMinNumDiffSequences;
	std::vector<std::string>*										m_vResults;
};

#endif /* STACKER_H_ */
