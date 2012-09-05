/*
 * Stacker.cpp
 *
 *  Created on: Apr 26, 2011
 *      Author: Douglas W Bryant Jr
 */

#include "Stacker.h"

Stacker::Stacker(StackParameters* p_dParams) {
	m_sReferenceNameToAnalyze 	= p_dParams->getReferenceNameToAnalyze();
	m_iMinReadCopyNumber 		= p_dParams->getMinReadCopyNumber();
	m_iMaxReadCopyNumber 		= p_dParams->getMaxReadCopyNumber();
	m_iMinNumDiffSequences 		= p_dParams->getMinNumDiffSequences();
	m_hStackerReferenceInfos	= new __gnu_cxx::hash_map<const std::string, StackerReferenceInfo>();
	m_vResults					= new std::vector<std::string>();
}

Stacker::~Stacker() {
	delete m_sReferenceNameToAnalyze;
	delete m_iMinReadCopyNumber;
	delete m_iMaxReadCopyNumber;
	delete m_iMinNumDiffSequences;
	delete m_hStackerReferenceInfos;
	delete m_vResults;
}

bool Stacker::initReferenceInfos(std::string p_sReferenceFileName) {
	__utility::outputStatus("Loading reference sequence(s) ...\n");

	std::ifstream l_dFastaReferenceFile;
	l_dFastaReferenceFile.open(p_sReferenceFileName.c_str());

	int l_iNumReferences = 0;
    std::string l_sReadLine = "", l_sReferenceHeader = "", l_sReferenceChunk = "";
    while(getline(l_dFastaReferenceFile, l_sReadLine)) {
        if(l_sReadLine[0] == '>') { l_sReferenceHeader = l_sReadLine.substr(1); }
        else {
            l_sReferenceChunk = l_sReadLine;
            while(l_dFastaReferenceFile.peek() != '>' && l_dFastaReferenceFile.peek() != EOF) {
                getline(l_dFastaReferenceFile, l_sReadLine);
                l_sReferenceChunk += l_sReadLine;
            }
            StackerReferenceInfo l_dStackerReferenceInfo;
            l_dStackerReferenceInfo.ReferenceHeader = new std::string(l_sReferenceHeader);
            l_dStackerReferenceInfo.DonarToAcceptorHash = new
            		__gnu_cxx::hash_map<int, __gnu_cxx::hash_map<int, std::vector<SplatResultInfo> >* >();
            int NumRefPos = l_sReferenceChunk.length();
            l_dStackerReferenceInfo.NumReferencePositions = new int(NumRefPos);
            m_hStackerReferenceInfos->insert(std::pair<const std::string, StackerReferenceInfo>
            	(l_sReferenceHeader, l_dStackerReferenceInfo));
            l_iNumReferences++;
        }
    }
    std::stringstream ss;
	std::string l_sNumReferences, l_sStatus;
	ss << l_iNumReferences; l_sNumReferences = ss.str(); ss.str("");
	if(l_iNumReferences == 1) l_sStatus = "\tDONE loading the single reference sequence.\n";
	else l_sStatus = "\tDONE loading the " + l_sNumReferences + " reference sequences.\n";
    __utility::outputStatus(l_sStatus);
    l_dFastaReferenceFile.close();
    return true;
}

bool Stacker::stackSplatOutput(std::string p_sSplatOutputFileName) {
	__utility::outputStatus("Stacking splat output ...\n");
	std::ifstream l_dSplatOutputFile;
	l_dSplatOutputFile.open(p_sSplatOutputFileName.c_str());
	__gnu_cxx::hash_map<const std::string, std::vector<SplatResultInfo> >* l_dReadSJHash = new __gnu_cxx::hash_map<const std::string, std::vector<SplatResultInfo> >();
	__gnu_cxx::hash_map<const std::string, std::vector<SplatResultInfo> >::iterator l_dReadSJHashIterator;

	std::string l_sReadLine;
	while(getline(l_dSplatOutputFile, l_sReadLine)) {
		std::stringstream l_dStringStream(l_sReadLine);
		std::istream_iterator<std::string> l_dSSIterator(l_dStringStream);
		std::istream_iterator<std::string> l_dEndIterator;
		std::vector<std::string> l_vTokens(l_dSSIterator, l_dEndIterator);

		std::string l_sRefName 			= l_vTokens.at(0);
		std::string l_sDonorAcceptor	= l_vTokens.at(1);
		int l_iLeftChunkLength 			= atoi(l_vTokens.at(2).c_str());
		int l_iRightChunkLength 		= atoi(l_vTokens.at(3).c_str());
		int l_iIntronLength 			= atoi(l_vTokens.at(4).c_str());
		int	l_iDonorSite				= atoi(l_vTokens.at(6).c_str());
		int l_iAcceptorSite				= atoi(l_vTokens.at(7).c_str());
		std::string l_sSequence 		= l_vTokens.at(9);
		int l_iNumReadCopies			= atoi(l_vTokens.at(10).c_str());

		SplatResultInfo l_dTempSplatRI;
		l_dTempSplatRI.ReferenceName			= l_sRefName;
		l_dTempSplatRI.DonorSite				= l_iDonorSite;
		l_dTempSplatRI.AcceptorSite				= l_iAcceptorSite;
		l_dTempSplatRI.FlankingDinucleotides 	= l_sDonorAcceptor;
		l_dTempSplatRI.LeftChunkLength 			= l_iLeftChunkLength;
		l_dTempSplatRI.RightChunkLength 		= l_iRightChunkLength;
		l_dTempSplatRI.IntronLength 			= l_iIntronLength;
		l_dTempSplatRI.Sequence 				= l_sSequence;
		l_dTempSplatRI.ReadCopyNum 				= l_iNumReadCopies;

		l_dReadSJHashIterator = l_dReadSJHash->find(l_sSequence);
		if(l_dReadSJHashIterator == l_dReadSJHash->end()) {
			std::vector<SplatResultInfo> l_dTempVector;
			l_dTempVector.push_back(l_dTempSplatRI);
			l_dReadSJHash->insert(std::pair<const std::string, std::vector<SplatResultInfo> >(l_sSequence, l_dTempVector));
		} else {
			l_dReadSJHashIterator->second.push_back(l_dTempSplatRI);
		}
	}

	// We now score and retain only the best scoring location for each read:
	std::vector<SplatResultInfo> l_dFinalSRIs;
	l_dReadSJHashIterator = l_dReadSJHash->begin();
	while(l_dReadSJHashIterator != l_dReadSJHash->end()) {
		std::vector<SplatResultInfo>* l_dSRIs = &(l_dReadSJHashIterator->second);
		std::vector<SplatResultInfo>::const_iterator l_dSRIsIterator = l_dSRIs->begin();
		SplatResultInfo l_dBestSRI;
		while(l_dSRIsIterator != l_dSRIs->end()) {
			float l_iScore = 0.0, l_iHighScore = 0.0;
			float l_fLeftChunkLength, l_fRightChunkLength, l_fFlankLengthRatio;
			if(l_dSRIsIterator->FlankingDinucleotides == "GT-AG" || l_dSRIsIterator->FlankingDinucleotides == "CT-AC") l_iScore += 20;
			else if(l_dSRIsIterator->FlankingDinucleotides == "GC-AG" || l_dSRIsIterator->FlankingDinucleotides == "CT-GC") l_iScore += 7;
			else if(l_dSRIsIterator->FlankingDinucleotides == "AT-AC" || l_dSRIsIterator->FlankingDinucleotides == "GT-AT") l_iScore += 4;
			l_fLeftChunkLength = (float)l_dSRIsIterator->LeftChunkLength; l_fRightChunkLength = (float)l_dSRIsIterator->RightChunkLength;
			l_fFlankLengthRatio = l_fLeftChunkLength < l_fRightChunkLength ? l_fLeftChunkLength / l_fRightChunkLength : l_fRightChunkLength / l_fLeftChunkLength;
			l_iScore += l_fFlankLengthRatio;

			if(l_iScore > l_iHighScore) {
				l_iHighScore = l_iScore;
				l_dBestSRI = *l_dSRIsIterator;
				l_dBestSRI.Score = l_iScore;
			}

			l_dSRIsIterator++;
		}

		l_dFinalSRIs.push_back(l_dBestSRI);
		l_dReadSJHashIterator++;
	}

	// Now stack results:
	std::vector<SplatResultInfo>::const_iterator l_dFinalSRIsIterator = l_dFinalSRIs.begin();
	while(l_dFinalSRIsIterator != l_dFinalSRIs.end()) {
		SplatResultInfo l_dTempSplatRI = *l_dFinalSRIsIterator;
		std::string l_sRefName = l_dFinalSRIsIterator->ReferenceName;
		int l_iDonorSite = l_dFinalSRIsIterator->DonorSite;
		int l_iAcceptorSite = l_dFinalSRIsIterator->AcceptorSite;

		__gnu_cxx::hash_map<const std::string, StackerReferenceInfo>::iterator l_dSRIit = m_hStackerReferenceInfos->find(l_sRefName);
		if(l_dSRIit != m_hStackerReferenceInfos->end()) {
			StackerReferenceInfo l_dTempSRI = l_dSRIit->second;
			__gnu_cxx::hash_map<int, __gnu_cxx::hash_map<int, std::vector<SplatResultInfo> >* >* l_dTempDAHash = l_dTempSRI.DonarToAcceptorHash;
			__gnu_cxx::hash_map<int, __gnu_cxx::hash_map<int, std::vector<SplatResultInfo> >* >::iterator l_dDAHit = l_dTempDAHash->find(l_iDonorSite);
			if(l_dDAHit == l_dTempDAHash->end()) {
				// Add donor and acceptor:
				__gnu_cxx::hash_map<int, std::vector<SplatResultInfo> >* l_dTempInternalHash = new __gnu_cxx::hash_map<int, std::vector<SplatResultInfo> >();
				std::vector<SplatResultInfo> l_dTempSplatResultInfosVector;
				l_dTempSplatResultInfosVector.push_back(l_dTempSplatRI);
				l_dTempInternalHash->insert(std::pair<int, std::vector<SplatResultInfo> >(l_iAcceptorSite, l_dTempSplatResultInfosVector));
				l_dSRIit->second.DonarToAcceptorHash->insert(std::pair<int, __gnu_cxx::hash_map<int, std::vector<SplatResultInfo> >* >(l_iDonorSite, l_dTempInternalHash));
			} else {
				// l_dDAHit != end, so donor site found. Now search for acceptor site ...
				__gnu_cxx::hash_map<int, std::vector<SplatResultInfo> >* l_dInternalHash = l_dDAHit->second;
				__gnu_cxx::hash_map<int, std::vector<SplatResultInfo> >::iterator l_dInternalHashit = l_dInternalHash->find(l_iAcceptorSite);
				if(l_dInternalHashit == l_dInternalHash->end()) {
					// Add acceptor:
					std::vector<SplatResultInfo> l_dTempSplatResultInfosVector;
					l_dTempSplatResultInfosVector.push_back(l_dTempSplatRI);
					l_dDAHit->second->insert(std::pair<int, std::vector<SplatResultInfo> >(l_iAcceptorSite, l_dTempSplatResultInfosVector));
				} else {
					// Add to list of results at existing acceptor:
					l_dInternalHashit->second.push_back(l_dTempSplatRI);
				}
			}
		}

		l_dFinalSRIsIterator++;
	}

	/*
	while(getline(l_dSplatOutputFile, l_sReadLine)) {
		std::stringstream l_dStringStream(l_sReadLine);
		std::istream_iterator<std::string> l_dSSIterator(l_dStringStream);
		std::istream_iterator<std::string> l_dEndIterator;
		std::vector<std::string> l_vTokens(l_dSSIterator, l_dEndIterator);

		std::string l_sRefName 			= l_vTokens.at(0);
		std::string l_sDonorAcceptor	= l_vTokens.at(1);
		int l_iLeftChunkLength 			= atoi(l_vTokens.at(2).c_str());
		int l_iRightChunkLength 		= atoi(l_vTokens.at(3).c_str());
		int l_iIntronLength 			= atoi(l_vTokens.at(4).c_str());
		int	l_iDonorSite				= atoi(l_vTokens.at(6).c_str());
		int l_iAcceptorSite				= atoi(l_vTokens.at(7).c_str());
		std::string l_sSequence 		= l_vTokens.at(9);
		int l_iNumReadCopies			= atoi(l_vTokens.at(10).c_str());

		// Search for reference name ...
		__gnu_cxx::hash_map<const std::string, StackerReferenceInfo>::iterator l_dSRIit = m_hStackerReferenceInfos->find(l_sRefName);
		if(l_dSRIit != m_hStackerReferenceInfos->end()) {
			// Reference found so set up splat results info to insert:
			SplatResultInfo l_dTempSplatRI;
			l_dTempSplatRI.ReferenceName			= l_sRefName;
			l_dTempSplatRI.DonorSite				= l_iDonorSite;
			l_dTempSplatRI.AcceptorSite				= l_iAcceptorSite;
			l_dTempSplatRI.FlankingDinucleotides 	= l_sDonorAcceptor;
			l_dTempSplatRI.LeftChunkLength 			= l_iLeftChunkLength;
			l_dTempSplatRI.RightChunkLength 		= l_iRightChunkLength;
			l_dTempSplatRI.IntronLength 			= l_iIntronLength;
			l_dTempSplatRI.Sequence 				= l_sSequence;
			l_dTempSplatRI.ReadCopyNum 				= l_iNumReadCopies;

			// Search for donor site within reference ...
			StackerReferenceInfo l_dTempSRI = l_dSRIit->second;
			__gnu_cxx::hash_map<int, __gnu_cxx::hash_map<int, std::vector<SplatResultInfo> >* >* l_dTempDAHash = l_dTempSRI.DonarToAcceptorHash;
			__gnu_cxx::hash_map<int, __gnu_cxx::hash_map<int, std::vector<SplatResultInfo> >* >::iterator l_dDAHit = l_dTempDAHash->find(l_iDonorSite);
			if(l_dDAHit == l_dTempDAHash->end()) {
				// Add donor and acceptor:
				__gnu_cxx::hash_map<int, std::vector<SplatResultInfo> >* l_dTempInternalHash = new __gnu_cxx::hash_map<int, std::vector<SplatResultInfo> >();
				std::vector<SplatResultInfo> l_dTempSplatResultInfosVector;
				l_dTempSplatResultInfosVector.push_back(l_dTempSplatRI);
				l_dTempInternalHash->insert(std::pair<int, std::vector<SplatResultInfo> >(l_iAcceptorSite, l_dTempSplatResultInfosVector));
				l_dSRIit->second.DonarToAcceptorHash->insert(std::pair<int, __gnu_cxx::hash_map<int, std::vector<SplatResultInfo> >* >(l_iDonorSite, l_dTempInternalHash));
			} else {
				// l_dDAHit != end, so donor site found. Now search for acceptor site ...
				__gnu_cxx::hash_map<int, std::vector<SplatResultInfo> >* l_dInternalHash = l_dDAHit->second;
				__gnu_cxx::hash_map<int, std::vector<SplatResultInfo> >::iterator l_dInternalHashit = l_dInternalHash->find(l_iAcceptorSite);
				if(l_dInternalHashit == l_dInternalHash->end()) {
					// Add acceptor:
					std::vector<SplatResultInfo> l_dTempSplatResultInfosVector;
					l_dTempSplatResultInfosVector.push_back(l_dTempSplatRI);
					l_dDAHit->second->insert(std::pair<int, std::vector<SplatResultInfo> >(l_iAcceptorSite, l_dTempSplatResultInfosVector));
				} else {
					// Add to list of results at existing acceptor:
					l_dInternalHashit->second.push_back(l_dTempSplatRI);
				}
			}
		}
	}
	*/

	delete l_dReadSJHash;
	l_dSplatOutputFile.close();
	__utility::outputStatus("\tDONE stacking splat output.\n");
	return true;
}

std::vector<std::string>* Stacker::getStackedResults() {
	if(m_vResults->size() == 0) {
		__gnu_cxx::hash_map<const std::string, StackerReferenceInfo>::iterator l_dSRIit =
				m_hStackerReferenceInfos->begin();
		while(l_dSRIit != m_hStackerReferenceInfos->end()) {
			StackerReferenceInfo l_dStackerRefInfo = l_dSRIit->second;

			std::string* ref_name = l_dSRIit->second.ReferenceHeader;

			__gnu_cxx::hash_map<int, __gnu_cxx::hash_map<int, std::vector<SplatResultInfo> >* >::iterator
				l_dDonarAcceptorHashIT;
			int* l_iRefLength = l_dStackerRefInfo.NumReferencePositions;
			for(int i = 1; i <= *l_iRefLength; i++) {
				l_dDonarAcceptorHashIT = l_dStackerRefInfo.DonarToAcceptorHash->find(i);
				if(l_dDonarAcceptorHashIT != l_dStackerRefInfo.DonarToAcceptorHash->end()) {
					__gnu_cxx::hash_map<int, std::vector<SplatResultInfo> >* l_dInnerHash = l_dDonarAcceptorHashIT->second;
					__gnu_cxx::hash_map<int, std::vector<SplatResultInfo> >::iterator l_dInnerHashIT = l_dInnerHash->begin();
					while(l_dInnerHashIT != l_dInnerHash->end()) {
						if(PassCuttoffs(&(l_dInnerHashIT->second)))
							m_vResults->push_back(Stacker::formatOutput(l_dInnerHashIT->second));
						l_dInnerHashIT++;
					}
				}
			}
			l_dSRIit++;
		}
		if(m_vResults->size() == 0) m_vResults->push_back("No splice junctions passed filter criteria.\n");
	}
	return m_vResults;
}

std::string Stacker::formatOutput(std::vector<SplatResultInfo> p_vSplatResultInfosVector) {
	std::sort(p_vSplatResultInfosVector.begin(), p_vSplatResultInfosVector.end());
	std::stringstream ss;
	std::string l_sRefName = p_vSplatResultInfosVector.at(0).ReferenceName;
	std::string l_sFankingDinucleotides = p_vSplatResultInfosVector.at(0).FlankingDinucleotides;
	int l_iDonarSite = p_vSplatResultInfosVector.at(0).DonorSite;
	std::string l_sDonarSite; ss << l_iDonarSite; l_sDonarSite = ss.str(); ss.str("");
	int l_iAcceptorSite = p_vSplatResultInfosVector.at(0).AcceptorSite;
	std::string l_sAcceptorSite; ss << l_iAcceptorSite; l_sAcceptorSite = ss.str(); ss.str("");
	int l_iIntronLength = p_vSplatResultInfosVector.at(0).IntronLength;
	std::string l_sIntronLength; ss << l_iIntronLength; l_sIntronLength = ss.str(); ss.str("");

	std::string l_sOutput = "@" + l_sRefName + "\t" + l_sFankingDinucleotides + "\t" + l_sDonarSite + "\t" +
			l_sAcceptorSite + "\t" + l_sIntronLength + "\n";

	for(int i = 0; i < p_vSplatResultInfosVector.size(); i++) {
		std::string l_sReadCopyNum; ss << p_vSplatResultInfosVector.at(i).ReadCopyNum;
			l_sReadCopyNum = ss.str(); ss.str("");
		std::string l_sSequence = p_vSplatResultInfosVector.at(i).Sequence;
		std::string l_sLeftChunkLength; ss << p_vSplatResultInfosVector.at(i).LeftChunkLength;
			l_sLeftChunkLength = ss.str(); ss.str("");
		std::string l_sRightChunkLength; ss << p_vSplatResultInfosVector.at(i).RightChunkLength;
			l_sRightChunkLength = ss.str(); ss.str("");

		l_sOutput += l_sReadCopyNum + "\t" + l_sSequence + "\t" + l_sLeftChunkLength + "\t" +
				l_sRightChunkLength + "\n";
	}

	return l_sOutput;
}

bool Stacker::PassCuttoffs(std::vector<SplatResultInfo>* p_vSplatResultInfos) {
	bool l_bPassed = false;
	int l_iNumPassCopyCuttoff 	= 0;
	for(int i = 0; i < p_vSplatResultInfos->size(); i++) {
		int l_iReadCopyNum = p_vSplatResultInfos->at(i).ReadCopyNum;
		if(l_iReadCopyNum >= *m_iMinReadCopyNumber &&
				(*m_iMaxReadCopyNumber == -1 || l_iReadCopyNum <= *m_iMaxReadCopyNumber))
			l_iNumPassCopyCuttoff++;
	}

	if(l_iNumPassCopyCuttoff >= *m_iMinNumDiffSequences)
		l_bPassed = true;

	return l_bPassed;
}
