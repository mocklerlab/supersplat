/*
 * File:   reads.cpp
 * Author: Douglas W Bryant Jr
 *
 * Created on December 28, 2008, 11:41 AM
 */

#include "Reads.h"

using namespace seqan;

Reads::Reads(char* p_isFastaReadsFile) {
    m_hShortReads   		= new __gnu_cxx::hash_map<const std::string, ShortReadStorageInfo>;
    m_iTotalReads   		= new int(0);
    m_iUniqueReads  		= new int(0);
    m_iErroredReads 		= new int(0);
    m_iReadsPassedFilter 	= new int(0);

    Reads::loadReadsData(p_isFastaReadsFile);
}

Reads::~Reads() {
    delete m_hShortReads;
    delete m_iTotalReads;
    delete m_iUniqueReads;
    delete m_iErroredReads;
}

bool Reads::loadReadsData(char* p_csReadsFile) {
	__utility::outputStatus("Loading reads file ...\n");

	MultiSeqFile multiSeqFile;
	open(multiSeqFile.concat, p_csReadsFile, OPEN_RDONLY);
	AutoSeqFormat format;
	guessFormat(multiSeqFile.concat, format);
	split(multiSeqFile, format);
	unsigned seqCount = length(multiSeqFile);
	CharString id;
	String<Dna5Q> seq;
	int l_iUniqueReads = 0, l_iErroredReads = 0;
	bool l_bForwardExistsInHash, l_bReverseExistsInHash;
	__gnu_cxx::hash_map<const std::string, ShortReadStorageInfo>::iterator l_hiShortReadsIterator;
	for(unsigned i = 0; i < seqCount; i++) {
		assignSeqId(id, multiSeqFile[i], format);
		assignSeq(seq, multiSeqFile[i], format);
		std::string header; assign(header, id);
		std::string sequence; assign(sequence, seq);

		std::transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper);
		for(unsigned i = 0; i < sequence.length(); i++) {
			if(sequence[i] != 'A' && sequence[i] != 'C' && sequence[i] != 'T' && sequence[i] != 'G') {
				l_iErroredReads++; break;
			}
		}

		m_hShortReads->find(sequence) == m_hShortReads->end() ? l_bForwardExistsInHash = false : l_bForwardExistsInHash = true;
		m_hShortReads->find(__utility::findReverseComplement(&sequence)) == m_hShortReads->end() ? l_bReverseExistsInHash = false : l_bReverseExistsInHash = true;
		if(!l_bForwardExistsInHash && !l_bReverseExistsInHash) {
			ShortReadStorageInfo l_dShortReadStorageInfoToInsert;
			l_dShortReadStorageInfoToInsert.ReadHeadersVector.push_back(header);
			m_hShortReads->insert(std::pair<std::string, ShortReadStorageInfo>(sequence, l_dShortReadStorageInfoToInsert));
			l_iUniqueReads++;
		} else if(l_bForwardExistsInHash && !l_bReverseExistsInHash) {
			l_hiShortReadsIterator = m_hShortReads->find(sequence);
			l_hiShortReadsIterator->second.ReadHeadersVector.push_back(header);
		} else if(!l_bForwardExistsInHash && l_bReverseExistsInHash) {
			l_hiShortReadsIterator = m_hShortReads->find(__utility::findReverseComplement(&sequence));
			l_hiShortReadsIterator->second.ReadHeadersVector.push_back(header);
		}
	}

	*m_iTotalReads 			= seqCount;
	*m_iUniqueReads 		= l_iUniqueReads;
	*m_iReadsPassedFilter 	= -1;
	*m_iErroredReads 		= l_iErroredReads;

	std::stringstream ss;
	std::string l_sTotalReads, l_sUniqueReads;
	ss << seqCount; l_sTotalReads = ss.str(); ss.str("");
	ss << l_iUniqueReads; l_sUniqueReads = ss.str(); ss.str("");
	__utility::outputStatus("\tDONE loading reads file. Loaded " + l_sTotalReads + " total reads comprising " + l_sUniqueReads + " unique sequences.\n");
    return true;
}

bool Reads::pruneReads(int p_iMinCopyNumber, int p_iMaxCopyNumber) {
	std::stringstream ss;
	std::string l_sMin, l_sMax;
	ss << p_iMinCopyNumber; l_sMin = ss.str(); ss.str("");
	ss << p_iMaxCopyNumber; l_sMax = ss.str(); ss.str("");
	std::string l_sStatus = "Pruning reads ";
	if(p_iMinCopyNumber != -1 && p_iMaxCopyNumber != -1) l_sStatus += "of copy number less than " + l_sMin + " and of copy number greater than " + l_sMax + " ...\n";
	else if(p_iMinCopyNumber != -1) l_sStatus += "of copy number less than " + l_sMin + " ...\n";
	else if(p_iMaxCopyNumber != -1) l_sStatus += "of copy number greater than " + l_sMax + " ...\n";
	__utility::outputStatus(l_sStatus);

	__gnu_cxx::hash_map<const std::string, ShortReadStorageInfo>::iterator it = m_hShortReads->begin();
	__gnu_cxx::hash_map<const std::string, ShortReadStorageInfo>::iterator bit;
	int l_iReadsRemoved = 0;
	int l_iReadsPassedFilter = 0;
	while(it != m_hShortReads->end()) {
		if((p_iMinCopyNumber != -1 && it->second.ReadHeadersVector.size() < p_iMinCopyNumber) || (p_iMaxCopyNumber != -1 && it->second.ReadHeadersVector.size() > p_iMaxCopyNumber)) {
			bit = it;
			it++;
			m_hShortReads->erase(bit);
			l_iReadsRemoved++;
		} else {
			it++;
			l_iReadsPassedFilter++;
		}
	}

	*m_iReadsPassedFilter = l_iReadsPassedFilter;
	int l_iToSplat = m_hShortReads->size();

	std::string l_sRemoved, l_sPassed, l_sToSplat;
	ss << l_iReadsRemoved; l_sRemoved = ss.str(); ss.str("");
	ss << l_iReadsPassedFilter; l_sPassed = ss.str(); ss.str("");
	ss << l_iToSplat; l_sToSplat = ss.str(); ss.str("");
	l_sStatus = "\tDONE pruning reads. " + l_sRemoved + " reads removed, " + l_sToSplat + " reads to splat.\n";
	__utility::outputStatus(l_sStatus);

	return true;
}
