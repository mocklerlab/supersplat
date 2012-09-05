/* 
 * File:   scheduler.cpp
 * Author: Douglas W Bryant Jr
 * 
 * Created on December 29, 2008, 2:53 PM
 */

#include "Scheduler.h"

Scheduler::Scheduler(std::ifstream& p_isFastaReferenceFile, std::ofstream& p_osOutputFile, Reads* p_dFastaReads, int p_iMinChunkSize, int p_iMinIntronSize, int p_iMaxIntronSize,
		int p_iMaxIndexChunkSize, int p_iNumInternalIntronBases, bool p_bCannonicalOnly, int p_iNumSimultaneousRefs, int p_iNumThreads,
		std::string l_sReferenceNameToIndex, bool p_bPrintHeaders) : m_osOutputFile(p_osOutputFile) {
    m_iNumReferences        	= new int(0);
    m_iMinChunkSize         	= new int(p_iMinChunkSize);
    m_vReferenceSequences   	= new std::vector<ReferenceStorageInfo>;
    m_dFastaReads           	= p_dFastaReads;
    m_iMinIntronSize        	= new int(p_iMinIntronSize);
    m_iMaxIntronSize        	= new int(p_iMaxIntronSize);
    m_iMaxIndexChunkSize    	= new int(p_iMaxIndexChunkSize);
    m_iNumInternalIntronBases 	= new int(p_iNumInternalIntronBases);
    m_iNumThreads				= new int(p_iNumThreads);
    m_bCannonicalOnly       	= new bool(p_bCannonicalOnly);
    m_bPrintHeaders				= new bool(p_bPrintHeaders);
    m_iNumSimultaneousRefs  	= new int(p_iNumSimultaneousRefs);

    if(l_sReferenceNameToIndex == "") Scheduler::loadReferenceData(p_isFastaReferenceFile);
    else Scheduler::loadReferenceData(p_isFastaReferenceFile, l_sReferenceNameToIndex);
}

Scheduler::~Scheduler() {
    delete m_iNumReferences;
    delete m_iMinChunkSize;
    delete m_vReferenceSequences;
    delete m_iMinIntronSize;
    delete m_iMaxIntronSize;
    delete m_iMaxIndexChunkSize;
    delete m_iNumInternalIntronBases;
    delete m_iNumThreads;
    delete m_bCannonicalOnly;
    delete m_bPrintHeaders;
    delete m_iNumSimultaneousRefs;
}

bool Scheduler::loadReferenceData(std::ifstream& p_isFastaReferenceFile) {
	__utility::outputStatus("Loading reference sequence(s) ...\n");

	int l_iNumReferences = 0;
    std::string l_sReadLine = "";
    std::string l_sReferenceHeader = "";
    std::string l_sReferenceChunk = "";
    while(getline(p_isFastaReferenceFile, l_sReadLine)) {
        if(l_sReadLine[0] == '>') {
            l_sReferenceHeader = l_sReadLine.substr(1);
        } else {
            l_sReferenceChunk = l_sReadLine;            

            while(p_isFastaReferenceFile.peek() != '>' && p_isFastaReferenceFile.peek() != EOF) {
                getline(p_isFastaReferenceFile, l_sReadLine);
                l_sReferenceChunk += l_sReadLine;
            }

            ReferenceStorageInfo l_dReferenceStorageInfo;
            l_dReferenceStorageInfo.Header = l_sReferenceHeader;
            l_dReferenceStorageInfo.ReferenceSequence = l_sReferenceChunk;
            m_vReferenceSequences->push_back(l_dReferenceStorageInfo);

            l_iNumReferences++;
        }
    }

    *m_iNumReferences = l_iNumReferences;

    std::stringstream ss;
	std::string l_sNumReferences, l_sStatus;
	ss << l_iNumReferences; l_sNumReferences = ss.str(); ss.str("");

	if(l_iNumReferences == 1) l_sStatus = "\tDONE loading the single reference sequence.\n";
	else l_sStatus = "\tDONE loading the " + l_sNumReferences + " reference sequences.\n";
    __utility::outputStatus(l_sStatus);

    return true;
}

bool Scheduler::loadReferenceData(std::ifstream& p_isFastaReferenceFile, std::string p_sReferenceName) {
	__utility::outputStatus("Loading the single reference " + p_sReferenceName + " ...\n");

    int l_iNumReferences = 0;
    std::string l_sReadLine = "";
    std::string l_sReferenceHeader = "";
    std::string l_sReferenceChunk = "";
    bool l_bFoundReference = false;
    while(getline(p_isFastaReferenceFile, l_sReadLine)) {
        if(l_sReadLine[0] == '>') {
            l_sReferenceHeader = l_sReadLine.substr(1);
        } else {
            l_sReferenceChunk = l_sReadLine;

            while(p_isFastaReferenceFile.peek() != '>' && p_isFastaReferenceFile.peek() != EOF) {
                getline(p_isFastaReferenceFile, l_sReadLine);
                l_sReferenceChunk += l_sReadLine;
            }

            if(l_sReferenceHeader == p_sReferenceName) {
            	l_bFoundReference = true;

				ReferenceStorageInfo l_dReferenceStorageInfo;
				l_dReferenceStorageInfo.Header = l_sReferenceHeader;
				l_dReferenceStorageInfo.ReferenceSequence = l_sReferenceChunk;
				m_vReferenceSequences->push_back(l_dReferenceStorageInfo);
            }
        }
    }

    if(l_bFoundReference) {
    	__utility::outputStatus("\tDONE loading reference " + p_sReferenceName + ".\n");
    	*m_iNumReferences = 1;
    } else {
    	__utility::outputStatus("\tERROR: Did not find reference " + p_sReferenceName + " - exiting.\n");
    	*m_iNumReferences = 0;
    }

    return l_bFoundReference;
}

void Scheduler::printData(std::vector<ShortReadFileInfo> p_dSplatResults, std::vector<std::string> p_vReadHeads) {
    std::vector<ShortReadFileInfo> l_dSplatResults = p_dSplatResults;
    std::vector<std::string> l_vReadHeads = p_vReadHeads;

    if(!p_dSplatResults.empty()) {
        // Print results to file
        for(std::vector<ShortReadFileInfo>::const_iterator it = p_dSplatResults.begin(); it != p_dSplatResults.end(); it++) {
            ShortReadFileInfo l_dTempShortReadFileInfo = *it;

            if(!(*m_bCannonicalOnly) || l_dTempShortReadFileInfo.FlankingSequence == "GT-AG" || l_dTempShortReadFileInfo.FlankingSequence == "GC-AG"
                    || l_dTempShortReadFileInfo.FlankingSequence == "AT-AC" || l_dTempShortReadFileInfo.FlankingSequence == "CT-AC"
                    || l_dTempShortReadFileInfo.FlankingSequence == "CT-GC" || l_dTempShortReadFileInfo.FlankingSequence == "GT-AT") {
                m_osOutputFile << l_dTempShortReadFileInfo.ChromosomeName << "\t";
                m_osOutputFile << l_dTempShortReadFileInfo.FlankingSequence << "\t";
                m_osOutputFile << l_dTempShortReadFileInfo.ChunkOneLength << "\t";
                m_osOutputFile << l_dTempShortReadFileInfo.ChunkTwoLength << "\t";
                m_osOutputFile << l_dTempShortReadFileInfo.IntronLength << "\t";
                m_osOutputFile << l_dTempShortReadFileInfo.ChunkOneBegin << "\t";
                m_osOutputFile << l_dTempShortReadFileInfo.ChunkOneEnd << "\t";
                m_osOutputFile << l_dTempShortReadFileInfo.ChunkTwoBegin << "\t";
                m_osOutputFile << l_dTempShortReadFileInfo.ChunkTwoEnd << "\t";
                m_osOutputFile << l_dTempShortReadFileInfo.ReadSequence << "\t";
                m_osOutputFile << l_vReadHeads.size() << "\t";

                if(*m_bPrintHeaders) {
					for(std::vector<std::string>::const_iterator iu = l_vReadHeads.begin(); iu != l_vReadHeads.end(); iu++) {
						if(iu != l_vReadHeads.begin()) m_osOutputFile << ",";
						m_osOutputFile << *iu;
					}
                }

				if(l_dTempShortReadFileInfo.DonorSequence != "" && l_dTempShortReadFileInfo.AcceptorSequence != "") {
					m_osOutputFile << "\t" << l_dTempShortReadFileInfo.DonorSequence << "\t";
					m_osOutputFile << l_dTempShortReadFileInfo.AcceptorSequence;
				} else if(l_dTempShortReadFileInfo.DonorSequence != "") {
					m_osOutputFile << "\t" << l_dTempShortReadFileInfo.DonorSequence;
				}

                m_osOutputFile << std::endl;
            }
        }
    }
}

bool Scheduler::run() {
	std::string* l_sCurrentReferenceSequence;
    int l_iNumReferences = m_vReferenceSequences->size(), l_iNumSimultaneousRefs = 0;

    if(*m_iNumSimultaneousRefs >= l_iNumReferences) l_iNumSimultaneousRefs = l_iNumReferences;
    else if(*m_iNumSimultaneousRefs <= 0) l_iNumSimultaneousRefs = 1;
    else l_iNumSimultaneousRefs = *m_iNumSimultaneousRefs;

    std::vector<ReferenceStorageInfo> l_dCurrentSimultaneousRefs;
    const __gnu_cxx::hash_map<const std::string, ShortReadStorageInfo>* reads = m_dFastaReads->getShortReads();
    __gnu_cxx::hash_map<const std::string, ShortReadStorageInfo>::const_iterator l_hiReadsIterator;
    std::vector<ReferenceStorageInfo>::iterator l_itRefsIterator = m_vReferenceSequences->begin();
    bool l_bReachedEndOfRefs = false;
    while(!l_bReachedEndOfRefs) {
        Reference* l_dCurrentReference;

        l_dCurrentSimultaneousRefs.clear();        
        while(!l_bReachedEndOfRefs && l_dCurrentSimultaneousRefs.size() < l_iNumSimultaneousRefs) {
            l_dCurrentSimultaneousRefs.push_back(*l_itRefsIterator);
            l_itRefsIterator++;
            if(l_itRefsIterator == m_vReferenceSequences->end()) l_bReachedEndOfRefs = true;
        }

        l_dCurrentReference = new Reference(&l_dCurrentSimultaneousRefs, *m_iMinChunkSize, *m_iMaxIndexChunkSize, *m_iNumInternalIntronBases);

        std::stringstream ss;
        int l_iNumReferences = l_dCurrentSimultaneousRefs.size(), l_iNumReads = reads->size();
		std::string l_sNumReferences, l_sNumReads, l_sStatus;
        ss << l_iNumReferences; l_sNumReferences = ss.str(); ss.str("");
        ss << l_iNumReads; l_sNumReads = ss.str(); ss.str("");
        if(l_iNumSimultaneousRefs == 1) l_sStatus = "Splatting " + l_sNumReads + " reads over reference " + l_dCurrentSimultaneousRefs.front().Header + " ...\n";
        else l_sStatus = "Splatting " + l_sNumReads + " reads over " + l_sNumReferences + "references ...\n";
        __utility::outputStatus(l_sStatus);

        int l_iTotalReadsPercent = 0, l_iTotalReadsProcessed = 0;
        int const l_iTenPercentOfTotalReads = reads->size() / 10;
        l_hiReadsIterator = reads->begin();
        std::vector<ShortReadFileInfo> p_dSplatResults;
        std::vector<std::string> l_vReadHeads;
        std::string l_sShortRead, l_sReverseComplement;
        int nthreads, tid, th_id;
        omp_set_num_threads(*m_iNumThreads);
		#define CHUNKSIZE 1
		#pragma omp parallel private(p_dSplatResults, l_vReadHeads, l_sShortRead, l_sReverseComplement, th_id)
        {
		#pragma omp for schedule(dynamic, CHUNKSIZE)
        for(int i = 0; i < reads->size(); i++) {
        	p_dSplatResults.clear();

			#pragma omp critical
        	{
        	l_sShortRead = l_hiReadsIterator->first;
        	l_vReadHeads = l_hiReadsIterator->second.ReadHeadersVector;
        	l_hiReadsIterator++;
        	}

        	l_dCurrentReference->splatRead(&l_sShortRead, p_dSplatResults, *m_iMinIntronSize, *m_iMaxIntronSize);

			l_sReverseComplement = __utility::findReverseComplement(&l_sShortRead);
			l_dCurrentReference->splatRead(&l_sReverseComplement, p_dSplatResults, *m_iMinIntronSize, *m_iMaxIntronSize);

			#pragma omp critical
			{
				l_iTotalReadsProcessed++;
				l_iTotalReadsPercent++;
				if(l_iTotalReadsPercent == l_iTenPercentOfTotalReads) {
					int l_iTotalPercentProcessed = l_iTotalReadsProcessed / l_iTenPercentOfTotalReads * 10;
					std::string l_sTotalPercentProcessed, l_sTotalReadsProcessed;
					ss << l_iTotalPercentProcessed; l_sTotalPercentProcessed = ss.str(); ss.str("");
					ss << l_iTotalReadsProcessed; l_sTotalReadsProcessed = ss.str(); ss.str("");
					l_sStatus = "\t" + l_sTotalPercentProcessed + "% of reads = " + l_sTotalReadsProcessed + " reads processed ...\n";
					__utility::outputStatus(l_sStatus);
					l_iTotalReadsPercent = 0;
				}

				Scheduler::printData(p_dSplatResults, l_vReadHeads);
			}
        }
        }

        if(l_iNumSimultaneousRefs == 1) l_sStatus = "\tDONE splatting all reads over reference " + l_dCurrentSimultaneousRefs.front().Header + ".\n";
		else l_sStatus = "\tDONE splatting all reads over " + l_sNumReferences + "references.\n";
        __utility::outputStatus(l_sStatus);

        delete l_dCurrentReference;
    }

    __utility::outputStatus("DONE.\n");
}
