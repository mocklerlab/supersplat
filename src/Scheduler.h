/* 
 * File:   scheduler.h
 * Author: Douglas W Bryant Jr
 *
 * Created on December 29, 2008, 2:53 PM
 */

#ifndef _SCHEDULER_H
#define	_SCHEDULER_H

#include <fstream>
#include <iostream>
#include <omp.h>
#include <sstream>
#include "stdafx.h"
#include "Reference.h"
#include "Reads.h"

class Scheduler
{
public:
    Scheduler(std::ifstream& p_isFastaReferenceFile, std::ofstream& p_osOutputFile, Reads* p_dFastaReads, int p_iMinChunkSize, int p_iMinIntronSize, int p_iMaxIntronSize, 
            int p_iMaxIndexChunkSize, int p_iNumInternalIntronBases, bool p_bCannonicalOnly, int p_iNumSimultaneousRefs, int p_iNumThreads,
            std::string p_sNameOfReferenceToIndex, bool p_bPrintHeaders);
    virtual ~Scheduler();
    // Public accessor functions.
    int getNumReferences() {return *m_iNumReferences;} ;
    bool run();
private:
    // Private functions.
    bool loadReferenceData(std::ifstream& p_isFastaReferenceFile) ;
    bool loadReferenceData(std::ifstream& p_isFastaReferenceFile, std::string p_sNameOfReferenceToIndex) ;
    void printData(std::vector<ShortReadFileInfo> p_dSplatResults, std::vector<std::string> p_vReadHeads) ;
    // Private member variables.
    Reads* m_dFastaReads;
    int* m_iNumReferences;
    int* m_iMinChunkSize;
    int* m_iMinIntronSize;
    int* m_iMaxIntronSize;
    int* m_iMaxIndexChunkSize;
    int* m_iNumInternalIntronBases;
    int* m_iNumThreads;
    int* m_iNumSimultaneousRefs;
    bool* m_bCannonicalOnly;
    bool* m_bPrintHeaders;
    std::vector<ReferenceStorageInfo>* m_vReferenceSequences;
    std::ofstream& m_osOutputFile;
};

#endif	/* _SCHEDULER_H */
