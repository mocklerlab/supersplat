/*******************************************************************************
supersplat
Douglas Bryant, Jr.
Weng-Keen Wong
Todd Mockler
Oregon State University,
Departments of Electrical Engineering and Computer Science and Botany and
Plant Pathology

This software is free for academic and public use. No warranty is offered nor
implied. Please contact Douglas Bryant with any questions, comments, or
concerns at <bryantjr AT eecs.oregonstate.edu>.

http://supersplat.cgrb.oregonstate.edu
******************************************************************************/

#ifndef _REFERENCE_H
#define	_REFERENCE_H

#include <iostream>
#include <stdio.h>
#include <sstream>
#include <algorithm>
#include "stdafx.h"
#include "math.h"

struct ReferenceIndexInfo
{
    int ReferenceNumber;
    long Position;
    const ReferenceStorageInfo* ReferenceInfo;
};

class Reference
{
public:
    Reference(std::vector<ReferenceStorageInfo>* p_vReferenceSequences, int p_iMinChunkLength, int p_iMaxIndexChunkSize, int p_iNumInternalIntronBases) ;
    virtual ~Reference() ;    
    // Public accessor functions  
    int getMinChunkLength() {return *m_ipMinChunkLength;} ;
    std::string getReferenceSequence() {return *m_spReferenceSequence;} ;
    //Public utility functions
    bool splatRead(std::string* p_spReadToSplat, std::vector<ShortReadFileInfo>& p_dSplatResults, const int p_iMinIntronLength, const int p_iMaxIntronLength) ;
private:
    // Private member variables
    const int* MAX_INDEXED_CHUNKSIZE;
    const int* NUM_FLANKING_BASES;
	const int* NUM_INTERNAL_INTRON_BASES;
    __gnu_cxx::hash_map<const std::string, std::vector<ReferenceIndexInfo> >* m_hReferenceKmerPositions;
    std::vector<ReferenceStorageInfo>* m_vReferenceSequences;
    std::vector<long>*  m_vReferenceOffsets;
    std::string*        m_spReferenceSequence;          // Reference genomic sequence.
    int*                m_ipMinChunkLength;             // Minimum length k-mer chunk.
    int*                m_ipReferenceSize;              // Size of the reference.
    // Private functions    
    bool indexAllPositions() ;
    int whichReference(long p_iPosition, int p_iChunkSize) ;
};

#endif	/* _REFERENCE_H */
