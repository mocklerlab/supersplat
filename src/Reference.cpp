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

#include "Reference.h"

Reference::Reference(std::vector<ReferenceStorageInfo>* p_vReferenceSequences, int p_iMinChunkLength, int p_iMaxIndexChunkSize, int p_iNumInternalIntronBases) {
    MAX_INDEXED_CHUNKSIZE           = new int(p_iMaxIndexChunkSize);
    NUM_FLANKING_BASES              = new int(2);
	NUM_INTERNAL_INTRON_BASES		= new int(p_iNumInternalIntronBases);
    m_hReferenceKmerPositions       = new __gnu_cxx::hash_map<const std::string, std::vector<ReferenceIndexInfo> >;
    m_ipMinChunkLength              = new int(p_iMinChunkLength);
    m_vReferenceSequences           = p_vReferenceSequences;
    m_vReferenceOffsets             = new std::vector<long> (p_vReferenceSequences->size());

    Reference::indexAllPositions();
}

Reference::~Reference() {
    delete MAX_INDEXED_CHUNKSIZE;
    delete NUM_FLANKING_BASES;
    delete m_hReferenceKmerPositions;
    delete m_ipMinChunkLength;
    delete m_vReferenceOffsets;
}

bool Reference::indexAllPositions() {
	std::stringstream ss;
	std::string l_sNumToIndex, l_sStatus;
	int l_iNumToIndex = m_vReferenceSequences->size();
	ss << l_iNumToIndex; l_sNumToIndex = ss.str(); ss.str("");
	if(l_iNumToIndex == 1) {
		std::string l_sReferenceName = m_vReferenceSequences->at(0).Header;
		l_sStatus = "Indexing reference " + l_sReferenceName + " ...\n";
	} else l_sStatus = "Indexing " + l_sNumToIndex + " references ...\n";
	__utility::outputStatus(l_sStatus);

	bool l_bError             = false;
    int l_iMinChunkLength     = *m_ipMinChunkLength;
    long l_iCurrentReferenceOffset = 0;
    std::vector<ReferenceStorageInfo>::const_iterator l_viReferencesVectorIterator;
    std::string l_sKMer;
    __gnu_cxx::hash_map<const std::string, std::vector<ReferenceIndexInfo> >::iterator l_hiHashIterator;
    
    int l_iNumRefsIndexed = 0;
    for(l_viReferencesVectorIterator = m_vReferenceSequences->begin(); l_viReferencesVectorIterator != m_vReferenceSequences->end(); l_viReferencesVectorIterator++) {
        const ReferenceStorageInfo* l_dpCurrentReferenceStorageInfo = &(*l_viReferencesVectorIterator);
        const std::string* l_spCurrentReferenceName = &(l_dpCurrentReferenceStorageInfo->Header);
        const std::string* l_spCurrentReferenceSequence = &(l_dpCurrentReferenceStorageInfo->ReferenceSequence);

        int l_iCurrentReferenceLength = l_spCurrentReferenceSequence->length();

        for(long i = 0; i < l_iCurrentReferenceLength - l_iMinChunkLength; i++) {
            // Take care of l_iMinChunkLength first, in case it's > MAX_CHUNK_SIZE_TO_INDEX.
            l_sKMer = l_spCurrentReferenceSequence->substr(i, l_iMinChunkLength);
            l_hiHashIterator = m_hReferenceKmerPositions->find(l_sKMer);

            ReferenceIndexInfo l_dTempInitialReferenceIndexInfo;
            l_dTempInitialReferenceIndexInfo.Position = i + l_iCurrentReferenceOffset;
            l_dTempInitialReferenceIndexInfo.ReferenceInfo = &(*l_viReferencesVectorIterator);
            l_dTempInitialReferenceIndexInfo.ReferenceNumber = l_iNumRefsIndexed;

            if(l_hiHashIterator == m_hReferenceKmerPositions->end()) {
                std::vector<ReferenceIndexInfo> l_vTempVector;
                l_vTempVector.push_back(l_dTempInitialReferenceIndexInfo);
                m_hReferenceKmerPositions->insert(std::pair<const std::string, std::vector<ReferenceIndexInfo> >(l_sKMer, l_vTempVector));
            } else {
                l_hiHashIterator->second.push_back(l_dTempInitialReferenceIndexInfo);
            }

            for(long j = l_iMinChunkLength + 1; j <= *MAX_INDEXED_CHUNKSIZE; j++) {
                if(l_iCurrentReferenceLength - i > j) {
                    l_sKMer = l_spCurrentReferenceSequence->substr(i, j);
                    l_hiHashIterator = m_hReferenceKmerPositions->find(l_sKMer);

                    ReferenceIndexInfo l_dTempReferenceIndexInfo;
                    l_dTempReferenceIndexInfo.Position = i + l_iCurrentReferenceOffset;
                    l_dTempReferenceIndexInfo.ReferenceInfo = &(*l_viReferencesVectorIterator);
                    l_dTempReferenceIndexInfo.ReferenceNumber = l_iNumRefsIndexed;

                    if(l_hiHashIterator == m_hReferenceKmerPositions->end()) {
                        std::vector<ReferenceIndexInfo> l_vTempVector;
                        l_vTempVector.push_back(l_dTempReferenceIndexInfo);
                        m_hReferenceKmerPositions->insert(std::pair<const std::string, std::vector<ReferenceIndexInfo> >(l_sKMer, l_vTempVector));
                    } else {
                        l_hiHashIterator->second.push_back(l_dTempReferenceIndexInfo);
                    }
                }
            }
        }

        m_vReferenceOffsets->at(l_iNumRefsIndexed) = l_iCurrentReferenceOffset;
        l_iCurrentReferenceOffset += (l_iCurrentReferenceLength - 1);

        l_iNumRefsIndexed++;
    }
    
    if(l_iNumToIndex == 1) {
    	std::string l_sReferenceName = m_vReferenceSequences->at(0).Header;
    	l_sStatus = "\tDONE indexing reference " + l_sReferenceName + ".\n";
    }
    else l_sStatus = "\tDONE indexing " + l_sNumToIndex + " references.\n";
    __utility::outputStatus(l_sStatus);

    return l_bError;
}

int Reference::whichReference(long p_iPosition, int p_iChunkSize) {
    std::vector<ReferenceStorageInfo>::const_iterator l_viReferencesVectorIterator = m_vReferenceSequences->begin();
    long l_iCurrentReferenceOffset = 0;
    int l_iFoundReference = -1;
    bool l_bFoundReference = false;

    int l_iCurrentReference = 0;
    const ReferenceStorageInfo* l_dpCurrentReferenceStorageInfo;
    int l_iCurrentReferenceSize;
    while(!l_bFoundReference && l_viReferencesVectorIterator != m_vReferenceSequences->end()) {
        l_dpCurrentReferenceStorageInfo = &(*l_viReferencesVectorIterator);
        l_iCurrentReferenceSize = l_dpCurrentReferenceStorageInfo->ReferenceSequence.size();

        if(l_iCurrentReferenceOffset <= p_iPosition && p_iPosition <= l_iCurrentReferenceOffset + l_iCurrentReferenceSize
                && p_iPosition + p_iChunkSize <= l_iCurrentReferenceOffset + l_iCurrentReferenceSize) {
            l_iFoundReference = l_iCurrentReference;
            l_bFoundReference = true;
        }

        l_iCurrentReferenceOffset += (l_iCurrentReferenceSize - 1);
        l_iCurrentReference++;
    }

    return l_iFoundReference;
}

bool Reference::splatRead(std::string* p_spReadToSplat, std::vector<ShortReadFileInfo>& p_dSplatResults, const int p_iMinIntronLength, const int p_iMaxIntronLength) {
    bool l_bError                       = false;
    const int l_iMinChunkLength         = *m_ipMinChunkLength;
    const int l_iReadToSplatLength      = p_spReadToSplat->length();
    const int l_iMinIntronLenth         = p_iMinIntronLength;
    const int l_iMaxIntronLength        = p_iMaxIntronLength;
    const std::string l_sReadToSplat    = *p_spReadToSplat;

    if(l_iReadToSplatLength < l_iMinChunkLength * 2) { // Read must have at least two chunks of at least length MinChunkLength in order to search for both chunks in the indexed reference.
        l_bError = true;
    } else {
        __gnu_cxx::hash_map<const std::string, std::vector<ReferenceIndexInfo> >::const_iterator l_hiChunkOneIterator;
        __gnu_cxx::hash_map<const std::string, std::vector<ReferenceIndexInfo> >::const_iterator l_hiChunkTwoIterator;

        std::string l_sChunkOne;
        std::string l_sChunkTwo;
        std::string l_sChunkOneKMer;
        std::string l_sChunkTwoKMer;

        const std::vector<ReferenceIndexInfo>* l_vChunkOnePositions;
        const std::vector<ReferenceIndexInfo>* l_vChunkTwoPositions;
        std::vector<ReferenceIndexInfo>::const_iterator l_dChunkOnePositionsIterator;
        std::vector<ReferenceIndexInfo>::const_iterator l_dChunkTwoPositionsIterator;

        int l_iChunkOneSize = l_iMinChunkLength;
        int l_iChunkTwoSize = l_iReadToSplatLength - l_iChunkOneSize;

        l_sChunkOne = l_sReadToSplat.substr(0, l_iChunkOneSize);
        l_sChunkTwo = l_sReadToSplat.substr(l_iChunkOneSize, l_iChunkTwoSize);
        while(l_iChunkTwoSize >= l_iMinChunkLength) {
            if(l_iChunkOneSize <= *MAX_INDEXED_CHUNKSIZE) l_sChunkOneKMer = l_sChunkOne;
            else l_sChunkOneKMer = l_sChunkOne.substr(0, *MAX_INDEXED_CHUNKSIZE);

            if(l_iChunkTwoSize <= *MAX_INDEXED_CHUNKSIZE) l_sChunkTwoKMer = l_sChunkTwo;
            else l_sChunkTwoKMer = l_sChunkTwo.substr(0, *MAX_INDEXED_CHUNKSIZE);

            l_hiChunkOneIterator = m_hReferenceKmerPositions->find(l_sChunkOneKMer);
            l_hiChunkTwoIterator = m_hReferenceKmerPositions->find(l_sChunkTwoKMer);
            if(l_hiChunkOneIterator != m_hReferenceKmerPositions->end() && l_hiChunkTwoIterator != m_hReferenceKmerPositions->end()) {
                l_vChunkOnePositions = &l_hiChunkOneIterator->second;
                l_vChunkTwoPositions = &l_hiChunkTwoIterator->second;

                ReferenceIndexInfo l_dTempReferenceIndexInfoChunkOne;
                ReferenceIndexInfo l_dTempReferenceIndexInfoChunkTwo;
                std::string l_sChunkOneReferenceName;
                std::string l_sChunkTwoReferenceName;
                long l_lChunkOnePosition;
                long l_lChunkTwoPosition;
                int l_iReferenceNumber;
                for(l_dChunkTwoPositionsIterator = l_vChunkTwoPositions->begin(); l_dChunkTwoPositionsIterator != l_vChunkTwoPositions->end(); l_dChunkTwoPositionsIterator++) {
                    l_dChunkOnePositionsIterator = l_vChunkOnePositions->begin();

                    l_dTempReferenceIndexInfoChunkOne = *l_dChunkOnePositionsIterator;
                    l_dTempReferenceIndexInfoChunkTwo = *l_dChunkTwoPositionsIterator;
                    l_sChunkOneReferenceName = l_dTempReferenceIndexInfoChunkOne.ReferenceInfo->Header;
                    l_sChunkTwoReferenceName = l_dTempReferenceIndexInfoChunkTwo.ReferenceInfo->Header;
                    l_lChunkOnePosition = l_dTempReferenceIndexInfoChunkOne.Position;
                    l_lChunkTwoPosition = l_dTempReferenceIndexInfoChunkTwo.Position;
                    l_iReferenceNumber = l_dTempReferenceIndexInfoChunkOne.ReferenceNumber;

                    while(l_lChunkTwoPosition - l_lChunkOnePosition - l_iChunkOneSize + 1 >= l_iMinIntronLenth && l_dChunkOnePositionsIterator != l_vChunkOnePositions->end()) {
                        if(l_lChunkTwoPosition - l_lChunkOnePosition + l_iChunkOneSize - 1 <= l_iMaxIntronLength) {
                            if(l_sChunkOneReferenceName == l_sChunkTwoReferenceName) {
                                const ReferenceStorageInfo* l_dCurrentReferenceStorageInfo = l_dTempReferenceIndexInfoChunkOne.ReferenceInfo;
                                const std::string* l_spCurrentReferenceSequence = &(l_dCurrentReferenceStorageInfo->ReferenceSequence);
                                const std::string* l_spCurrentReferenceName = &(l_dCurrentReferenceStorageInfo->Header);

                                long l_lChunkOneStart = l_lChunkOnePosition - m_vReferenceOffsets->at(l_iReferenceNumber);
                                long l_lChunkTwoStart = l_lChunkTwoPosition - m_vReferenceOffsets->at(l_iReferenceNumber);

                                if(l_lChunkTwoStart + l_iChunkTwoSize <= l_spCurrentReferenceSequence->size()) {
                                    bool l_bChunkOneMatches;
                                    bool l_bChunkTwoMatches;

                                    if(l_iChunkOneSize <= *MAX_INDEXED_CHUNKSIZE) {
                                        l_bChunkOneMatches = true;
                                    } else {
                                        l_spCurrentReferenceSequence->substr(l_lChunkOneStart, l_iChunkOneSize) == l_sChunkOne ? l_bChunkOneMatches = true : l_bChunkOneMatches = false;
                                    }

                                    if(l_iChunkTwoSize <= *MAX_INDEXED_CHUNKSIZE) {
                                        l_bChunkTwoMatches = true;
                                    } else {
                                        l_spCurrentReferenceSequence->substr(l_lChunkTwoStart, l_iChunkTwoSize) == l_sChunkTwo ? l_bChunkTwoMatches = true : l_bChunkTwoMatches = false;
                                    }

                                    if(l_bChunkOneMatches && l_bChunkTwoMatches) {
                                        ShortReadFileInfo l_dShortReadInfoToReturn;

                                        l_dShortReadInfoToReturn.ChromosomeName     = *l_spCurrentReferenceName;

										std::string l_sLeftFlankingSequence         = l_spCurrentReferenceSequence->substr(l_lChunkOneStart + l_iChunkOneSize, *NUM_FLANKING_BASES);
                                        std::string l_sRightFlankingSequence        = l_spCurrentReferenceSequence->substr(l_lChunkTwoStart - *NUM_FLANKING_BASES, *NUM_FLANKING_BASES);
										int l_iIntronLength 						= l_lChunkTwoStart - l_lChunkOneStart - l_iChunkOneSize + 1;

										std::string l_sDonorSequence				= "";
										std::string l_sAcceptorSequence				= "";
										if(*NUM_INTERNAL_INTRON_BASES > 0 && *NUM_INTERNAL_INTRON_BASES <= l_iIntronLength) {
											std::string l_sDonorSequenceJunction		= "(" + l_sLeftFlankingSequence + ")";
											std::string l_sAcceptorSequenceJunction		= "(" + l_sRightFlankingSequence + ")";
											l_sDonorSequence							= l_sDonorSequenceJunction + l_spCurrentReferenceSequence->substr(l_lChunkOneStart + l_iChunkOneSize + *NUM_FLANKING_BASES, *NUM_INTERNAL_INTRON_BASES);
											l_sAcceptorSequence							= l_spCurrentReferenceSequence->substr(l_lChunkTwoStart - *NUM_FLANKING_BASES - *NUM_INTERNAL_INTRON_BASES, *NUM_INTERNAL_INTRON_BASES) + l_sAcceptorSequenceJunction;
										} else if(*NUM_INTERNAL_INTRON_BASES > l_iIntronLength) {
											l_sDonorSequence = "(Intron too short.)";
										}

                                        l_dShortReadInfoToReturn.FlankingSequence   = l_sLeftFlankingSequence + "-" + l_sRightFlankingSequence;
                                        l_dShortReadInfoToReturn.ChunkOneLength     = l_iChunkOneSize;
                                        l_dShortReadInfoToReturn.ChunkTwoLength     = l_iChunkTwoSize;
                                        l_dShortReadInfoToReturn.IntronLength       = l_iIntronLength;
                                        l_dShortReadInfoToReturn.ChunkOneBegin      = l_lChunkOneStart + 1; // The following four (+ 1)s are to enforce 1-indexing in output.
                                        l_dShortReadInfoToReturn.ChunkOneEnd        = l_lChunkOneStart + l_iChunkOneSize - 1 + 1;
                                        l_dShortReadInfoToReturn.ChunkTwoBegin      = l_lChunkTwoStart + 1;
                                        l_dShortReadInfoToReturn.ChunkTwoEnd        = l_lChunkTwoStart + l_iChunkTwoSize - 1 + 1;
                                        l_dShortReadInfoToReturn.ReadSequence       = l_sReadToSplat;
										l_dShortReadInfoToReturn.DonorSequence 		= l_sDonorSequence;
										l_dShortReadInfoToReturn.AcceptorSequence	= l_sAcceptorSequence;

                                        p_dSplatResults.push_back(l_dShortReadInfoToReturn);
                                    }
                                }
                            }
                        }

                        l_dChunkOnePositionsIterator++;

                        if(l_dChunkOnePositionsIterator != l_vChunkOnePositions->end()) {
                            l_dTempReferenceIndexInfoChunkOne = *l_dChunkOnePositionsIterator;
                            l_sChunkOneReferenceName = l_dTempReferenceIndexInfoChunkOne.ReferenceInfo->Header;
                            l_lChunkOnePosition = l_dTempReferenceIndexInfoChunkOne.Position;
                            l_iReferenceNumber = l_dTempReferenceIndexInfoChunkOne.ReferenceNumber;
                        }
                    }
                }
            }

            l_sChunkOne = l_sReadToSplat.substr(0, ++l_iChunkOneSize);
            l_sChunkTwo = l_sReadToSplat.substr(l_iChunkOneSize, --l_iChunkTwoSize);
        }
    }

    return l_bError;
}
