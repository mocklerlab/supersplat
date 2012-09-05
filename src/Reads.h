/*
 * File:   reads.h
 * Author: Douglas W Bryant Jr
 *
 * Created on December 28, 2008, 11:41 AM
 */

#ifndef _READS_H
#define	_READS_H

#include "stdafx.h"
#include <fstream>
#include <algorithm>
#include <seqan/file.h>

class Reads
{
public:
	Reads(char* p_isFastaReadsFile);
    virtual ~Reads();
    // Public accessor functions.
    const __gnu_cxx::hash_map<const std::string, ShortReadStorageInfo>* getShortReads() {return m_hShortReads;} ;
    int getTotalReads() {return *m_iTotalReads;} ;
    int getUniqueReads() {return *m_iUniqueReads;} ;
    int getErroredReads() {return *m_iErroredReads;} ;
    // Public utility functions.
    bool pruneReads(int p_iMinCopyNumber, int p_iMaxCopyNumber) ;
private:
    // Private global variables.
    __gnu_cxx::hash_map<const std::string, ShortReadStorageInfo>* 	m_hShortReads;
    int*                                                   			m_iTotalReads;
    int*                                                    		m_iUniqueReads;
    int*                                                   			m_iErroredReads;
    int*															m_iReadsPassedFilter;
    // Private functions.
    bool loadReadsData(char* p_csReadsFile) ;

};

#endif	/* _READS_H */
