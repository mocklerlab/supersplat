/* 
 * File:   stdafx.h
 * Author: Douglas W Bryant Jr
 *
 * Created on December 12, 2008, 10:47 PM
 */

#ifndef _STDAFX_H
#define	_STDAFX_H

#include <iostream>
#include <ext/hash_map>
#include <time.h>
#include <string.h>

namespace __gnu_cxx {
template<> struct hash<const std::string> {
	size_t operator()(const std::string& s) const {
		size_t h = 0;
		std::string::const_iterator p, p_end;
		for (p = s.begin(), p_end = s.end(); p != p_end; ++p) {
			h = 31 * h + (*p);
		}
		return h;
	}
};
}

struct ReferenceStorageInfo {
	std::string Header;
	std::string ReferenceSequence;
};

struct ShortReadStorageInfo {
	std::vector<std::string> ReadHeadersVector;
};

struct ShortReadFileInfo {
    std::string ChromosomeName;
    std::string FlankingSequence;
    int ChunkOneLength;
    int ChunkTwoLength;
    int IntronLength;
    int ChunkOneBegin;
    int ChunkOneEnd;
    int ChunkTwoBegin;
    int ChunkTwoEnd;
    std::string ReadSequence;
    int NumberOfReads;
    std::vector<std::string> ReadHeads;
	std::string DonorSequence;
	std::string AcceptorSequence;
};

namespace __utility {
static std::string findReverseComplement(std::string* p_sInputRead) {
	const std::string l_sInputRead = *p_sInputRead;
	std::string l_sReverseComplementRead = "";
	std::string::const_iterator l_siReadIterator = l_sInputRead.end();

	while (l_siReadIterator != l_sInputRead.begin()) {
		l_siReadIterator--;
		switch (*l_siReadIterator) {
		case 'a':
		case 'A':
			l_sReverseComplementRead += "T";
			break;
		case 'c':
		case 'C':
			l_sReverseComplementRead += "G";
			break;
		case 't':
		case 'T':
			l_sReverseComplementRead += "A";
			break;
		case 'g':
		case 'G':
			l_sReverseComplementRead += "C";
			break;
		default:
			//throw exception
			break;
		}
	}
	return l_sReverseComplementRead;
}

static void outputStatus(std::string p_sPrintMe) {
	time_t tim=time(NULL);
	char *s=ctime(&tim);
	s[strlen(s)-1]=0;        // remove \n

	std::cout << "[" << s << "] " << p_sPrintMe;
	return;
}
}

#endif	/* _STDAFX_H */
