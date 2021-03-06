/* 
 * File:   SuperSplat.cpp
 * Author: Douglas W Bryant Jr
 *
 * Created on December 12, 2008, 10:58 AM
 */

#include <fstream>
#include <iostream>
#include <stdio.h>
#include "stdafx.h"
#include "Reference.h"
#include "Reads.h"
#include "Scheduler.h"

//int supersplat_splat(int argc, char** argv) ;

void outputUsage() {
    std::cout << "Usage: supersplat" << std::endl;
    std::cout << "\t-r <Fasta-formatted reference file.> (REQD)" << std::endl;
    std::cout << "\t-f <Fasta-formatted reads file.> (REQD)" << std::endl;
    std::cout << "\t-c <Minimum read chunk size - must be greater than 3.> (REQD)" << std::endl;
    std::cout << "\t-n <Minimum intron size - must be greater than 0.> (REQD)" << std::endl;
    std::cout << "\t-x <Maximum intron size - must be greater than minimum intron size.> (REQD)" << std::endl;
    std::cout << "\t-i <Maximum index chunk size - must be greater than minimum read chunk size.> (Optional, default = 12)" << std::endl;
    std::cout << "\t-t <Number of threads used for splatting.> (Optional, default = 1)" << std::endl;
    std::cout << "\t-s <Number of references to index simultaneously, 0 for all at once.> (Optional, default = 1)" << std::endl;
    std::cout << "\t-o <0 for all results, 1 for only canonical results.> (Optional, default = 0)" << std::endl;
    std::cout << "\t-e <Output file name.> (Optional, default = name_of_reads_file.supersplat)" << std::endl;
    std::cout << "Additional parameters:" << std::endl;
    std::cout << "\t-B <Number of internal flanking intron bases to output.> (Optional, default 0 = no output for this parameter)" << std::endl;
    std::cout << "\t-R <Name of reference over which to restrict analysis> (Optional, default = no such restriction, applies only for -s=1)" << std::endl;
    std::cout << "\t-H <Controls printing of read headers to output, 0 suppresses headers> (Optional, default = 1)" << std::endl;
    std::cout << "\t-N <Minimum read copy number> (Optional, restricts splatting to reads duplicated at least this many times)" << std::endl;
    std::cout << "\t-X <Maximum read copy number> (Optional, restricts splatting to reads duplicated at most this many times)" <<std::endl << std::endl;
}

int supersplat_splat(int argc, char** argv) {
    std::cout << std::endl << "supersplat version 2.0 - build date: April 18, 2011." << std::endl << std::endl;

    if(argc < 11) {
        outputUsage();
    } else {
        int l_iMinChunkSize 					= -1;
        int l_iMinIntronSize 					= -1;
        int l_iMaxIntronSize 					= -1;
        int l_iMaxIndexChunkSize 				= 12;
        int l_iNumInternalIntronBases 			= 0;
        int l_iNumThreads 						= 1;
        int l_iNumSimultaneousRefs 				= 0;
        int l_iMinReadCopyNumber				= -1;
        int l_iMaxReadCopyNumber				= -1;
        bool l_bCannonicalOnly 					= false;
        bool l_bPrintHeaders 					= true;
        std::string l_sOutputFileNameCmdLine	= "";
        std::string l_sReferenceNameToIndex 	= "";
        char* l_sFastaReadsFile;
        char* l_sFastaReferenceFile;
        std::ifstream l_dFastaReferenceFile;
        std::ofstream l_dOutputFile;

        char* parameter;
        for(int i = 1; i < argc; i++) {
            parameter = argv[i];
            if(parameter[0] == '-') {
                switch(parameter[1]) {
                    case 'r': l_sFastaReferenceFile = argv[++i]; break;
                    case 'f': l_sFastaReadsFile = argv[++i]; break;
                    case 'e': l_sOutputFileNameCmdLine = argv[++i]; break;
                    case 'c': l_iMinChunkSize = atoi(argv[++i]); break;
                    case 'n': l_iMinIntronSize = atoi(argv[++i]); break;
                    case 'x': l_iMaxIntronSize = atoi(argv[++i]); break;
                    case 'i': l_iMaxIndexChunkSize = atoi(argv[++i]); break;
                    case 't': l_iNumThreads = atoi(argv[++i]); break;
                    case 'o': l_bCannonicalOnly = atoi(argv[++i]); break;
                    case 's': l_iNumSimultaneousRefs = atoi(argv[++i]); break;
                    case 'B': l_iNumInternalIntronBases = atoi(argv[++i]); break;
                    case 'R': l_sReferenceNameToIndex = argv[++i]; break;
                    case 'H': l_bPrintHeaders = atoi(argv[++i]); break;
                    case 'N': l_iMinReadCopyNumber = atoi(argv[++i]); break;
                    case 'X': l_iMaxReadCopyNumber = atoi(argv[++i]); break;
                    default:
                        std::cout << "Error in " << argv[i] << " " << argv[i+1] << std::endl << std::endl;
                        break;
                }
            }
        }

        if(l_sFastaReadsFile == "" || l_sFastaReferenceFile == "" || l_iMinChunkSize < 3 || l_iMaxIndexChunkSize < l_iMinChunkSize
        		|| l_iMinIntronSize < 1 || l_iMaxIntronSize < l_iMinIntronSize || (l_sReferenceNameToIndex != "" && l_iNumSimultaneousRefs > 1)) {
            outputUsage();
        } else {
        	std::string l_sOutputFileName = "";
			std::string l_stdsFastaReadsFile(l_sFastaReadsFile);

			if(l_sOutputFileNameCmdLine == "") l_sOutputFileName = l_stdsFastaReadsFile + ".supersplat";
			else l_sOutputFileName = l_sOutputFileNameCmdLine;

            l_dFastaReferenceFile.open(l_sFastaReferenceFile);
            l_dOutputFile.open(l_sOutputFileName.c_str(), std::ios::app | std::ios::out);

            Reads* l_dReads = new Reads(l_sFastaReadsFile);
            if(l_iMinReadCopyNumber != -1 || l_iMaxReadCopyNumber != -1) l_dReads->pruneReads(l_iMinReadCopyNumber, l_iMaxReadCopyNumber);

            Scheduler* l_dScheduler = new Scheduler(l_dFastaReferenceFile, l_dOutputFile, l_dReads, l_iMinChunkSize, l_iMinIntronSize, l_iMaxIntronSize,
            		l_iMaxIndexChunkSize, l_iNumInternalIntronBases, l_bCannonicalOnly, l_iNumSimultaneousRefs, l_iNumThreads, l_sReferenceNameToIndex, l_bPrintHeaders);

            if(l_dScheduler->getNumReferences() > 0) l_dScheduler->run();

            delete l_dScheduler;
            delete l_dReads;
            l_dFastaReferenceFile.close();
            l_dOutputFile.close();
        }
    }

    return (EXIT_SUCCESS);
}
