/*
 * Stack.cpp
 *
 *  Created on: Apr 22, 2011
 *      Author: Douglas W Bryant Jr
 */

#include "Stack.h"

Stack::Stack(int argc, char** argv) {
	m_dParams = new StackParameters();
	if(argc < 9) { Stack::outputUsage(); return; }
	else {
		char* parameter;
		for(int i = 1; i < argc; i++) {
			parameter = argv[i];
			if(parameter[0] == '-') {
				switch(parameter[1]) {
					case 'r': m_dParams->setReferenceFileName(argv[++i]); break;
					case 's': m_dParams->setSplatOutputFileName(argv[++i]); break;
					case 'c': m_dParams->setMinReadCopyNumber(atoi(argv[++i])); break;
					case 'd': m_dParams->setMaxReadCopyNumber(atoi(argv[++i])); break;
					case 'n': m_dParams->setMinNumDiffSequences(atoi(argv[++i])); break;
					case 'o': m_dParams->setStackOutputFileName(argv[++i]); break;
					case 'R': m_dParams->setReferenceNameToAnalyze(argv[++i]); break;
					default: std::cout << "Error in " << argv[i] << " " << argv[i+1] << std::endl << std::endl;
				}
			}
		}
	}
	Stack::run();
}

Stack::~Stack() {
	delete m_dParams;
}

bool Stack::run() {
	bool l_bSuccess = true;
	if(!m_dParams->verifyParameters()) {
		Stack::outputUsage();
		l_bSuccess = false;
	} else {
		Stacker* stacker = new Stacker(m_dParams);
		bool success = stacker->initReferenceInfos(*(m_dParams->getReferenceFileName()));
		success = stacker->stackSplatOutput(*(m_dParams->getSplatOutputFileName()));
		printResults(stacker->getStackedResults());
		delete stacker;
	}

	return l_bSuccess;
}

void Stack::printResults(std::vector<std::string>* p_vResultsVector) {
	if(!p_vResultsVector->empty()) {
		std::ofstream l_dOutputFile;
		std::string l_sFileName = *(m_dParams->getStackOutputFileName());
		l_dOutputFile.open(l_sFileName.c_str(), std::ios::app | std::ios::out);
		for(int i = 0; i < p_vResultsVector->size(); i++) { l_dOutputFile << p_vResultsVector->at(i); }
		l_dOutputFile.close();
	}
}

void Stack::outputUsage() {
	std::cout << std::endl;
    std::cout << "Usage:   supersplat stack <options>" << std::endl;
    std::cout << "Options: -r <Fasta-formatted reference file.> (REQD)" << std::endl;
    std::cout << "         -s <Splat output file> (REQD)" << std::endl;
    std::cout << "         -n <Minimum number of unique hits.> (REQD)" << std::endl;
    std::cout << "         -c <Minimum read copy number per unique hit.> (REQD)" << std::endl;
    std::cout << "         -d <Maximum read copy number per unique hit.> (Optional, default = no such maximum)" << std::endl;
    std::cout << "         -o <Output file name.> (Optional, default = name_of_splat_output.stack)" << std::endl;
    std::cout << "Additional options:" << std::endl;
    std::cout << "         -R <Name of reference over which to restrict analysis> (Optional, default = no such restriction)" << std::endl;
}
