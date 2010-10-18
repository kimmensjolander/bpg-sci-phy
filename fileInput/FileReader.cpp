#include "FileReader.h"
#include <fstream>
#include <iostream>

using namespace std;

namespace sciphy
{

  FileReader::FileReader(char* filename) 
  {
    fileIn = new ifstream;

    fileIn->open(filename);
    if (!fileIn) {
      cerr << "Error: Unable to open input file " << filename << endl;
      exit(1);
    } else {
      cout << "Opening file " << filename << " to read" << endl;
    }
  }

  FileReader::~FileReader()
  {
    delete fileIn;
  }

  int
  FileReader::readAlignment(Alignment* alignPtr)
  {
    string line;
    string seq = "";
    string name = "";
    string title = "";
    int lineNum = 0;

    Sequence tempSeq;

    while(! fileIn->eof() ) {
      getline(*fileIn, line);
      lineNum++;

      //cout << "Read: " << line << endl;

      // we've found the def line
      if (line[0] == '>') {
	if (seq != "") {
          //cout << seq << endl;
	  tempSeq = Sequence::Sequence(seq, name, title);
	  alignPtr->addSequenceObj(&tempSeq);
	  seq = name = title = "";
	}

	int found = line.find_first_of(' ');
	name = line.substr(1,found-1);
	title = line.substr(found+1, line.length()-name.length()-1);
	//cout << "sequence found: " << "name=" << name << " title=" << title << endl;
      } else if (lineNum==1) {
	cerr << "Error: The input file is not in FASTA format." << endl;
	exit(1);
      } else {
	seq.append(line);
      }
    }

    // process the last seq
    tempSeq = Sequence::Sequence(seq, name, title);
    alignPtr->addSequenceObj(&tempSeq);

    return 1;
  }

}
