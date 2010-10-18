#include <sstream>
#include <iostream>
#include "Alignment.h"
#include "../general/param.h"
#include "../general/userparams.h"
#include <fstream>

using namespace std;

namespace sciphy
{
  Alignment::Alignment()
  {
    numSeqs = 0;
    alnLength = 0;
  }
 
  Alignment::~Alignment()
  {
    clearSeqArray();
  }

  void
  Alignment::addSequenceObj(Sequence* seq) 
  {
    numSeqs += 1;
    vector<int>* seqcode = seq->getSequence();
    if (alnLength==0)
      alnLength = seqcode->size();

    seqArray.push_back(*seqcode);
    names.push_back((*seq).getName());
    titles.push_back((*seq).getTitle());
  }

  void
  Alignment::clearSeqArray()
  {
    for (int i=0; i<(int)seqArray.size(); i++) {
      seqArray[i].clear();
    }
    seqArray.clear();
  }

  string
  Alignment::getName(int index)
  {
    return names[index];
  }

  string
  Alignment::getTitle(int index)
  {
    return titles[index];
  }

  bool
  Alignment::isConserved(int index)
  {
    return conservedFlags[index];
  }

  void
  Alignment::calcConservation() 
  {
    for (int j=0; j<alnLength; j++) {
      bool conserved = true;

      int firstAA = seqArray[0][j]; //[seq][col]

      for (int i=1; i<numSeqs; i++) {
	int AA = seqArray[i][j];

	if (AA == GAP) // gap character
	  continue;

	if (firstAA == GAP) {
	  firstAA = AA;
	  continue;
	}

	if (AA != firstAA) {
	  conserved = false;
	  break;
	}
      }
      conservedFlags.push_back(conserved);
    }
  }

  void
  Alignment::writeSequence(int i, ofstream* ptrToFile) {
    vector<int> seq = seqArray[i];

    int width = 80;
    int j, count=0, resi;
    char aa;
    (*ptrToFile) << ">" << names[i] << " " << titles[i] << endl;
    for (int j=0; j < seq.size(); j++) {
      resi = seq[j];
      aa = userParameters->resChar(resi);
      (*ptrToFile) << aa;
      count++;
      if (count == width) {
	(*ptrToFile) << endl;
	count = 0;
      }
    }

    int remainder = seq.size() / width;
    if (remainder != 0) 
      (*ptrToFile) << endl;
  }

}
