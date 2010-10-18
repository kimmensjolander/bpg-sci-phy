#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include <vector>
#include <string>
#include <iomanip>
#include "Sequence.h"
#include "../general/Array2D.h"

using namespace std;

namespace sciphy
{
  typedef std::vector<vector <int> > SeqArray;

  class Alignment
  {
  public:
    Alignment();
    ~Alignment();

    void addSequenceObj(Sequence* seq);
    void clearSeqArray();

    vector<int>* getSequence(int index) {return &seqArray[index];}
    int residue(int i, int j) { return seqArray[i][j]; }
    void writeSequence(int i, ofstream* ptrToFile);

    string getName(int index);
    string getTitle(int index);
    int getLength() { return alnLength;}
    int getNumSeqs() { return numSeqs; }
    void calcConservation();
    bool isConserved(int index);

  private:

    int numSeqs;
    int alnLength;
    SeqArray seqArray;
    vector<string> names;
    vector<string> titles;
    vector<bool> conservedFlags;

    //bitset conservedFlags;
  };
  
}

#endif

