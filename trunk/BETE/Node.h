#ifndef NODE_H
#define NODE_H

#include <vector>
#include "../general/Array2D.h"
#include "../general/param.h"
#include "../general/userparams.h"
#include "../alignment/Alignment.h"
#include "../general/dirichlet.h"
#include "../fileInput/OutputFile.h"

using namespace std;

namespace sciphy
{

  extern float* DirichletPriors;
  extern Array2D<double> SingleResidueProbMatrix;
  extern Alignment* alignPtr;

  struct Column
  {
    int pos;
    int residueCount;
    int counts[ALPHABET_SIZE];
    double probs[ALPHABET_SIZE];
  };

  class Node
  {
  public:
    Node(int nodeId, int seqIndex);
    Node(int nodeId, Node* node1, Node* node2);
    ~Node();

    bool isLeafNode();
    bool isValid() { return validFlag; }
    void setValid(bool flag) { validFlag = flag; }
    double getEncodingCost();

    int getId() { return id; }
    int getNumSeqs() { return numSeqs; };
    int getNumCols() { return numCols; };
    vector<int>* getSeqIndexes() { return &_seqIndexes; };
    void removeNode(Node* node);
    void getLeafSeqIndexes(vector<int>* indexes); 
    Column& getCol(int i) { return columns[i]; }
    bool getSkipFlag(int i) {return skipFlags[i]; }
    void clearMemory();
    void printNewickTree(ofstream* ptrToFile);
   
    /* Attributes */
    Node *prev;  
    Node *next; // for a linked list
    Node *left; // for tree structure
    Node *right; // a leaf node will have both left=right=0


  private:
    int id; // a unique identifier for the node
    int numSeqs; // number of sequences in this node
    int numCols; // profile width, not necessarily the alignment length (if skip_inserts is enabled) 
    bool validFlag; // indicate if the node is valid
    float* seqWeights; // vector of sequence weights
    float totalWeight; // total weight (NIC is currently used)
    double encodingCost;

    Column* columns;
    bool* skipFlags; // vector of alignment length, indicate if the position is an insert state and should be skipped

    // index of the sequence in the input MSA, so we don't have the store the actual sequence info
    vector<int> _seqIndexes;

    void createProfile();
    void calcTotalWeight();
    void calcRelativeWeights();
  };

  inline
  bool Node::isLeafNode()
  {
    if (left==0 && right==0) {
      return true;
    } else {
      return false;
    }
  }
  
}

#endif





