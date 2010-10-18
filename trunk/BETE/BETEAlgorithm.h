#ifndef BETEALGORITHM_H
#define BETEALGORITHM_H

#include <vector>
#include <queue>
#include "Node.h"
#include "AgglomerationStep.h"

namespace sciphy
{
  extern Alignment* alignPtr;

  struct NodePair
  {
    Node* node1;
    Node* node2;
    float dist;
  };

  // sort queue such that pair with lower dist will rank on top
  class CompareDist {
    public:
    bool operator()(NodePair& n1, NodePair& n2) // Returns true if n1 has larger distance
    {
      if (n1.dist > n2.dist)
	return true;
      else
	return false;
    }
  };

  typedef priority_queue<NodePair, vector<NodePair>, CompareDist> pqueue;

  class BETEAlgorithm
  {

  public:
    BETEAlgorithm();
    ~BETEAlgorithm();
    
    void run();
    
    int getNodeCount() { return nodeCount; }
    int getCurrentStep() { return currentStep; }
    int getNumClasses() { return numClasses; }

  private:
    Node* nodeList; // a linked list used for agglomeration
    int currentStep;
    int numSeqs;
    int numClasses;
    int nodeCount;
    vector<bool> validNodeFlags;
    vector<AgglomerationStep> agglomerationSteps;
    pqueue nodeQueue;

    void init();
    void initializeNodes();
    float calcPairDistance(Node* node1, Node* node2);
    double calcEncodingCost();
  };

}

#endif
