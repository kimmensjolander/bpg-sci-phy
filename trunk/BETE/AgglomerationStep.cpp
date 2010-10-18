#include "AgglomerationStep.h"

#include <iostream>

using namespace std;

namespace sciphy
{
  AgglomerationStep::AgglomerationStep()
  {
  }

  AgglomerationStep::AgglomerationStep(int stepNum, Node* nodeList, double cost)
  {
    step = stepNum;
    encodingCost = cost;

    for (Node* n=nodeList; n; n=n->next) {
      vector<int> x;
      n->getLeafSeqIndexes(&x);
      groups.push_back(x);
    }

  }

  AgglomerationStep::~AgglomerationStep()
  {
    for (int i=0; i<groups.size(); i++) {
      groups[i].clear();
    }
    groups.clear();
  }

}
