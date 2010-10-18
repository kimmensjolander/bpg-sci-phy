#ifndef AGGLOMERATIONSTEP_H
#define AGGLOMERATIONSTEP_H

#include <vector>
#include "Node.h"

using namespace std;

namespace sciphy
{

  class AgglomerationStep
  {
  public:
    AgglomerationStep();
    AgglomerationStep(int nodeId, Node* nodeList, double cost);
    ~AgglomerationStep();

    int getStep() { return step; }
    double getEncodingCost() { return encodingCost; }
    vector< vector<int> >*  getGroups() { return &groups; }

  private:
    int step;
    double encodingCost;
    vector< vector<int> > groups; // groups for the current tree forest
  };

}

#endif
