#ifndef USERPARAMETERS_H
#define USERPARAMETERS_H

#include <vector>
#include <string>
#include <iostream>
#include "sciphy_version.h"

using namespace std;

namespace sciphy
{

class UserParameters
{
  
public:
  UserParameters();
  ~UserParameters();

  int parseParams(vector<string>* args);

  int resIndex(char c);
  char resChar(int i) { return aminoAcidCodes[i]; }
  string getAminoAcids() { return aminoAcidCodes; }

  void printHelp();

  string getInfileName() { return infileName; }
  string getOutfileName() { return outfileName; }
  string getMixtureName() { return mixtureName; }
  string getTreeFileName() { return treeFileName; }

  string getTotalWeight() { return totalWeight; }

  void setRelativeWeight(string val) { relativeWeight = val; }
  string getRelativeWeight() { return relativeWeight; }

  void setMinOverlap(int val) { minOverlap = val; }
  int getMinOverlap() { return minOverlap; }

  bool skipInserts() { return bSkipInserts; }

private:
  string aminoAcidCodes;
  string mixtureName;
  string relativeWeight;
  string totalWeight;
  string infileName;
  string outfileName;
  string treeFileName;
  int   minOverlap;
  bool   bSkipInserts;
  bool isInteger(char* number);
};

}

#endif
