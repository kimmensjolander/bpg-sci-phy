#include <iostream>
#include "general/param.h"
#include "fileInput/FileReader.h"
#include "alignment/Alignment.h"
#include "general/Utility.h"
#include "general/UserParameters.h"
#include "general/DirichletMixture.h"
#include "BETE/BETEAlgorithm.h"

namespace sciphy
{
  Alignment* alignPtr;
  UserParameters* userParameters;
  Utility* utilityObject;
  //extern int n_reuse;
  //extern int n_new;
}

using namespace std;
using namespace sciphy;

int main(int argc, char **argv) 
{

  userParameters = new UserParameters();
  alignPtr = new Alignment();
  utilityObject = new Utility();

  if (argc > 1) {    
    time_t start, end;
    double dif;
    start = time (NULL);

    vector<string> args;
    for (int i = 1; i < argc; ++i)
      args.push_back(argv[i]);
    
    int result = userParameters->parseParams(&args);
    if (!result) {
      exit(1);
    }

    string infileName = userParameters->getInfileName();
    
    FileReader* fileReader = new FileReader( (char*) infileName.c_str()); 
    fileReader->readAlignment(alignPtr);
    
    alignPtr->calcConservation();

    BETEAlgorithm* bete = new BETEAlgorithm();
    bete->run();

    end = time (NULL);
    dif = difftime(end, start);
    cout << "\nRunning time: " << dif << " seconds\n";

    //cout << "# of new " << n_new << endl;
    //cout << "# of reuse " << n_reuse << endl;
  } else {
    userParameters->printHelp();
  }

  return 0;

}
