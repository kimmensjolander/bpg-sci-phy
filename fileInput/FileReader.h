#ifndef FILEREADER_H
#define FILEREADER_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include "../alignment/Alignment.h"
#include "../alignment/Sequence.h"

using namespace std;

namespace sciphy
{
  
  class FileReader
  {
  public:
    FileReader(char* filename);
    ~FileReader();
    int readAlignment(Alignment *alignPtr);

  private:

    string filename;
    ifstream* fileIn;
  };

}

#endif
