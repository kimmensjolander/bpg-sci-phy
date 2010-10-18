#include "OutputFile.h"
#include "../general/utils.h"

using namespace std;

namespace sciphy
{

  OutputFile::OutputFile()
  {
  }

  OutputFile::~OutputFile()
  {
    // If it is open, close it and say that a file has been created!!!!!
    if(file.get()) {
      file->close();
      utilityObject->info("file created:   [%s]\n", name.c_str());                         
    }
  }

  bool
  OutputFile::openFile(string& fileName, string& fileType)
  {
    if (fileName.empty())
    {
      return false;
    }

    file.reset(new ofstream(fileName.c_str(), ofstream::trunc));
                
    if(!file->is_open()) 
    {
      utilityObject->error("Cannot open output file [%s]\n", fileName.c_str()); 
      return false;
    }

    name = fileName; 
    typeOfFile = fileType;
    
    return true;
  }

  bool
  OutputFile::isOpen()
  {
    return file->is_open();
  }

  ofstream*
  OutputFile::getPtrToFile()
  {
    return file.get();
  }

}
