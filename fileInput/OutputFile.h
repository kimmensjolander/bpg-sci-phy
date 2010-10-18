#ifndef OUTPUTFILE_H
#define OUTPUTFILE_H
#include <memory>
#include <fstream>

using namespace std;

namespace sciphy 
{
  
  class OutputFile
  {
  public:
    OutputFile();
    ~OutputFile();
    bool openFile(string& fileName, string& fileType);
    bool isOpen();
    ofstream* getPtrToFile();

  private:
    auto_ptr<ofstream> file;
    string name;
    string typeOfFile;
  };

}

#endif
