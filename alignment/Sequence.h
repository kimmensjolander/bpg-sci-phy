#ifndef SEQUENCE_H
#define SEQUENCE_H

#include <vector>
#include <string>
#include "../general/userparams.h"

using namespace std;

namespace sciphy
{
  class Sequence
  {
  public:
    Sequence();
    Sequence(string& seq, string& name, string& title);
    ~Sequence();

    void encodeSequence();
    void printSequence();
    vector<int>* getSequence();
    string getName();
    string getTitle();

  private:

    void copyStringIntoVector(vector<char>* _vectorTo, string* _stringFrom);
 
    vector<char> _sequence;
    vector<int> _encodedSequence;
    string _name;
    string _title;
  };

}

#endif

