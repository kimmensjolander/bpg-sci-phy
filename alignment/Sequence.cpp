#include "Sequence.h"
#include <iostream>

using namespace std;

namespace sciphy
{

  Sequence::Sequence() 
  {
    _name = "";
    _title = "";
  }

  Sequence::~Sequence() {
    
  }

  Sequence::Sequence(string& seq, string& name, string& title)
  {
    copyStringIntoVector(&_sequence, &seq);
    encodeSequence();
    _name = name;
    _title = title;
  }
  
  void
  Sequence::encodeSequence()
  {
    std::vector<char>::iterator it;

    for (it = _sequence.begin(); it!= _sequence.end(); it++) {
      //cout << "residue " << *it << " index " << userParameters->resIndex(*it) << endl;
      _encodedSequence.push_back(userParameters->resIndex(*it));
    }
  }

  void
  Sequence::copyStringIntoVector(vector<char>* _vectorTo, string* _stringFrom)
  {
    _vectorTo->clear();
    for (int i=0; i< (int)_stringFrom->size(); i++) {
      _vectorTo->push_back(_stringFrom->at(i));
    }
  }

    
  std::vector<int>* Sequence::getSequence()
  {
     return &_encodedSequence;
   }

  std::string Sequence::getName()
  {
     return _name;
  }

  std::string Sequence::getTitle()
  {
    return _title;
  }
  
}
