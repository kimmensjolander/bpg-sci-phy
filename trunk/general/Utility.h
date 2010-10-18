#ifndef UTILITY_H
#define UTILITY_H

#include <stdio.h>
#include <stdarg.h>

namespace sciphy
{
  class Utility
  {
  public:
    Utility();
    
    void error( char *msg,...);
    void warning( char *msg,...);
    void info( char *msg,...);
    void beQuiet(bool b) { quiet=b; }

  private:
    bool quiet;
  };

}

#endif
