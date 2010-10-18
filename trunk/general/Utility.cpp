#include "Utility.h"

namespace sciphy
{

  Utility::Utility()
  {
    quiet=false;
  }

  /**
   * @param msg 
   */
  void
  Utility::error(char *msg,...)
  {
    va_list ap;
    
    va_start(ap, msg);
    fprintf(stderr, "\n\nERROR: ");
    vfprintf(stderr, msg, ap);
    fprintf(stderr, "\n\n");
    va_end(ap);
  }

  /**
   * 
   * @param msg 
   */
  void
  Utility::warning( char *msg,...)
  {
    va_list ap;
    
    va_start(ap, msg);
    fprintf(stderr, "\n\nWARNING: ");
    vfprintf(stderr, msg, ap);
    fprintf(stderr, "\n\n");
    va_end(ap);
  }
  
  /**
   * 
   * @param msg 
   */
  void Utility::info( char *msg,...)
  {
    va_list ap;
    va_start(ap, msg);
    
    if(! quiet) {
      fprintf(stdout, "\n");
      vfprintf(stdout, msg, ap);
      va_end(ap);
    }
  }
  
}
