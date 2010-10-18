#include "UserParameters.h"

namespace sciphy
{
  UserParameters::UserParameters() 
  {
    // set all default params here
    // . is the gap character for insert state in SAM alignment
    aminoAcidCodes = "ACDEFGHIKLMNPQRSTVWY-BUXZacdefghiklmnpqrstvwy.";
    mixtureName = "recode4";
    minOverlap = 10;
    relativeWeight = "pw";
    totalWeight = "NIC";
    infileName = "";
    outfileName = "";
    treeFileName = "";
    bSkipInserts = false;
  }

  UserParameters::~UserParameters()
  {
  }

  // returns the index of the chracter c in string aminoAcidCodes
  int
  UserParameters::resIndex(char c) 
  {
    register int i;

    for (i = 0; aminoAcidCodes[i] && aminoAcidCodes[i]!=c; i++)
      ;

    if (aminoAcidCodes[i]) {
      return (i);
    } else {
      return -1;
    }
  }

  int
  UserParameters::parseParams(vector<string>* args)
  {
    string str;

    for (int i=0; i<args->size(); i++) {
      str = args->at(i);

      if (str == "-h") {
	printHelp();
	exit(0);
      } else if (str == "-i") {
	i++;
	str = args->at(i);
	infileName = str;
      } else if (str == "-o") {
	i++;
	str = args->at(i);
	outfileName = str;
      } else if (str == "-m") {
        i++;
        str = args->at(i);
        mixtureName = str;
	if (mixtureName != "recode4" && mixtureName != "blocks9") {
	  cerr << "Error: Invalid mixture. Choose between recode4 and blocks9." << endl;
	  return 0;
	}
      } else if (str == "-rw") {
	i++;
	str = args->at(i);
	setRelativeWeight(str);
        if (0) { // add code to check options
	  cerr << "Error: Invalid argument for -rw option for relative weights. Valid options are none or pw (position-weighted)." << endl;
	  return 0;
	}
      } else if (str == "-mino") {
	i++;
	str = args->at(i);
	char* val = (char*) str.c_str();
	if (isInteger(val) && (int)val <=100 && (int)val >=0) {
	  setMinOverlap((int) val);
	} else {
	  cerr << "Error: Invalid number for -mino option. The minimum overlap has to be a number between 0 and 100." << endl;
	  return 0;
	}
      } else if (str == "-tree") {
	int found = infileName.find_first_of(".");
	if (found) {
	  treeFileName = infileName.substr(0, found);
	} else {
	  treeFileName = infileName;
	}
	treeFileName.append(".ph");
      } else if (str == "-skip_inserts") {
	bSkipInserts = true;
      } else {
	cerr << "Error: Unknown option " << str << endl;
	return 0;
      }
    }

    if (infileName == "") {
      cerr << "Error: No input file is specified." << endl;
      return 0;
    } else if (outfileName == "") {
      cerr << "Error: No output file is specified." << endl;
      return 0;
    }

    return 1;
  }

  void 
  UserParameters::printHelp() {
    string title = "SCIPHY - Subfamily Classification In PHYlogenomics";
    string reference = "PLoS Computational Biology, 3, e160, 1526-1538 (2007)";
    string usage = 
      "Usage: sciphy [options] -i aln_file -o out_file\n\n"
      "Available options are:\n"
      "   -h               : print help\n"
      "   -m               : dirichlet mixture: recode4, blocks9 (Default=recode4)\n"
      "   -rw              : relative weighting method: pw, none. pw is Henikoff position-weighted scheme. (Default=pw)\n"
      "   -skip_inserts    : skip columns representing insert states, in which residues are represented in lower case in alignment in UCSC a2m format.\n"
      "   -mino            : minimum percent overlap between two profiles for them to be joined. (Default=10)\n"
      "   -tree            : output tree (Newick format) in a file\n"
      ;

    cout << "------------------------------------------------------------------------------\n";
    cout << title << endl;
    cout << "VERSION: " << SCIPHY_VERSION << endl;
    cout << reference << endl;
    cout << "------------------------------------------------------------------------------\n";
    cout << usage << endl;
   }

  bool
  UserParameters::isInteger(char *number )
  {
    int len = strlen( number );
    bool isnumber = true;
    int i = 0;

    if( len == 0 || len > 9 ) // if the user did not enter an input
      {                         // or if the input is too big to be
	return false;         // a "long" value.
      }

    while( i < len && isnumber ) // scanning the the input
      {
	if( isdigit( number[i] ) == 0 ) // check if there is a none digit character in the input
	  {
            if( number[i] == '+' || number[i] == '-' ) // if so we verify if it is "+" or "-"
	      {
                if(  i + 1 > len - 1 || i - 1 >= 0 ) // check the position of "+" or "-"
		  {                                    // this signs could only be infront of
		    isnumber = false;                 // the number,no where else.
		  }
	      }
            if( number[i] != '+' && number[i] != '-' ) // if it's not "+" or "-" than
	      {                                          // the expression is not a number
	      isnumber = false;
	      }
	  }
	i++;
      }
    
    return isnumber;
  }

}
