#ifndef ARRAY2D_H
#define ARRAY2D_H

#include <vector>

namespace sciphy
{

using namespace std;

  template<class T>
  class Array2D
  {
  public:
    Array2D(): dimRow(0), dimCol(0) {;}

    /* This dynamic 2D array is implemented column-wise, useful for accessing columns in a profile */
    Array2D(int nRow, int nCol)
    {
      dimRow = nRow;
      dimCol = nCol;

      for (int i=0; i<nCol; i++) {
	vector<T> x(nRow);
	int y = x.size();
	elements.push_back(x);
      }
    } 

    ~Array2D() 
    {
      clearArray();
    }

    void resizeRect(int nRow, int nCol)
    {
      clearArray();

      dimRow = nRow;
      dimCol = nCol;

      for (int i=0; i<nCol; i++) {
	vector<T> x(nRow);
	int y = x.size();
	elements.push_back(x);
      }
    } 
 
    T getAt(int iRow, int iCol)
    {
      return elements[iCol][iRow];
    }

    void
    setCol(vector<T>& col, int iRow, const T& val) 
    {
       col[iRow] = val;
    }

    void setAt(int iRow, int iCol, const T& val)
    {
      elements[iCol][iRow] = val;
    }

    void setAllArrayToVal(const T& val)
    {
      for (int i=0; i<dimCol; i++) {
        vector<T>& tmpCol = elements[i];
	for (int j=0; j<dimRow; j++) {
	  tmpCol[j] = val;
	}
      }
    }

    void clearArray()
    {
      int size = elements.size();
      for (int i=0; i<size; i++) {
	elements[i].clear();
      }
      elements.clear();
      dimRow = 0;
      dimCol = 0;
    }

    int getRowSize() { return dimRow; }
    int getColSize() { return dimCol; }
 
    vector<T>& operator[](int x)
    {
      return elements[x];
    }

  private:
    unsigned int dimRow;
    unsigned int dimCol;
    vector< vector<T> > elements;

  };

}

#endif
