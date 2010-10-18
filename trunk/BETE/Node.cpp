#include "Node.h"
#include <iostream>
#include "math.h"

using namespace std;

namespace sciphy
{
  // construct a node from a single sequence
  Node::Node(int nodeId, int seqIndex)
  {
    id = nodeId;
    numSeqs = 1;
    validFlag = true;
    seqWeights = new float[numSeqs];

    _seqIndexes.push_back(seqIndex);
    seqWeights[0] = 1;

    left = right = 0;
    prev = next = NULL; 
    
    int i, j, res;
    vector<int>* tempSeq;

    numCols = alignPtr->getLength();
    skipFlags = new bool[numCols];

    tempSeq = alignPtr->getSequence(seqIndex);
    vector<int> mySeq;

    if (userParameters->skipInserts()) {
      for (j=0; j<numCols; j++) {
	skipFlags[j] = false;
	res = tempSeq->at(j);
	// skip lower case residues (representing inserts in UCSC a2m alignment format. Note gap is represented by ., instead of - in other regular columns)
	if (res > 24) {
	  skipFlags[j] = true;
	  continue;
	}
	mySeq.push_back(res);
       }
      numCols = mySeq.size();
    } else {
      for (j=0; j<numCols; j++) {
	skipFlags[j] = false;
	res = tempSeq->at(j);
	// convert lower case to upper case
	if (res > 24) 
	  res = res - 24;
	mySeq.push_back(res);
      }
    }

    // get the raw counts of each residue type in every column

    columns = new Column[numCols];
    Column* col;

    for (j=0; j<numCols; j++) {
      col = &columns[j];

      for (i=0; i<ALPHABET_SIZE; i++)
	col->counts[i] = 0;

      res = mySeq.at(j); // encoded residue

      if (res >= GAP) {
	col->residueCount = 0;
	for (i=0; i<ALPHABET_SIZE; i++) 
	  col->probs[i] = DirichletPriors[i];
      } else {
	col->counts[res] = 1;
	col->residueCount = 1;
	for (i=0; i<ALPHABET_SIZE; i++) 
	  col->probs[i] = SingleResidueProbMatrix.getAt(res,i);
      }
    }

   }

  // join 2 nodes
  Node::Node(int nodeId, Node* node1, Node* node2)
  {
    id = nodeId;
    numSeqs = node1->getNumSeqs() + node2->getNumSeqs();
    numCols = node1->getNumCols();
    validFlag = true;
    left = node1;
    right = node2;

    int i, j;
    vector<int>* p = node1->getSeqIndexes();
    for (i=0; i<node1->getNumSeqs(); i++)
      _seqIndexes.push_back(p->at(i));

    p = node2->getSeqIndexes();
    for (i=0; i<node2->getNumSeqs(); i++)
      _seqIndexes.push_back(p->at(i));

    int alnLen = alignPtr->getLength();
    skipFlags = new bool[alnLen];
    for (j=0; j<alnLen; j++)
      skipFlags[j] = node1->getSkipFlag(j);

    columns = new Column[numCols];
    Column* col;

    /* add up the counts */
    Column c1, c2;
    for (j=0; j<numCols; j++) {
      c1 = node1->getCol(j);
      c2 = node2->getCol(j);
      col = &columns[j];
      
      for (i=0; i<ALPHABET_SIZE; i++) 
	col->counts[i] = c1.counts[i] + c2.counts[i];

      col->residueCount = c1.residueCount + c2.residueCount;
    }

    createProfile();
  }
 
  Node::~Node()
  {
    _seqIndexes.clear();
    delete [] columns;
    delete [] seqWeights;
    delete [] skipFlags;
  }

  void
  Node::removeNode(Node* node)
  {
    if (node->prev == NULL) { 
      node->next->prev = NULL;
    } else if (node->next == NULL) {
      node->prev->next = NULL;
    } else {
      node->prev->next = node->next;
      node->next->prev = node->prev;
    }
  }

  void
  Node::clearMemory() {
    delete [] columns;
    delete [] seqWeights;
    delete [] skipFlags;
  }

  void
  Node::createProfile() {

    seqWeights = new float[numSeqs];
    float weightedCounts[ALPHABET_SIZE];
    double tempMix[ALPHABET_SIZE];
    Column* col;
    int i, j;

    // multiple sequences, calculate the posterior probability
    // using weighted counts

    calcTotalWeight();

    string rwMethod = userParameters->getRelativeWeight();

    if (rwMethod == "none") {
      float weight = totalWeight/numSeqs;

      for (j=0; j<numCols; j++) {
	col = &columns[j];
	
	if (col->residueCount == 0) {
	  // there are only gaps in this column in the subtree, use background probabilities
	  for (i=0; i<ALPHABET_SIZE; i++) 
	    col->probs[i] = DirichletPriors[i];
	} else {
	  for (i=0; i<ALPHABET_SIZE; i++)
	    weightedCounts[i] = col->counts[i]*weight;

	  dirichlet->calcPosteriorProb(weightedCounts, col->probs);
	}
      }

    } else if (rwMethod == "pw") {
      for (i=0; i<numSeqs; i++)
	seqWeights[i] = 0; // init
      
      calcRelativeWeights();

      vector<int> alignColumn;
      int idx, res;

      j = 0;

      for (int pos=0; pos < alignPtr->getLength(); pos++) {
	if (skipFlags[pos]) 
	  continue;

       	col = &columns[j];

	if (col->residueCount == 0) {
	  // there are only gaps in this column in the subtree, use background probabilities
	  for (i=0; i<ALPHABET_SIZE; i++) 
	    col->probs[i] = DirichletPriors[i];
	} else {
	  for (i=0; i<ALPHABET_SIZE; i++)
	    weightedCounts[i] = 0;
      
	  for (i=0; i<numSeqs; i++) {
	    idx = _seqIndexes[i];
	    res = alignPtr->residue(idx, pos);
	    if (res > GAP) continue;
	    weightedCounts[res] += seqWeights[i];
	  }

	  dirichlet->calcPosteriorProb(weightedCounts, col->probs);
	}

	j++;
      }

    }

    //cout << "Finished creating profile" << endl;
  }

  void
  Node::calcTotalWeight()
  {
    int i, j;

    /* calculate the total weight */
    int mostFreqCount;
    Column col;

    float totp = 0;
    int nc = 0; // # of valid columns, I don't think we should use numCols [HH]

    for (j=0; j<numCols; j++) {
      mostFreqCount = 0; 

      col = columns[j];
      if (col.residueCount==0)
	continue;

      for (i=0; i<ALPHABET_SIZE; i++) {
	if (col.counts[i] > mostFreqCount)
	  mostFreqCount = col.counts[i];
      }
      
      // calculate the frequency of the most frequent amino acid
      totp += mostFreqCount/col.residueCount;
      nc++;
    }

    float pcons = totp/nc;
    float NIC = pow(numSeqs, 1-pcons); // NIC = N^(1-Pcons)
    //cout << "NIC=" << NIC << endl;

    totalWeight = NIC;
  }

  void
  Node::calcRelativeWeights()
  {
    int i, j, idx, res;

    float totw = 0, resw;
    int residueTypeCount;

    vector<int> alignColumn;
    Column col;

    j = 0;
    for (int pos=0; pos < alignPtr->getLength(); pos++) {
      if (skipFlags[pos])
	continue;

      residueTypeCount = 0;

      col = columns[j];
      j++;

      for (i=0; i<ALPHABET_SIZE; i++) {
	if (col.counts[i] > 0)
	  residueTypeCount += 1;
      }

      if (residueTypeCount == 0)
	continue;

      for (i=0; i<numSeqs; i++) {
	idx = _seqIndexes[i];
	res = alignPtr->residue(idx, pos);;

	if (res >= GAP)
	  continue; // ignore gaps or other non-standard AA
	// a consequence of this is that longer seq will get more weight [HH]

	resw = (float) 1/(residueTypeCount*col.counts[res]);
	seqWeights[i] += resw;
	totw += resw;
      }

    }

    /* Normalize the weights */
    //cout << "Finished calculating relative sequence weights." << endl;
    float factor = totalWeight/totw;
    for (i=0; i<numSeqs; i++) {
      seqWeights[i] *= factor;
      // cout << i << ":" << seqWeights[i] << endl;
    }
    
  }

  void
  Node::getLeafSeqIndexes(vector<int>* indexes) 
  {
    if (isLeafNode()) {
      indexes->push_back(_seqIndexes.at(0));
    } else {
      left->getLeafSeqIndexes(indexes);
      right->getLeafSeqIndexes(indexes);
    }
  }

  void
  Node::printNewickTree(ofstream* ptrToFile)
  {
    int idx;
    if (isLeafNode()) {
      idx = _seqIndexes.at(0);
      string name = alignPtr->getName(idx);
      (*ptrToFile) << name;
    } else {
      (*ptrToFile) << "(" << endl;
      left->printNewickTree(ptrToFile);
      (*ptrToFile) << "," << endl;
      right->printNewickTree(ptrToFile);
      (*ptrToFile) << ")" << endl;
    }
  }

}
