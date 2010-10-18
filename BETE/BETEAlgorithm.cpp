#include "BETEAlgorithm.h"
#include "../general/userparams.h"
#include <cmath>
#include <iostream>
#include "../fileInput/OutputFile.h"

using namespace std;

namespace sciphy {

  DirichletMixture* dirichlet;
  float* DirichletPriors;
  Array2D<double> SingleResidueProbMatrix(ALPHABET_SIZE, ALPHABET_SIZE);
  
  BETEAlgorithm::BETEAlgorithm()
  {
    currentStep = 0;
    nodeCount = 0;
    nodeList = NULL;
  }

  BETEAlgorithm::~BETEAlgorithm()
  {
    Node* curr = nodeList, *temp;
    while (curr != NULL) {
      temp = curr;
      curr = curr->next;
      delete(temp);
    }

    delete DirichletPriors;
    delete dirichlet;
    SingleResidueProbMatrix.clearArray();
  }


  /*
   * This is the main program, which does all the initializations and agglomeration
   */
  void
  BETEAlgorithm::run() 
  {
    init();
    initializeNodes();
    currentStep++;

    double encodingCost = calcEncodingCost();
    cout << "Step " << currentStep << ": EncodingCost = " << encodingCost << endl;
    
    AgglomerationStep step(currentStep, nodeList, encodingCost);
    agglomerationSteps.push_back(step);

 
    NodePair bestpair;

    float dist;
    // calculate upper triangle of the distance matrix
    for (Node* node1=nodeList; node1; node1=node1->next) {
      for (Node* node2=node1->next; node2; node2=node2->next) {

	dist = calcPairDistance(node1, node2);
	//cout << node1->getId() << " : " << node2->getId() << " dist=" << dist << endl;
	  
	NodePair p;
	p.node1 = node1;
	p.node2 = node2;
	p.dist = dist;

	// push the pair with their distance into the priority queue
	// note that pairs will be stored in the queue according to the dist value
	// as a priority queue is a heap sort
	nodeQueue.push(p);
      }
    }

    while (numClasses > 1) {
      currentStep++;

    GET:
      bestpair = nodeQueue.top();
      //cout << "pair: "<< bestpair.node1->getId() << "," << bestpair.node2->getId() << endl;
      nodeQueue.pop();

      /* 
       * Some of the nodes have been joined, so they are no long valid, but they might still be distance
       * pair. Skip them.
       */
      if (! bestpair.node1->isValid() || ! bestpair.node2->isValid()) {
	goto GET;
      }

      cout << "Found best scoring pair: (" << bestpair.node1->getId() << "," << bestpair.node2->getId() << "). Distance = " << bestpair.dist << endl;

      cout << "Creating a new node: Id=" << nodeCount << endl;
      Node* newnode = new Node(nodeCount, bestpair.node1, bestpair.node2);
      nodeCount++;
      numClasses--;

      // add the new node to the list
      newnode->next = nodeList;
      nodeList->prev = newnode;
      nodeList = newnode;

      // remove 2 nodes that have just been joined from list and mark them as unusable 
      bestpair.node1->setValid(false);
      bestpair.node2->setValid(false);

      // release memory in profile columns
      bestpair.node1->clearMemory();
      bestpair.node2->clearMemory();

      nodeList->removeNode(bestpair.node1);
      nodeList->removeNode(bestpair.node2);

      // Calculating the encoding cost for the new partition
      double encodingCost = calcEncodingCost();

      cout << "\nStep " << currentStep << ": EncodingCost = " << encodingCost << endl;

      AgglomerationStep step(currentStep, nodeList, encodingCost);
      agglomerationSteps.push_back(step);

      if (numClasses == 1) 
	break;

      //cout << "Calculating distance between the new node and all other nodes" << endl;
      for (Node* node2=nodeList; node2; node2=node2->next) {
	if (newnode == node2) 
	  continue;
	
	dist = calcPairDistance(newnode, node2);
	//cout << ": node2 " << node2->getId() << " dist=" << dist << endl;

	NodePair p;
	p.node1 = newnode;
	p.node2 = node2;
	p.dist = dist;

	nodeQueue.push(p);
      }
    } //while

    /*
     * Now let's check the encoding costs for all the partitions
     */

    float minCost = FLT_MAX; 
    AgglomerationStep minStep;
    vector<AgglomerationStep>::iterator p;
    for (p=agglomerationSteps.begin(); p!=agglomerationSteps.end(); p++) {
      if (p->getEncodingCost() < minCost) {
	minCost = p->getEncodingCost();
	minStep = *p;
      }
    }

    /*
     * Print the partition corresponding to the minimum encoding cost tree
     */

    string outfileName = userParameters->getOutfileName();
    OutputFile outfile;
    string outfileFormat = "text";
    outfile.openFile(outfileName, outfileFormat);
    ofstream* ptrToOutfile = outfile.getPtrToFile();

    vector< vector<int> >* groups = minStep.getGroups();
    for (int i=0; i< groups->size(); i++) {
      vector<int> members = groups->at(i);
      (*ptrToOutfile) << "%subfamily " << i+1 << endl;
      for (int j=0; j<members.size(); j++) {
	alignPtr->writeSequence(members[j], ptrToOutfile);
      }
    }

    /*
     * Print the minimum encoding cost tree
     */
    
    string treeFileName = userParameters->getTreeFileName();
    if (treeFileName != "") {
      OutputFile treeFile;
      string fileFormat = "Newick";
      treeFile.openFile(treeFileName, fileFormat);
      ofstream* ptrToTreeFile = treeFile.getPtrToFile();
      cout << "Write Newick tree to file " << treeFileName << endl;
      nodeList->printNewickTree(ptrToTreeFile);
    }
  }

  void
  BETEAlgorithm::init() 
  {
    dirichlet = new DirichletMixture(userParameters->getMixtureName());
    DirichletPriors = new float[ALPHABET_SIZE];
    dirichlet->getDirichletPriors(DirichletPriors);
    
    //for (int i=0; i<ALPHABET_SIZE; i++) 
    //cout << DirichletPriors[i] << endl;
    

    /*
     * When we create a profile from a single sequence, there is only one residue
     * with weight=1. we can precompute the posterior probability distribution
     * for all 20 amino acids, which can then be reused later. 
     */ 

    int i, j;
    float counts[ALPHABET_SIZE];
    double prob[ALPHABET_SIZE];

    for (i=0; i<ALPHABET_SIZE; i++) {
      // create the count vector for one residue i
      for (j=0; j<ALPHABET_SIZE; j++) 
	counts[j] = 0;
      counts[i] = 1;

      dirichlet->calcPosteriorProb(counts, prob);
      
      for (j=0; j<ALPHABET_SIZE; j++) {
	SingleResidueProbMatrix.setAt(i, j, prob[j]);
	//cout << SingleResidueProbMatrix.getAt(i,j) << ",";
      }
      //cout << endl;
    }

  }
  
  /*
   * Read in all sequences in the input alignment and create a node 
   * from each sequence
   */
  void
  BETEAlgorithm::initializeNodes()
  {
    numSeqs = alignPtr->getNumSeqs();
    numClasses = numSeqs;
    
    cout << "Initializing " << numSeqs << " nodes..." << endl;

    int firstSeq = 0;
    for (int i=0; i<numSeqs; i++) {
      //cout << "Node # " << nodeCount << endl;

      Node* newnode = new Node(nodeCount, firstSeq+i);

      nodeCount++; 

      if (nodeList == NULL) 
	nodeList = newnode;
      else {
	nodeList->prev = newnode;
	newnode->next = nodeList;
	nodeList = newnode;
      }
    }

    cout << "Finished initialization.\n";
  }

  float
  BETEAlgorithm::calcPairDistance(Node* node1, Node* node2)
  {

    int i,j;
    int ncols = node1->getNumCols();
    Column col1, col2;

    /* calculate TRE (total relative entropy) distance
     * TRE = sum_{c=1..l} {i=1..20} p_c^i*log p_c^i/q_c^i + q_c^i*log q_c^i/p_c^i
     */

    // need to consider overlap between 2 profiles
    int overlap = 0;
    
    for (j=0; j<ncols; j++) {
      col1 = node1->getCol(j);
      col2 = node2->getCol(j);

      // both columns cannot be all gaps
      if (col1.residueCount > 0 && col2.residueCount > 0) {
	overlap++;
      }
    }

    int minOverlap = userParameters->getMinOverlap();
    minOverlap = (int) (minOverlap*ncols)/100;
    if ( overlap < minOverlap) {
      //cout << "There is not enough overlap between 2 nodes, thus no distance will be calculated. A very large distance (30) is assigned." << endl;
      return 30;
    }

    // number of non-gap positions
    float dist = 0, p1, p2, x;
    int ngpos = 0;
    j = 0;
    for (int pos=0; pos < alignPtr->getLength(); pos++) {
      // conserved position should not affect the tree topology (i.e., not contribute to distance)
      if (alignPtr->isConserved(pos))
	continue;
      
      col1 = node1->getCol(j);
      col2 = node2->getCol(j);
      j++;

      if (col1.residueCount > 0 && col2.residueCount > 0) {
	ngpos++;

	for (i=0; i<ALPHABET_SIZE; i++) {
	  p1 = col1.probs[i];
	  p2 = col2.probs[i];

	  if ( abs(p1-p2) < 1e-138)
	    continue;

	  x = log(p1)-log(p2);

	  //dist += p1*(log(p1)-log(p2)) + p2*(log(p2)-log(p1));
	  dist += p1*x - p2*x;
	}
      }

    }

    return dist/ngpos;
  }


  double
  BETEAlgorithm::calcEncodingCost()
  {
    /* unconstrained cost to encode the current set of subtrees under 
     * a Dirichlet mixture density
     * EncodingCost = N*log2(S) - sum_p log2 Prob(n_{p,1} n_{p,2} ...n_{p,S}|mixture)
     * p is the postion in the alignment
     * 1, 2, ... S is the index of the subtrees
     */
    
    int i, j, nonzero;
    int numCols = nodeList->getNumCols();
    
    double logp, log_ucc = 0;
    Column col;

    /*
     *
     * unconstrained: product of sums
     * P(n_{c,1}) * P(n_{c,2}) * ...
     * assuming different subfamilies have been independently drawn from the dirichlet distribution
     * Different from the case of calcalating TRE, conserved positions should be considered -
     * such groups of sequences are more likely to correspond to a single subfamily
     */

    for (j=0; j<numCols; j++) {
      for (Node* node=nodeList; node; node=node->next) {
	col = node->getCol(j);
	
	if (col.residueCount > 0) {
	  /* use raw counts */
	  dirichlet->logProbOrderedCounts(col.counts, logp);
	  log_ucc += logp; // product of sum
	}
      }
    }

    //cout << "Final log_ucc = " << log_ucc << endl;
    double cost = numSeqs * log2(numClasses) - log_ucc;
    return cost;
  }
  
}
