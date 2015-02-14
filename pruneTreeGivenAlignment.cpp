//
// File: pruneTreeGivenAlignment.cpp
// Created by: Bastien Boussau
// Created on: Friday, January 20 2011 16:48
//

/*
Copyright or © or Copr. CNRS

This software is a computer program whose purpose is to estimate
phylogenies and evolutionary parameters from a dataset according to
the maximum likelihood principle.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

// From the STL:
#include <iostream>
#include <iomanip>
#include <string.h>

using namespace std;

#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Numeric/DataTable.h>
#include <Bpp/Numeric/Random/RandomTools.h>

// From SeqLib:
#include <Bpp/Seq/Alphabet.all>
#include <Bpp/Seq/Container.all>
#include <Bpp/Seq/Io.all>
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/SequenceTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>

// From PhylLib:
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/Distance.all>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Io/PhylipDistanceMatrixFormat.h>
#include <Bpp/Phyl/Io/Nhx.h>
#include <Bpp/Phyl/Io/Newick.h>
using namespace bpp;

void help()
{
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
  (*ApplicationTools::message << "pruneTreeGivenAlignment parameter1_name=parameter1_value").endLine();
  (*ApplicationTools::message << "      parameter2_name=parameter2_value ... param=option_file").endLine();
  (*ApplicationTools::message).endLine();
  (*ApplicationTools::message << "  Refer to the Bio++ Program Suite Manual for a list of available options.").endLine();
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
}




class Index {
  public:
    double distance;
    unsigned int i1, i2;

  public:
    Index(double dist, unsigned int i, unsigned int j) : distance(dist), i1(i), i2(j) {}

  public:
    bool operator==(const Index& index) const { return distance == index.distance; }
    bool operator<(const Index& index) const { return distance < index.distance; }
};

class Test {
  private:
    unsigned int pos_;

  public:
    Test(unsigned int pos) : pos_(pos) {}
    
  public:
    bool operator()(const Index& index) { return index.i1 == pos_ || index.i2 == pos_; }
};


/************************************************************************
 * Function to removes leaves from a tree.
 ************************************************************************/

void dropLeaves(TreeTemplate<Node> & tree, const std::vector<string> &spToDrop) {
    for (unsigned int i = 0 ; i < spToDrop.size() ; i++) {
        TreeTemplateTools::dropLeaf(tree, spToDrop[i]);
    }
    return;
}



/*Compilation:
 *dynamic linking:
 g++ -lbpp-core -lbpp-seq -lbpp-phyl -o pruneTreeGivenAlignment pruneTreeGivenAlignment.cpp
 *static linking:
 g++ pruneTreeGivenAlignment.cpp -o pruneTreeGivenAlignment /usr/local/lib/libbpp-core.a /usr/local/lib/libbpp-seq.a /usr/local/lib/libbpp-phyl.a
 */



int main(int args, char ** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*           Bio++ Phylogenetic Pruner, version 0.2              *" << endl;
  cout << "* Author: B. Boussau                        Last Modif. 20/01/11 *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;
  
  if(args == 1)
  {
    help();
    return 0;
  }
  
    try {
        
        BppApplication pruneTreeGivenAlignment(args, argv, "pruneTreeGivenAlignment");
        pruneTreeGivenAlignment.startTimer();
        
        //Get sequences:
        Alphabet* alphabet      = SequenceApplicationTools::getAlphabet(pruneTreeGivenAlignment.getParams());
        SequenceContainer* seqs = SequenceApplicationTools::getSequenceContainer(alphabet, pruneTreeGivenAlignment.getParams());
        
        string inputMethod = ApplicationTools::getStringParameter("input.method", pruneTreeGivenAlignment.getParams(), "tree");
        ApplicationTools::displayResult("Input method", inputMethod);
        
        TreeTemplate< Node > * tree =0;
        tree = new TreeTemplate<Node>(*(PhylogeneticsApplicationTools::getTree(pruneTreeGivenAlignment.getParams())));
        
        std::vector <std::string> allNames = tree->getLeavesNames();
        
        vector<string> seqNames = seqs->getSequencesNames();
        std::vector <std::string> seqToDrop;
        VectorTools::diff(allNames, seqNames, seqToDrop);
        
        dropLeaves(*tree, seqToDrop);

        string outTreePath = ApplicationTools::getAFilePath("output.tree.file", pruneTreeGivenAlignment.getParams(), true, false);

	Fasta fasta;
	string name;
	for (unsigned int i = 0 ; i < seqNames.size() ; i++) {
	  name = seqNames[i] + ".fa";
	  filebuf fb;
	  fb.open (name.c_str(),ios::out);
	  ostream os(&fb);
	  
	  fasta.writeSequence(os, seqs->getSequence(seqNames[i]) );
 
	}

        Newick newick;
        newick.write (*tree , outTreePath);
        
        pruneTreeGivenAlignment.done();
    }
  catch (exception& e)
  {
    cout << endl;
    cout << "_____________________________________________________" << endl;
    cout << "ERROR!!!" << endl;
    cout << e.what() << endl;
    return 1;
  }

  return 0;
}

