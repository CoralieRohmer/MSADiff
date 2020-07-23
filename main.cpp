#include "MSA.h"
#include <fstream>
#include <cstring>
#include <string>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <string>
#include <vector>



using namespace std;



int main(int argc, char** argv){
	MSA msa_fasta = MSA();
	msa_fasta.parser_fasta(argv[1]);
	string consensus(msa_fasta.consensus_IG(80));
	cout<<consensus<<endl;
	MSA msa_masked(msa_fasta.apply_mask_MSA(consensus));
	msa_masked.printing();
	//~ MSA result = msa_fasta.get_compacted_quasi(1);
	//~ std::cout << msa_fasta.consensus_IG(80) << '\n';
  return 0;
}
