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
	MSA result = msa_fasta.get_compacted_quasi(1);
	result.printing();
	//std::cout << msa_fasta.consensus(10) << '\n';
  return 0;
}
