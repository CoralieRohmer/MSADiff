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
	MSA msa_masked(msa_fasta.apply_mask_MSA(consensus));
	// msa_masked.printing();
	// cin.get();
	MSA result = msa_masked.get_compacted();
	// result.printing();
		// cout<< "calculates_distance_matrix"<<endl;

	auto matrix=result.calculates_distance_matrix();
			// cout<< "OK"<<endl;

	//result.printing();
  //result.get_diploid();

	uint score(0);
	cout<< "get_best_consensus_from_prefix"<<endl;
	string ringo(get_best_consensus_from_prefix("T",matrix,score));
	cout<<ringo<<endl;

	string lestringo(get_best_consensus(matrix,1));
	cout<<lestringo<<endl;

	vector<string> haplo_snip;
	haplo_snip.push_back("AACA");
	haplo_snip.push_back("GACC");

	vector<string> haplo = msa_masked.haplotype_merge(haplo_snip);
	std::cout << "Haplotype1 : "<< haplo[0] << '\n';
	std::cout << "Haplotype2 : "<< haplo[1] << '\n';

	/* MSA msa_stack = MSA();
	msa_stack.parser_fasta("test/multi_msa_stack.fasta");
	msa_fasta.count_agreements();
	msa_fasta.compare_consensus(); */
  return 0;
}
