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
	MSA result = msa_masked.get_compacted();
	result.calculates_distance_matrix();

	//result.printing();
  //result.get_diploid();

	vector<vector<uint>>neo = {
    {0,0,0,0,0,1,0,0},
	{0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0},
	{1,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0}
};

	string ringo(get_best_consensus_from_prefix("A",neo,2));
	cout<<ringo<<endl;

	/* MSA msa_stack = MSA();
	msa_stack.parser_fasta("test/multi_msa_stack.fasta");
	msa_fasta.count_agreements();
	msa_fasta.compare_consensus(); */
  return 0;
}
