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
	string const file_fasta(argv[1]);
	MSA msa_fasta = MSA();

	ifstream flux(file_fasta.c_str());
	string seq("");
	string name_seq("");
	if(flux){   
		srand(time(NULL));
		string ligne("");
		while(getline(flux, ligne)){
			if (ligne[0] == '>'){
				name_seq = ligne;
				if (seq != ""){
					msa_fasta.add_sequence(seq);
				}
				seq = "";
			}
			else{
				seq += ligne;
			}
		}
		msa_fasta.add_sequence(seq);
		msa_fasta.printing();
		MSA result = msa_fasta.get_compacted_quasi(2);
		result.printing();
	}
	else{
		cout << "ERREUR: Impossible d'ouvrir le fichier." << endl;
		exit(0);
	}

  return 0;
}
