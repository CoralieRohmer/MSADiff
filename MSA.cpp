#include "MSA.h"
#include <algorithm>

using namespace std;



void MSA::add_sequence(const string& str){
	text.push_back(str);
	length=str.size();
	++lines;
}



void MSA::concat_sequence(const string& str,int i){
	if((int)text.size()<=i){
		text.resize(i+1);
		lines=i+1;
	}
	text[i]+=str;
	length=text[i].size();
}



bool MSA::check_sane() const{
	bool result(true);
	if(lines!=(int)text.size()){
		cout<<"Incorrect number of lines:	"<<text.size()<<" instead of:	"<<lines<<endl;
		result=false;
	}
	for(int i(0);i<(int)text.size();++i){
		if((int)text[i].size()!=length){
			result=false;
			cout<<"Incorrect length line:	"<<i<<" length:	"<<text[i].size()<<" instead of:	"<<length<<endl;

		}
		for(int j(0);j<(int)text[i].size();++j){
			if(alphabet.find(text[i][j])== std::string::npos){
				cout<<"Incorrect char:	"<<text[i][j]<<endl;
			}
		}
	}
	return result;
}



void MSA::printing() const{
	check_sane();
	for(int i(0);i<(int)text.size();++i){
		cout<<text[i]<<endl;
	}
}



bool MSA::perfect_column(int indice)const {
	char c(text[0][indice]);
	for(int i(1);i<(int)text.size();++i){
		if(c!=text[i][indice]){
			return false;
		}
	}
	return true;
}


char MSA::major_nuc(int indice)const{
	int scores[alphabet.size()]={0};
	char result('Z');
	for(int i(0);i<(int)text.size();++i){
		scores[alphabet.find(text[i][indice])]++;
	}
	int max_score(0);
	for(int i(0);i<(int)alphabet.size();++i){
		if(scores[i]>max_score){
			result=alphabet[i];
			max_score=scores[i];
		}
	}
	return result;
}



bool MSA::perfect_column_quasi(int indice,int max_errors)const {
	char c(major_nuc(indice));
	int errors(0);
	for(int i(0);i<(int)text.size();++i){
		if(c!=text[i][indice]){
			errors++;
		}
	}
	return errors<=max_errors;
}



MSA MSA::get_compacted() const{
	MSA result;
	for(int i(0);i<length;++i){
		if(not perfect_column(i)){
			for(int j(0);j<(int)text.size();++j){
				result.concat_sequence(text[j].substr(i,1),j);
			}
		}else{
		}
	}
	cout<<"Compaction done, length before:	"<<length<<" after:	"<<result.length<<endl;
	return result;
}



MSA MSA::get_compacted_quasi(int max_errors ) const{
	MSA result;
	for(int i(0);i<length;++i){
		if(not perfect_column_quasi(i,max_errors)){
			for(int j(0);j<(int)text.size();++j){
				result.concat_sequence(text[j].substr(i,1),j);
			}
		}else{
		}
	}
	cout<<"Compaction done, length before:	"<<length<<" after:	"<<result.length<<endl;
	return result;
}



void MSA::parser_fasta(string file){
	ifstream flux(file.c_str());
	string seq("");
	if(flux){
		srand(time(NULL));
		string ligne("");
		while(getline(flux, ligne)){
			if (ligne[0] == '>'){
				if (seq != ""){
					add_sequence(seq);
				}
				seq = "";
			}
			else{
				seq += ligne;
			}
		}
		add_sequence(seq);
	}
	else{
		cout << "ERROR: Unable to open the file.." << endl;
		exit(0);
	}
}



string MSA::consensus(int threshold) const{
	string consensus_seq("");
	for(int i(0);i<length;++i){
		int scores[alphabet.size()]={0};

		for(int j(0);j<(int)text.size();++j){
			scores[alphabet.find(text[j][i])]++;
		}

		int max_score(0);
		char conserved_nuc(' ');
		bool equality(false);
		for(int j(0);j<(int)alphabet.size();++j){
			if(scores[j]>max_score){
				conserved_nuc=alphabet[j];
				max_score=scores[j];
				equality = false;
			}
			else{
				if (scores[j] == max_score){
					equality = true;
				}
			}
		}

		if ( (max_score*100)/(int)text.size() < threshold || equality) {
			conserved_nuc = 'N';
		}
		consensus_seq += conserved_nuc;
	}
	return consensus_seq;
}



string MSA::consensus_IG(int threshold) const{
	string consensus_seq("");
	for(int i(0);i<length;++i){
		int scores[alphabet.size()]={0};

		for(int j(0);j<(int)text.size();++j){
			scores[alphabet.find(text[j][i])]++;
		}
		int totale_score(0);
		string totale_conserved_nuc("");
		while (totale_score < threshold){
			int max_score(0);
			string conserved_nuc("");
			bool equality(false);

			for(int j(0);j<(int)alphabet.size();++j){
				bool already_treated_nuc(false);
				for (int k = 0; k < (int)totale_conserved_nuc.size(); k++) {
					if (alphabet[j] == totale_conserved_nuc[k]){
						already_treated_nuc = true;
					}
				}
				if (!already_treated_nuc){
					if (scores[j]>max_score){
						conserved_nuc=alphabet[j];
						max_score=scores[j];
						equality = false;
					}
					else{
						if (scores[j] == max_score){
							conserved_nuc += alphabet[j];
							equality = true;
						}
					}
				}
			}

			if (equality){
				totale_score += ((max_score*100)/(int)text.size())*conserved_nuc.size();
			}
			else{
				totale_score += (max_score*100)/(int)text.size();
			}

			totale_conserved_nuc += conserved_nuc;

		}
		consensus_seq += conversion_IG(totale_conserved_nuc);
	}
	return consensus_seq;
}



char MSA::conversion_IG(string nucleotides) const{
	char result(' ');
	if (nucleotides.size() == 1) {
		result = nucleotides[0];
	}
	else{
		//std::cout << nucleotides << '\n';
		bool gap(false);
		int i(0);
		while( !gap && i < (int)nucleotides.size()) {
			if ('-' == nucleotides[i]){
				gap = true;
			}
			i++;
		}
		string bases[11] = {"AG","CT","CG","AT","GT","AC","CGT","AGT","ACT","ACG","ACGT"};
		string IUPAC_code =	"R     Y    S    W    K   M     B     D     H     V     N";
		bool code_found(false);
		i=0;

		while(!code_found && i < 11) {
			int nuc_find(0);
			for (int j = 0; j < (int)nucleotides.size(); j++) {
				if (nucleotides.size() == bases[i].size()) {
					for (int k = 0; k < (int)bases[i].size(); k++) {
						if (nucleotides[j] == bases[i][k]){
							nuc_find++;
						}
					}
				}
			}
			int nucleotide_to_be_found = nucleotides.size();
			if (gap){
				nucleotide_to_be_found--;
			}
			if(nuc_find == nucleotide_to_be_found){
				result = IUPAC_code[i];
				code_found = true;
			}
			i++;
		}

		if (gap) {
			result = tolower(result);
		}
	}

	return result;
}



bool char_in_string(char c, const string& str){
	for(int i(0);i<(int)str.size();++i){
		if(str[i]==c){
			return true;
		}
	}
	return false;
}



string apply_mask(string seq, const string& mask){
	for(int i(0);i<(int)seq.size();++i){
		char c(mask[i]);
		//WE TRUST THE CONSENSUS FOR THIS NUC
		if(char_in_string(c,"ACTG-")){
			seq[i]=c;
		}else{
			switch (seq[i]){
				case 'A':
				if(char_in_string(c,"RWMDHVN")){
					//WE DO NOTHING, THIS SHOULD BE THE RIGHT NUC
				}else{
					//WE HAVE NO IDEA
					seq[i]='N';
				}
				break;
				case 'C':
				if(char_in_string(c,"YSMBHVN")){
					//WE DO NOTHING, THIS SHOULD BE THE RIGHT NUC
				}else{
					//WE HAVE NO IDEA
					seq[i]='N';
				}
				break;
				case 'G':
				if(char_in_string(c,"RSKBDVN")){
					//WE DO NOTHING, THIS SHOULD BE THE RIGHT NUC
				}else{
					//WE HAVE NO IDEA
					seq[i]='N';
				}
				break;
				case 'T':
				if(char_in_string(c,"YWKBDHN")){
					//WE DO NOTHING, THIS SHOULD BE THE RIGHT NUC
				}else{
					//WE HAVE NO IDEA
					seq[i]='N';
				}
				break;
				case '-':
				if(islower(c)){
					//WE DO NOTHING, THIS SHOULD BE THE RIGHT NUC
				}else{
					//WE HAVE NO IDEA
					seq[i]='N';
				}
			}
		}
	}
	return seq;
}
