#include "MSA.h"
#include <algorithm>
#include <unordered_map>


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




MSA MSA::apply_mask_MSA(const string& mask) const{
	MSA result;
	string line;
	for(int i(0);i<(int)text.size();++i){
		line=text[i];
		line=apply_mask(line,mask);
		result.add_sequence(line);
	}
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
				name.push_back(ligne.substr(1, ligne.size()-1));
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
		//std::cout << totale_conserved_nuc << '\n';
	}
	return consensus_seq;
}



string MSA::conversion_IG(string nucleotides_and_gap) const{
	string result("");
	bool gap(false);
	string nucleotides("");
	for (int i = 0; i < (int)nucleotides_and_gap.size(); i++) {
		if ('-' == nucleotides_and_gap[i]){
			gap = true;
		}
		else{
			nucleotides += nucleotides_and_gap[i];
		}
	}

	if (nucleotides.size() <= 1) {
		if (nucleotides != ""){
				result = nucleotides[0];
		}
		else{
			result = "-";
		}
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
		string IUPAC_code =	"RYSWKMBDHVN";
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

			if(nuc_find == (int)nucleotides.size()){
				result = IUPAC_code[i];
				code_found = true;
			}
			i++;
		}
	}

	if (gap) {
		result = tolower(result[0]);
	}

	//std::cout <<  "[" << nucleotides_and_gap << "] [" << nucleotides << "] [" << result << "]" << '\n';
	return result;
}




string remove_gap(const string& str){
	string res;
	for(int i(0);i<(int)str.size();++i){
		if(str[i]!='-'){
			res+=str[i];
		}
	}
	return res;
}



bool char_in_string(char c, const string& str){
	c=toupper(c);
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
					if(char_in_string(c,"RWMDHVN")){//THIS ALSO MATCH IF c is LOWERCASE
						//WE DO NOTHING, THIS SHOULD BE THE RIGHT NUC
					}else{
						//WE HAVE NO IDEA
						seq[i]='n';
					}
				break;
				case 'C':
					if(char_in_string(c,"YSMBHVN")){
						//WE DO NOTHING, THIS SHOULD BE THE RIGHT NUC
					}else{
						//WE HAVE NO IDEA
						seq[i]='n';
					}
				break;
				case 'G':
					if(char_in_string(c,"RSKBDVN")){
						//WE DO NOTHING, THIS SHOULD BE THE RIGHT NUC
					}else{
						//WE HAVE NO IDEA
						seq[i]='n';
					}
				break;
				case 'T':
					if(char_in_string(c,"YWKBDHN")){
						//WE DO NOTHING, THIS SHOULD BE THE RIGHT NUC
					}else{
						//WE HAVE NO IDEA
						seq[i]='n';
					}
				break;
				case '-':
					if(islower(c)){
						//WE DO NOTHING, THIS SHOULD BE THE RIGHT NUC
					}else{
						//WE HAVE NO IDEA
						seq[i]='n';
					}
				break;
			}
		}
	}
	return seq;
}



void MSA::count_agreements() const{
	int size = (int)text.size();
	int matrice[size][size];
	//Initialize the matrix
	for(int i=0; i<size ; i++){
		for(int j=0; j<size ; j++){
			if (i == j) {
				matrice[i][j]=length;
			}
			else{
      	matrice[i][j]=0;
		 	}
		}
	}

	//fill the matrix
	for(int i(0);i<length;++i){
		for(int j(0);j<size;++j){
			if (text[j][i] != 'N'){
				for (int k = j+1; k<size; k++) {
					if (text[j][i] == text[k][i]){
						matrice[j][k]++;
					}
				}
			}
		}
	}

	//Matrix display
	for(int i=0; i<size ; i++){
		std::cout << name[i][0] << name[i][1] << name[i][2] <<'\t';
	}
	std::cout << '\n';

	for(int i=0; i<size ; i++){
		for(int j=0; j<size ; j++){
       cout<<matrice[i][j]<<"\t";
		}
		std::cout << "\t"  << name[i] << '\n';
	}

	//Top three
	int score[3];
	string name_score[3];
	score[0]=score[1]=score[2]=0;
	for(int i=0; i<size ; i++){
		for(int j=i+1; j<size ; j++){
       if(matrice[i][j] > score[0]){
				 score[2] = score[1];
				 score[1] = score[0];
				 score[0] = matrice[i][j];

				 name_score[2] = name_score[1];
				 name_score[1] = name_score[0];
				 name_score[0] = name[i] + " & " + name[j];
			 }
			 else{
				 if (matrice[i][j] > score[1]) {
				 	score[2] = score[1];
				 	score[1] = matrice[i][j];

					name_score[2] = name_score[1];
					name_score[1] = name[i] + " & " + name[j];
				 }
				 else{
					 if (matrice[i][j] > score[2]) {
					 	score[2] = matrice[i][j];

						name_score[2] = name[i] + " & " + name[j];
					 }

				 }
			 }
		}
	}

	std::cout << "\nTop Three:" << '\n';
	std::cout << score[0] << " " << name_score[0] << '\n';
	std::cout << score[1] << " " << name_score[1] << '\n';
	std::cout << score[2] << " " << name_score[2] << '\n';
}

void MSA::compare_consensus() const{
	std::cout << "\nConsensus sequence comparison:"<< '\n';
	int unique_base[lines];
	string unique_base_position[lines];
	for (int i = 0; i < lines; i++) {
		unique_base[i] = 0;
		unique_base_position[i] = "";
	}


	for(int i(0);i<lines;++i){
		//std::cout << "\npos" << i << '\n';
		for(int j(0);j<length;++j){
			int k(0);
			bool unique(true);
			while ( unique && k < lines) {
				if (k != i && text[i][j] == text[k][j] ) {
					unique = false;
				}
				k++;
			}
			if (unique){
				unique_base[i]++;
				unique_base_position[i] += to_string(j) + ",";
			}
		}
	}

	std::cout << "name\tunique\tposition" << '\n';
	for(int i=0; i<lines ; i++){
		std::cout << name[i] << "\t" << unique_base[i] \
		<< "\t" << unique_base_position[i] << '\n';
	}
}


int char2int(char c){
	switch(c){
		case 'A':return 0;
		case 'a': return 0;
		case 'C':return 1;
		case 'c': return 1;
		case 'G':return 2;
		case 'g': return 2;
		case 'T':return 3;
		case 't': return 3;
	}
	return -1;
}


struct next_nuc{
	char first_nuc;
	char second_nuc;
	int weight;
};

bool compareWeight(const next_nuc &a, const next_nuc &b)
{
    return a.weight > b.weight;
}

void MSA::get_diploid(){
	string haplotype1,haplotype2;
	//FOREACH POSITION
	for(int i(0);i<length-1;++i){
		cout<<"Next nucleotide"<<endl;
		unordered_map<string,uint> two_mer_count;
		//FOREACH SEQUENCE
		for(int j(0);j<lines;++j){
			two_mer_count[text[j].substr(i,2)]++;
		}
		//we select the heaviest next nuc for A,C,G,T
		next_nuc nextA={'A','@',0},nextC={'C','@',0},nextG={'G','@',0},nextT={'T','@',0};
		for(auto& it: two_mer_count) {
			//cout<<"bimer "<<it.first<<endl;
			//cout<<"bimer "<<it.first[0]<<" "<<it.first[1]<<endl;
			switch (it.first[0]){
				case 'A':
					if(it.second>nextA.weight){
						cout<<"After A its probably a "<<it.first[1]<<" instead of a "<<nextA.second_nuc<<" "<<nextA.weight<<" "<<it.second<<endl;
						nextA.weight=it.second;
						nextA.second_nuc=it.first[1];
					}
					break;
				case 'C':
					if(it.second>nextC.weight){
						cout<<"After C its probably a "<<it.first[1]<<" instead of a "<<nextC.second_nuc<<" "<<nextC.weight<<" "<<it.second<<endl;
						nextC.weight=it.second;
						nextC.second_nuc=it.first[1];
					}
					break;
				case 'G':
					if(it.second>nextG.weight){
						cout<<"After G its probably a "<<it.first[1]<<" instead of a "<<nextG.second_nuc<<" "<<nextG.weight<<" "<<it.second<<endl;
						nextG.weight=it.second;
						nextG.second_nuc=it.first[1];
					}
					break;
				case 'T':
					if(it.second>nextT.weight){
						cout<<"After T its probably a "<<it.first[1]<<" instead of a "<<nextT.second_nuc<<" "<<nextT.weight<<" "<<it.second<<endl;
						nextT.weight=it.second;
						nextT.second_nuc=it.first[1];
					}
					break;
				default:
				cout<<"Do not know what to do with this "<<it.first[0]<<endl;
				break;
			}
		}

		if(haplotype1.empty()){
			vector<next_nuc> bimers={nextA,nextT,nextC,nextG};
			sort(bimers.begin(),bimers.end(),compareWeight);
			haplotype1+=bimers[0].first_nuc;
			haplotype1+=bimers[0].second_nuc;
			haplotype2+=bimers[1].first_nuc;
			haplotype2+=bimers[1].second_nuc;
			cout<<bimers[0].weight<<" "<<bimers[1].weight<<endl;
			cout<<"init:	"<<haplotype1<<" "<<haplotype2<<endl;
		}else {
			cout<<"A->"<<nextA.second_nuc<<" C->"<<nextC.second_nuc<<" G->"<<nextG.second_nuc<<" T->"<<nextT.second_nuc<<endl;
			switch(haplotype1[haplotype1.size()-1]){
				case 'A':
					haplotype1+=nextA.second_nuc;
					break;
				case 'C':
					haplotype1+=nextC.second_nuc;
					break;
				case 'G':
					haplotype1+=nextG.second_nuc;
					break;
				case 'T':
					haplotype1+=nextT.second_nuc;
					break;
				default:
				cout << "weird nuc "<<haplotype1[haplotype1.size()-1]<<endl;
				cout<<"haplotype1: "<<haplotype1<<endl;
				cout<<"haplotype2: "<<haplotype2<<endl;
				haplotype1.clear();
				haplotype2.clear();
				break;
			}
			switch(haplotype2[haplotype2.size()-1]){
				case 'A':
					haplotype2+=nextA.second_nuc;
					break;
				case 'C':
					haplotype2+=nextC.second_nuc;
					break;
				case 'G':
					haplotype2+=nextG.second_nuc;
					break;
				case 'T':
					haplotype2+=nextT.second_nuc;
					break;
				default:
				cout << "weird nuc "<<haplotype2[haplotype2.size()-1]<<endl;
				cout<<"haplotype1: "<<haplotype1<<endl;
				cout<<"haplotype2: "<<haplotype2<<endl;
				haplotype1.clear();
				haplotype2.clear();
				break;
			}
			if(haplotype1[haplotype1.size()-1]==haplotype2[haplotype2.size()-1]){
				cout<<"haplotypes convergence "<<endl;
				cout<<"haplotype1: "<<haplotype1<<endl;
				cout<<"haplotype2: "<<haplotype2<<endl;
				haplotype1.clear();
				haplotype2.clear();
				cin.get();
			}

		}
	}
	cout<<"Haplotype1: " <<haplotype1<<endl;
	cout<<"Haplotype2: " <<haplotype2<<endl;
}

int match_column_nucleotide(char nuc){
	switch ( nuc )
	{
	case 'A':
	  return 0;
	case 'C':
	  return 1;
	case 'G':
		return 2;
	case 'T':
		return 3;
	default:
	  return -1;
	}
}

void MSA::calculates_distance_matrix(){
	int matrice[length*4][length*4];
	for (int i = 0; i < length*4; i++) {
		for (int j = 0; j < length*4; j++) {
			matrice[i][j]=0;
		}
	}
	for (int i = 0; i < lines; i++) {
		for (int j = 0; j < length; j++) {
			for (int k = j+1; k < length; k++) {
				int idNucJ=match_column_nucleotide(text[i][j])+4*j;
				int idNucK=match_column_nucleotide(text[i][k]) +4*k;
				matrice[idNucJ][idNucK] ++;
				matrice[idNucK][idNucJ] ++;
			}
		}
	}
}
