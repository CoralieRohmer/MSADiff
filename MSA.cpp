#include "MSA.h"
#include <algorithm>
#include <unordered_map>



using namespace std;

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
		case '-': return 4;
	}
	return -1;
}

char int2char(int i){
	switch(i){
		case 0 :return 'A';
		case 1 :return 'C';
		case 2 :return 'G';
		case 3 :return 'T';
		case 4 : return '-';
	}
	return -1;
}

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

double MSA::pb_colonne(int i, int ploidy, double pb_error){
	double pb_haplo((double)1/ploidy);
	vector<int> occ(5,0);
	for(int j(0);j<lines;++j){
		occ[char2int(text[j][i])]++;
	}
	vector<int> max_occ(ploidy,0);
	vector<char> max_nuc(ploidy,0);
	//search majority nucleotides
	for (int j = 0; j < 5; j++) {
		int min(max_occ[0]);
		int min_k(0);
		for (int k = 0; k < ploidy; k++) {
			//search for the smallest maximum that exceeds our occurrence (to manage polyploidy)
			if (occ[j] > max_occ[k]) {
				if (max_occ[k] < min) {
					min=max_occ[k];
					min_k=k;
				}
			}
		}
		if (occ[j] > max_occ[min_k]) {
			max_occ[min_k]=occ[j];
			max_nuc[min_k]=int2char(j);
		}
	}

	double pb_column(1);
	for (int j = 0; j < lines; j++) {
		if (*find(max_nuc.begin(), max_nuc.end(), text[j][i]) == text[j][i]) {
			pb_column*=pb_haplo*(1-pb_error) + (ploidy-1)*pb_haplo*pb_error;
		}
		else{
			pb_column*=ploidy*pb_haplo*pb_error;
		}
	}
	return pb_column;
}

MSA MSA::cleaning(double pb_error){
	MSA result;
	//SEARCH COLUMNS TO BE CORRECTED
	vector<pair<int,char>> column_to_clean;
	for (int i = 0; i < length; i++) {
		vector<int> occ(5,0);
		double pb_haplo = pb_colonne(i, 1, pb_error);
		double pb_diploide = pb_colonne(i, 2, pb_error);
		if (pb_haplo > pb_diploide) {
			for (int j = 0; j < lines; j++) {
				occ[char2int(text[j][i])]++;
			}
			int max_occ(0);
			int maj_nuc(0);
			for (int j = 0; j < 5; j++) {
				if (occ[j] > max_occ) {
					max_occ = occ[j];
					maj_nuc = j;
				}
			}
			if (max_occ != lines) {
				column_to_clean.push_back({i,int2char(maj_nuc)});
			}
		}
	}
	//CORRECTION OF COLUMNS
	string line;
	int nb_column(column_to_clean.size());
	std::cout << "Cleaning steps: " << nb_column <<" sequencing error(s) found and corrected." << '\n';
	for (int i = 0; i < lines; i++) {
		line = text[i];
		for (int j = 0; j < nb_column; j++) {
			line[column_to_clean[j].first] = column_to_clean[j].second;
		}
		result.add_sequence(line);
	}
	return result;
}

MSA MSA::apply_mask_MSA(const string& mask) const{
	MSA result;
	string line;
	for(int i(0);i<(int)text.size();++i){
		line=text[i];
		std::cout << line << '\n';
		line=apply_mask(line,mask);
		std::cout << line << "\n\n";
		result.add_sequence(line);
	}
	return result;
}

string apply_mask(string seq, const string& mask){
	for(int i(0);i<(int)seq.size();++i){
		char c(mask[i]);
		//WE TRUST THE CONSENSUS FOR THIS NUC
		if(char_in_string(c,"ACTG-")){
			seq[i]=c;
		}
	}
	return seq;
}

bool char_in_string(char c, const string& str){
	for(int i(0);i<(int)str.size();++i){
		if(str[i]==c){
			return true;
		}
	}
	return false;
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
	}
	return consensus_seq;
}



vector<uint> MSA::shuffle_msa(){
	vector<pair<double,uint>> ratios;
	for(int i(0);i<length;++i){
		vector<int> scores(5,0);
		for(int j(0);j<lines;++j){
			switch(text[j][i]){
				case 'A':
					scores[0]++;
					break;
				case 'C':
					scores[1]++;
					break;
				case 'G':
					scores[2]++;
					break;
				case 'T':
					scores[3]++;
					break;
				case '-':
					scores[4]++;
					break;
			}
		}
		double cmax(0),cmin(0);
		for(uint j(0);j<5;++j){
			if(scores[j]>cmax){
				if(cmax>cmin){
					cmin=cmax;
				}
				cmax=scores[j];
			}else{
				if(scores[j]>cmin){
					cmin=scores[j];
				}
			}
		}
		ratios.push_back({cmax/cmin,ratios.size()});
	}
	sort(ratios.begin(),ratios.end());
	vector<string> new_text(text);
	vector<uint> result;
	for(uint i(0);i<ratios.size();++i){
		for(int j(0);j<lines;++j){
			new_text[j][i]=text[j][ratios[i].second];
		}
		result.push_back(ratios[i].second);
	}
	text=new_text;
	return result;
}

vector<uint> MSA::order_colonnes(double pb_error){
	vector<pair<double,uint>> pb_colonnes;
	//CALCULATES THE PROBA FOR EACH COLUMN
	for(int i(0);i<length;++i){
		pb_colonnes.push_back({1-pb_colonne(i,2,pb_error),pb_colonnes.size()});
	}

	//SORT THE COLUMNS ACCORDING THEIR PROBA
	sort(pb_colonnes.begin(),pb_colonnes.end());
	vector<string> new_text(text);
	vector<uint> result;
	for(uint i(0);i<pb_colonnes.size();++i){
		for(int j(0);j<lines;++j){
			new_text[j][i]=text[j][pb_colonnes[i].second];
		}
		result.push_back(pb_colonnes[i].second);
	}
	text=new_text;
	return result;
}

void MSA::reverse_shuffle_msa(vector<uint>& result){
	vector<string> new_text(text);
	for(uint i(0);i<result.size();++i){
		for(int j(0);j<lines;++j){
			new_text[j][result[i]]=text[j][i];
		}
	}
	text=new_text;
}



void reverse_shuffle_haplotype(vector<uint>& result,string& haplotype){
	string new_haplotype=haplotype;
	for(uint i(0);i<result.size();++i){
		new_haplotype[result[i]]=haplotype[i];
	}
	haplotype=new_haplotype;
}



void reverse_shuffle_haplotypes(vector<uint>& result,vector<string>& haplotypes){
	for(uint i(0);i<haplotypes.size();++i){
		reverse_shuffle_haplotype(result,haplotypes[i]);
	}
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





struct next_nuc{
	char first_nuc;
	char second_nuc;
	int weight;
};



bool compareWeight(const next_nuc &a, const next_nuc &b)
{
    return a.weight > b.weight;
}


uint score_sequence(const string& sequence, vector<vector<uint>> matrix){
	uint score=0;
	for(uint i=0; i<sequence.length();++i){
		uint nuc_from(char2int(sequence[i]));
		for(uint j=0; j<sequence.size();++j){
			uint nuc_to(char2int(sequence[j]));
			score+=matrix[5*i+nuc_from][5*j+nuc_to];
		}
	}
	return score;
}



string get_best_consensus_from_prefix(const string& prefix, vector<vector<uint>> matrix, uint& score){
	uint nuc_number(matrix.size()/5);
	string result=prefix;
	string acgt="ACGT-";
	string best_local_solution;
	for(uint i(prefix.size());i<nuc_number;++i){
		for(uint i_nuc(0);i_nuc<5;++i_nuc){
			uint local_score(score_sequence(result+acgt[i_nuc],matrix));
			if(local_score>score){
				best_local_solution=result+acgt[i_nuc];
				score=local_score;
			}
		}
		result=best_local_solution;
	}
	std::cout << result << '\n';
	return result;
}


string get_best_consensus_aux(const string& prefix, vector<vector<uint>> matrix, uint prefix_size, uint& score){
	string acgt="ACGT-";
	string result;
	for(uint i_nuc(0);i_nuc<5;++i_nuc){
		string local_result;
		uint local_score(0);
		string local_prefix(prefix+acgt[i_nuc]);
		if(local_prefix.size()==prefix_size){
			local_result=get_best_consensus_from_prefix(local_prefix,matrix,local_score);
			cout<<"Found "<<local_result<<" scoring: "<<local_score<<endl;
			if(local_score>score){
				score=local_score;
				result=local_result;
			}
		}else{
			local_result=get_best_consensus_aux(local_prefix,matrix,prefix_size,local_score);
			if(local_score>score){
				score=local_score;
				result=local_result;
			}
		}
	}
	return result;
}


vector<pair<uint64_t,string>> get_best_consensus_vector(const string& prefix, vector<vector<uint>> matrix, uint prefix_size, uint& score){
	string acgt="ACGT-";
	vector<pair<uint64_t,string>> result;
	for(uint i_nuc(0);i_nuc<5;++i_nuc){
		string local_result;
		uint local_score(0);
		string local_prefix(prefix+acgt[i_nuc]);
		if(local_prefix.size()==prefix_size){
			local_result=get_best_consensus_from_prefix(local_prefix,matrix,local_score);
			cout<<"Found "<<local_result<<" scoring: "<<local_score<<endl;
			score=local_score;
			if(not local_result.empty()){
				result.push_back({local_score,local_result});
			}
		}else{
			local_result=get_best_consensus_aux(local_prefix,matrix,prefix_size,local_score);
			score=local_score;
			if(not local_result.empty()){
				result.push_back({local_score,local_result});
			}
		}
	}
	return result;
}



vector<string> get_best_consensus_vector(vector<vector<uint>> matrix, uint prefix_size,uint number_consensus){
	uint score(0);
	auto consensus=get_best_consensus_vector("", matrix, prefix_size, score);
	sort(consensus.begin(),consensus.end());
	reverse(consensus.begin(),consensus.end());
	vector<string> result;
	cout<<"sorted result:	"<<consensus.size()<<endl;
	for(uint i(0);i<min((uint)consensus.size(),number_consensus);++i){
		cout<<"score:	"<<consensus[i].first<<endl;
		cout<<"seq:	"<<consensus[i].second<<endl;
		result.push_back(consensus[i].second);

	}
	return result;
}




string get_best_consensus(vector<vector<uint>> matrix, uint prefix_size){
	uint score(0);
	string result=get_best_consensus_aux("", matrix, prefix_size, score);
	cout<<"score:	"<<score<<endl;
	return result;
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
					if((int)it.second>nextA.weight){
						cout<<"After A its probably a "<<it.first[1]<<" instead of a "<<nextA.second_nuc<<" "<<nextA.weight<<" "<<it.second<<endl;
						nextA.weight=it.second;
						nextA.second_nuc=it.first[1];
					}
					break;
				case 'C':
					if((int)it.second>nextC.weight){
						cout<<"After C its probably a "<<it.first[1]<<" instead of a "<<nextC.second_nuc<<" "<<nextC.weight<<" "<<it.second<<endl;
						nextC.weight=it.second;
						nextC.second_nuc=it.first[1];
					}
					break;
				case 'G':
					if((int)it.second>nextG.weight){
						cout<<"After G its probably a "<<it.first[1]<<" instead of a "<<nextG.second_nuc<<" "<<nextG.weight<<" "<<it.second<<endl;
						nextG.weight=it.second;
						nextG.second_nuc=it.first[1];
					}
					break;
				case 'T':
					if((int)it.second>nextT.weight){
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

vector<vector<uint>> MSA::calculates_distance_matrix(){
	vector<uint> colonne(length*5,0);
	vector<vector<uint>> matrix(length*5,colonne);
	for (int i = 0; i < lines; i++) {
		for (int j = 0; j < length; j++) {
			for (int k = j+1; k < length; k++) {
				int idNucJ=char2int(text[i][j])+5*j;
				int idNucK=char2int(text[i][k])+5*k;
				matrix[idNucJ][idNucK] ++;
				matrix[idNucK][idNucJ] ++;
			}
		}
	}
	return matrix;
}

vector<string> MSA::haplotype_merge(vector<string> haplo_snip){
	vector<string> haplo(haplo_snip.size());
	int snip_to_add=0;
	for(int i(0);i<length;++i){
		bool snip = false;
		char previous = 'N';
		int j=0;
		while (!snip && j<(int)text.size()) {
			if(previous!= 'N' && text[j][i] !=  previous){
				snip=true;
			}
			previous=text[j][i];
			j++;
		}
		if(snip){
			for (int k = 0; k < (int)haplo_snip.size(); k++) {
				haplo[k] += haplo_snip[k][snip_to_add];
			}
			snip_to_add++;
		}
		else{
			for (int k = 0; k < (int)haplo_snip.size(); k++) {
				haplo[k] += text[0][i];
			}
		}
	}
	return haplo;
}

vector<uint> create_one_cluster(int i,vector<vector<uint>> matrix,double pb_error, int length){
	vector<uint> cluster;
	vector<uint> treated;
	treated.push_back(i);
	while(treated.size() > 0){
		uint nb_read_treated=treated.size();
		for (uint j = 0; j < nb_read_treated; j++) {
			int read=treated[treated.size()-1];
			cluster.push_back(read);
			treated.pop_back();
			for (uint k = 0; k < matrix.size(); k++) {
				if (matrix[k][read] < 2.5*pb_error*length){
					if (find(cluster.begin(), cluster.end(), k) == cluster.end() &&
							find(treated.begin(), treated.end(), k) == treated.end()) {
						treated.push_back(k);
					}
				}
			}
		}
	}
	return cluster;
}

void MSA::clustering(double pb_error){
	vector<uint> colonne(lines,0);
	vector<vector<uint>> matrix(lines,colonne);
	for (int i = 0; i < length; i++) {
		for (int j = 0; j < lines-1; j++) {
			for (int k = j+1; k < lines; k++) {
				if (text[j][i] != text[k][i]) {
					matrix[j][k]++;
					matrix[k][j]++;
				}
			}
		}
	}

	vector<vector<uint>> clusters;
	bool all_treated(0);
	int i(0);
	while (!all_treated){
		clusters.push_back(create_one_cluster(i,matrix,pb_error,length));
		bool next(1);
		while(next && i < lines){
			i++;
			uint j(0);
			bool bfind(0);
			while (!bfind && j < clusters.size()) {
				if (find(clusters[j].begin(), clusters[j].end(), i) != clusters[j].end()) {
					bfind=1;
				}
				j++;
			}
			if (!bfind) {
				next=0;
			}
		}
		if (i >= lines) {
			all_treated=1;
		}
	}
	std::cout << "\n=======Clusters==============" << '\n';
	for (uint i = 0; i < clusters.size(); i++) {
		std::cout << "Cluster n°" << i << "\n{ ";
		for (uint j = 0; j < clusters[i].size(); j++) {
			std::cout << clusters[i][j] << ' ';
		}
		std::cout << '}' << '\n';
	}

}
