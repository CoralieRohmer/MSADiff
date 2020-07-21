#include "MSA.h"



using namespace std;



void MSA::add_sequence(const string& str){
	text.push_back(str);
	length=str.size();
	++lines;
}



void MSA::concat_sequence(const string& str,int i){
	if((int)text.size()<=i){
		//~ cout<<"I had to resize"<<endl;
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



bool MSA::perfect_column_quasi(int indice,int max_errors)const {
	char c(text[0][indice]);
	int errors(0);
	for(int i(1);i<(int)text.size();++i){
		if(c!=text[i][indice]){
			errors++;
		}
	}
	return errors<=max_errors;
}



MSA MSA::get_compacted() const{
	MSA result;
	for(int i(0);i<length;++i){
		//~ cout<<"pos"<<i<<endl;
		if(not perfect_column(i)){
			//~ cout<<"diff found"<<endl;
			for(int j(0);j<(int)text.size();++j){
				//~ cout<<text[j].substr(i,1)<<endl;
				result.concat_sequence(text[j].substr(i,1),j);
			}
		}else{
			//~ cout<<"BORING "<<endl;
		}
	}
	cout<<"Compaction done, length before:	"<<length<<" after:	"<<result.length<<endl;
	return result;
}



MSA MSA::get_compacted_quasi(int max_errors ) const{
	MSA result;
	for(int i(0);i<length;++i){
		//~ cout<<"pos"<<i<<endl;
		if(not perfect_column_quasi(i,max_errors)){
			//~ cout<<"diff found"<<endl;
			for(int j(0);j<(int)text.size();++j){
				//~ cout<<text[j].substr(i,1)<<endl;
				result.concat_sequence(text[j].substr(i,1),j);
			}
		}else{
			//~ cout<<"BORING "<<endl;
		}
	}
	cout<<"Compaction done, length before:	"<<length<<" after:	"<<result.length<<endl;
	return result;
}



