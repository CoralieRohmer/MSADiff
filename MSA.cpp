#include "MSA.h"



using namespace std;



void MSA::add_sequence(const string& str){
	text.push_back(str);
	length=str.size();
	++lines;
}



void MSA::concat_sequence(const string& str,int i){
	if((int)text.size()<=i){
		text.resize(i);
		lines=i;
	}
	text[i]+=(str);
	length=text[i].size();
}



bool MSA::check_sane()const{
	bool result(true);
	if(lines!=(int)text.size()){
		cout<<"Incorrect number of lines"<<endl;
		result=false;
	}
	for(int i(0);i<(int)text.size();++i){
		if((int)text[i].size()!=length){
			result=false;
			cout<<"Incorrect length"<<endl;

		}
		for(int j(0);j<(int)text[i].size();++j){
			if( alphabet.find(text[i][j])== std::string::npos){
				cout<<"Incorrect char"<<endl;
			}
		}
	}
	return result;
}



void MSA::print()const{
	check_sane()
	for(int i(0);i<(int)text.size();++i){
		cout<<text[i]<<endl;
	}
}


