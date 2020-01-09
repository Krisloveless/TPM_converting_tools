#include <iostream>
#include <string>
#include <map>
#include <fstream>
#include <vector>
#include <sstream>
#include <iomanip>
#include <filesystem>
#include <tclap/CmdLine.h>
#include <thread>
#include <mutex>


//cc TpmC.cpp -std=c++17 -o TpmC -lstdc++ -lstdc++fs -lpthread -I ../software/tclap/include/
using namespace std;
namespace fs = std::filesystem;
mutex mtx;

vector<string> split(string origin,string reg){
    vector<string> res;
    size_t start = 0;
    size_t found = origin.find(reg);
    string s;
    while(found != string::npos){
        s = origin.substr(start,found-start);
        res.push_back(s);
        s.clear();
        start = found+reg.length();
        found = origin.find(reg,start);
    }
    //npos indicate until the end of the line : -1
    res.push_back(origin.substr(start,string::npos));
    return res;
}

class dictfile{
    private:
        map<string,int> gene;
        ifstream file;
    public:
        dictfile(string filepath) :  file(filepath) {}
        void proc(){
            string line,tmp,key;
            int start,end;
            vector<string> container;
            stringstream ss;
            cerr << "processing dictionary .." << endl;
            while(getline(file,line)){
                if(line.substr(0,2) == "#!"){
                    continue;
                }
                container = split(line,"\t");
                //"ENSG00000243485"
                //cout << split(split(container[8],";")[0]," ")[1] <<endl;
                if(container[2] != "exon"){
                    continue;
                }
                ss << container[3];
                ss >> start;
                ss.clear();
                ss << container[4];
                ss >> end;
                ss.clear();
                tmp = split(split(container[8],";")[0]," ")[1];
                //cout << key.substr(1,key.length()-2);
                key = tmp.substr(1,tmp.length()-2);
                if(gene.count(key) == 0){
                    gene.insert({key,end-start});
                }
                else{
                    gene[key] = gene[key] + end - start;
                }
                line.clear();
                tmp.clear();
                key.clear();
            }
            file.close();
            cerr << "end of processing dictionary .." << endl;
        }
        void show(){
            for(auto& i : gene){
                cout << i.first << ":" << i.second << endl;
            }
        }
        map<string,int> getdict(){
            return gene;
        }

};

class counting{
    private:
        ifstream ifile;
        ofstream ofile;
        string name;
    public:
        counting(string filepath,string output) : ifile(filepath),ofile(output),name(output) {}
        void proc(map<string,int> dict){
            double total = 0.0;
            double tmp,inter,value;
            string line,key;
            stringstream ss;
            cerr << "Starting TPM for " << name << endl;
            while(getline(ifile,line)){
                if(line.substr(0,2) == "__"){
                    continue;
                }
                key = split(line,"\t")[0];
                if(dict.count(key) == 0){
                    cerr << key <<"not found" <<endl;
                    return;
                }
                inter = dict[key] / 1000.0;
                ss << split(line,"\t")[1];
                ss >> tmp;
                tmp = tmp / inter;
                ss.clear();
                //cout << tmp << endl;
                total += tmp;
                line.clear();
            }
            total = total / 1000000.0;
            //cout << "total :" << total << endl;
            ifile.clear();
            ifile.seekg(0,ios::beg);
            //cerr << "end of calculating total length .. for " << name << endl;
            //cerr << "calculating TPM .. for " << name << endl;
            while(getline(ifile,line)){
                if(line.substr(0,2) == "__"){
                    continue;
                }
                ss << split(line,"\t")[1];
                ss >> value;
                ss.clear();
                key = split(line,"\t")[0];
                if(dict.count(key) == 0){
                    cerr << key <<"not found" <<endl;
                    return;
                }
                inter = dict[key]/1000;
                value = value / inter;
                value = value / total;
                ofile << fixed;
                ofile << setprecision(3);
                ofile << key << "\t" << value << endl;
                key.clear();
                line.clear();
            }
            ifile.close();
            ofile.close();
            cerr << "end of calculating TPM for "<< name << endl;
        }
};

void process(vector<fs::path> &file, dictfile &d){
  fs::path p;
  stringstream name;
  string out = "";
  while(1){
    mtx.lock();
    if(file.size() == 0){
      mtx.unlock();
      return;
    }
    p = file.back();
    //idiot question
    file.pop_back();
    //cerr << fs::current_path() << endl;
    out += "TPMout/";
    name << p.filename();
    out += name.str().substr(1,name.str().length()-2);
    out += ".out";
    counting c(p,out);
    name.clear();
    out.clear();
    name.str("");
    mtx.unlock();
    c.proc(d.getdict());

  }
}

void x(int c){
  return;
}

int main(int argc,char *argv[]){
  string file,gtf;
  int ts;
  try{
    TCLAP::CmdLine cmd("A program for counting TPM, input a folder containing htseq-count output", ' ', "1");
    TCLAP::ValueArg<string> nameArg("f","folder","Folder containing htseq-count output",true,".","string");
    TCLAP::ValueArg<int> threadArg("t","thread","How many threads to run",false,1,"int");
    TCLAP::ValueArg<string> gtfArg("g","gtf","The gtf file for",true," ","string");
    cmd.add(nameArg);
    cmd.add(gtfArg);
    cmd.add(threadArg);
    cmd.parse(argc,argv);
    file = nameArg.getValue();
    ts = threadArg.getValue();
    gtf = gtfArg.getValue();
    //cout << "thread : " << threads << ", file : " << file << endl;
  }catch(TCLAP::ArgException &e){
    cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
  }
  vector<fs::path> filenames;
  for(auto &p : fs::directory_iterator(file)){
    filenames.push_back(p.path());
  }
  dictfile d(gtf);
  d.proc();
  fs::create_directory("TPMout");
  thread *cores[ts];
  //thread x(process,ref(filenames),ref(d));

  for(int i=0;i<ts;i++){
    cores[i] = new thread(process,ref(filenames),ref(d));
  }
  for(int i=0;i<ts;i++){
    cores[i]->join();
  }

  for(int i=0;i<ts;i++){
    delete cores[i];
  }


  /*thread a(x,10);
  a.join();
  */
  return 0;

  /*
    if(argc != 4){
        cerr << "input error : please input a gene.gtf, file.count and a name" << endl;
        return 1;
    }
    dictfile d(argv[1]);
    d.proc();
    counting c(argv[2],argv[3]);
    c.proc(d.getdict());
    return 0;
    */
}
