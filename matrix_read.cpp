#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <vector>
#include <map>
#include <cstdlib>
#include <boost/filesystem.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <tclap/CmdLine.h>
//#include <boost/algorithm/string_regex.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/lexical_cast.hpp>
// cc -std=c++17 matrix_read.cpp -o matrix_read -I ../software/169boost/include/ -I /home/kris/software/tclap/include -L ../software/169boost/lib/ -lstdc++ -lboost_iostreams -lboost_filesystem

using namespace std;
namespace fs = boost::filesystem;
namespace al = boost::algorithm;
namespace io = boost::iostreams;
/*
 barcode : row
 feature : col
barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz
*/
void unzip(string path,stringstream& ss){
  ifstream file(path);
  io::filtering_streambuf<io::input> in;
  in.push(io::gzip_decompressor());
  in.push(file);
  io::copy(in,ss);
  /*string line;
  getline(ss,line);
  cout << line;
  */
}

void zip(string path,stringstream& ss){
  ofstream file(path);
  io::filtering_streambuf<io::output> out;
  out.push(io::gzip_compressor());
  out.push(file);
  ostream os(&out);
  os << ss.str();
}

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

string prematrix(stringstream& ss){
  string line;
  getline(ss,line);
  while(boost::starts_with(line,"%")){
    line.clear();
    getline(ss,line);
  }
  return line;
}

class dictfile{
    private:
        map<string,pair<int,string> > gene;
        ifstream file;
    public:
        dictfile(string filepath) :  file(filepath) {}
        void proc(){
            string line,tmp,tmp1,name,key;
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
                ss.str("");
                ss.clear();
                ss << container[4];
                ss >> end;
                ss.str("");
                ss.clear();
                tmp = split(split(container[8],";")[0]," ")[1];
                tmp1 = split(split(container[8],";")[5]," ")[2];
                //cout << key.substr(1,key.length()-2);
                key = tmp.substr(1,tmp.length()-2);
                name = tmp1.substr(1,tmp1.length()-2);
                if(gene.count(key) == 0){
                    gene.insert({key,{end-start,name}});
                }
                else{
                    gene[key].first = gene[key].first + end - start;
                }
                line.clear();
                tmp1.clear();
                tmp.clear();
                key.clear();
                name.clear();
            }
            file.close();
            cerr << "end of processing dictionary .." << endl;
        }
        void show(){
            for(auto& i : gene){
                cout << i.first << ":" << i.second.first << "symbol : " << i.second.second << endl;
            }
        }
        map<string,pair<int,string> > getdict(){
            return gene;
        }
};

class processor{
  private:
    vector<string> barcodes;
    vector<string> features;
    double *tpm;
    int *mtx;
    int sR,sC;
  public:
    processor(stringstream& barcode,stringstream& feature,stringstream& matrix,int R,int C){
      sR = R;
      sC = C;
      tpm = new double[R * C];
      mtx = new int[R * C];
      string line;
      vector<string> tmp;
      while(getline(barcode,line)){
        barcodes.push_back(line);
        line.clear();
      }
      line.clear();
      while (getline(feature,line)) {
        tmp = split(line,"\t");
        features.push_back(tmp[0]);
        line.clear();
        tmp.clear();
      }
      line.clear();
      tmp.clear();
      //int re = 0;
      while(getline(matrix,line)){
        tmp = split(line," ");
        int index = (boost::lexical_cast<int>(tmp[1]) - 1) * C + boost::lexical_cast<int>(tmp[0]) - 1;
        int value = boost::lexical_cast<int>(tmp[2]);
        //cout << "mtx[" << index << "]="<<value<<endl;
        //if(mtx[index] != 0){
        //  re++;
        //}
        mtx[index] = value;
        line.clear();
        tmp.clear();
      }
      //cout << "re" << re <<endl;
    }
    ~processor(){
      delete [] mtx;
      delete [] tpm;
    }

    void TPM(map<string,pair<int,string> > dict){
        cout << "Start counting TPM.." << endl;
        for(int i=0;i<sR;i++){
          double inter,tmp;
          double total = 0.0;
          for(int j=0;j<sC;j++){
            if(dict.count(feature(j)) == 0){
              cerr << "no feature for " << feature(j) << endl;
            }
            inter = dict[feature(j)].first / 1000.0;
            tmp = mtx[i * sC + j] / inter;
            total += tmp;
            tpm[i * sC + j] = tmp;
          }
          total = total / 1000000.0;
          for(int k=0;k<sC;k++){
            if(dict.count(feature(k)) == 0){
              cerr << "no feature for " << feature(k) << endl;
            }
            tpm[i * sC + k] = tpm[i * sC + k] / total;
          }
        }
        cout << "End of counting TPM.." << endl;
    }

    void save(string name,bool gswitch,map<string,pair<int,string> > dict){
      if(!fs::exists("tpmoutput")){
        fs::create_directory("tpmoutput");
      }
      fs::current_path("tpmoutput");
      string mdname;
      mdname = name + ".features.gz";
      stringstream ss;
      for(int i=0;i<features.size();i++){
        ss << features[i];
        if(i != (features.size()-1)){
          ss << "\t";
        }
        else{
          ss << "\n";
        }
      }
      zip(mdname,ss);
      mdname.clear();
      mdname = name + ".mtx.gz";
      ss.clear();
      ss.str("");
      for(int i=0;i<sR;i++){
        for(int j=0;j<sC;j++){
          ss << mtx_value(i,j);
          if(j != (sC-1)){
            ss << "\t";
          }
          else{
            ss << "\n";
          }
        }
      }
      zip(mdname,ss);
      mdname.clear();
      mdname = name + ".tpm.gz";
      ss.str("");
      ss.clear();
      for(int i=0;i<sR;i++){
        for(int j=0;j<sC;j++){
          ss << fixed;
          ss << setprecision(3);
          ss << tpm_value(i,j);
          if(j != (sC-1)){
            ss << "\t";
          }
          else{
            ss << "\n";
          }
        }
      }
      zip(mdname,ss);
      if(gswitch){
        mdname.clear();
        mdname = name + ".gs.gz";
        ss.str("");
        ss.clear();
        for(int i=0;i<features.size();i++){
          ss << dict[features[i]].second;
          if(i != (features.size()-1)){
            ss << "\t";
          }
          else{
            ss << "\n";
          }
        }
        zip(mdname,ss);
      }
    }

    string feature(int index){
      return features[index];
    }
    string barcode(int index){
      return barcodes[index];
    }
    int mtx_value(int row,int col){
      return mtx[row * sC + col];
    }
    double tpm_value(int row,int col){
      return tpm[row * sC + col];
    }
    void show(){
      int mnz = 0;
      int tnz = 0;
      for(int i=0;i<sR*sC;i++){
        if(mtx[i] != 0){
          mnz++;
        }
        if(tpm[i] != 0){
          tnz++;
        }
      }
      cout << "barcodes :" << barcodes.size() << endl;
      cout << "features :" << features.size() << endl;
      cout << "mtx all :" << sR * sC << endl;
      cout << "mtx non-zero :" << mnz << endl;
      cout << "tpm all :" << sR * sC << endl;
      cout << "tpm non-zero :" << tnz << endl;
    }

};

/*
decompressing
ifstream file("hello");
io::filtering_streambuf<io::input> in;
in.push(io::gzip_decompressor());
in.push(file);
io::copy(in,cout);
*/

void flow(string path,string dictpath,string name,bool gswitch = false){
  fs::path cwd;
  cwd = fs::current_path();
  fs::current_path(path);
  stringstream bs,fs,ms;
  if(!fs::exists("barcodes.tsv.gz")){
    cerr << "barcodes.tsv.gz not exists." <<endl;
    abort();
  }
  if(!fs::exists("features.tsv.gz")){
    cerr << "features.tsv.gz not exists." <<endl;
    abort();
  }
  if(!fs::exists("matrix.mtx.gz")){
    cerr << "matrix.mtx.gz not exists." <<endl;
    abort();
  }
  unzip("barcodes.tsv.gz",bs);
  unzip("features.tsv.gz",fs);
  unzip("matrix.mtx.gz",ms);
  string header;
  header = prematrix(ms);
  vector<string> v;
  v = split(header," ");
  //cout << v[0] << " " <<v[1];
  int C = boost::lexical_cast<int>(v[0]);
  int R = boost::lexical_cast<int>(v[1]);
  dictfile d(dictpath);
  d.proc();
  processor p(bs,fs,ms,R,C);
  p.show();
  p.TPM(d.getdict());
  fs::current_path(cwd);
  p.save(name,gswitch,d.getdict());
  p.show();
}

int main(int argc, char const *argv[]) {
   ///home/kris/single-cell-data/output/case1/c1/outs/filtered_feature_bc_matrix
   //33538 2489 2412326
   string file,gtf,name;
   bool gswitch;
   try{
     TCLAP::CmdLine cmd("A program for counting TPM, input a folder from 10xGenomics containing filter_matrix", ' ', "1");
     TCLAP::ValueArg<string> pathArg("f","folder","Folder containing matrix, barcode, features.",true,".","string");
     TCLAP::ValueArg<string> gtfArg("g","gtf","The gtf file for",true," ","string");
     TCLAP::ValueArg<string> nameArg("n","name","The prefix for the output file",true,".","string");
     TCLAP::SwitchArg gsArg("s","symbol","Output gene symbol",cmd,false);
     cmd.add(pathArg);
     cmd.add(gtfArg);
     cmd.add(nameArg);
     cmd.parse(argc,argv);
     gswitch = gsArg.getValue();
     file = pathArg.getValue();
     gtf = gtfArg.getValue();
     name = nameArg.getValue();
   }catch(TCLAP::ArgException &e){
     cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
   }
   flow(file,gtf,name,gswitch);
   return 0;
}
