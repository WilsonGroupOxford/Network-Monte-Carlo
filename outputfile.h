//Output file
#ifndef NETMC_OUTPUTFILE_H
#define NETMC_OUTPUTFILE_H

#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>

using namespace std;

//Generic output file
class OutputFile{

public:

    //Data members
    ofstream file;
    int spacing, currIndent;
    string indent,dashed;

    //Make non-copyable
    OutputFile(const OutputFile&);
    OutputFile& operator=(const OutputFile&);

    //Constructor
    OutputFile(): file("a.out",ios::in|ios::trunc){
        initVariables();
    }
    explicit OutputFile(const string name): file(name,ios::in|ios::trunc){
        initVariables();
    }

    //Destructor
    ~OutputFile(){file.close();}

    //Initialise file variables
    void initVariables(int precision=8, int indentSize=4, int sepSize=60, int spaceSize=20){
        file << fixed << showpoint << setprecision(precision);
        indent="";
        dashed="";
        for(int i=0; i<indentSize; ++i) indent+=" ";
        for(int i=0; i<sepSize; ++i) dashed+="-";
        spacing=spaceSize;
        currIndent=0;
    }

    //Member functions
    void datetime(string message=""){
        time_t now=time(0);
        char* dt=ctime(&now);
        string time(dt);
        file<<message<<time;
    }

    void separator(){
        file<<dashed<<endl;
    }

    template <typename T>
    void write(T val){
        for(int i=0; i<currIndent; ++i) file<<indent;
        file<<val<<endl;
    }

    template <typename T, typename U>
    void write(T val0, U val1){
        for(int i=0; i<currIndent; ++i) file<<indent;
        file<<val0<<" "<<val1<<endl;
    }

    template <typename T>
    void writeRowVector(T vec){
        for(int i=0; i<vec.n; ++i){
            file<<setw(spacing)<<left<<vec[i];
        }
        file<<endl;
    }

};

//Logfile
class Logfile: public OutputFile{

public:

    //Data members
    chrono::high_resolution_clock::time_point tStart,tEnd,t0,t1;

    //Constructor
    explicit Logfile(const string name):OutputFile(name){
        tStart=chrono::high_resolution_clock::now();
        t0=tStart;
    };

    //Destructor
    ~Logfile(){
        tEnd=chrono::high_resolution_clock::now();
        double dt=chrono::duration_cast<chrono::seconds>(tEnd-tStart).count();
        dt/=60.0;
        file<<"Duration: "<<dt<<" minutes"<<endl;
    }

    //Member functions
    void criticalError(string message){
        file<<"Critical error: "<< message<<endl;
        throw(message);
    }

    double timeElapsed(){
        t1=chrono::high_resolution_clock::now();
        double dt=chrono::duration_cast<chrono::milliseconds>(t1-t0).count();
        t0=t1;
        return dt/1000.0;
    }

};

#endif //NETMC_OUTPUTFILE_H
