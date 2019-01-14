//Output file
#ifndef NETMC_OUTPUTFILE_H
#define NETMC_OUTPUTFILE_H

#include <iostream>
#include <fstream>
#include <chrono>

using namespace std;

//Generic output file
class OutputFile{

public:

    //Data members
    ofstream file;
    chrono::high_resolution_clock::time_point t0,t1;

    //Make non-copyable
    OutputFile(const OutputFile&);
    OutputFile& operator=(const OutputFile&);

    //Constructor
    OutputFile(): file("a.out",ios::in|ios::trunc){}
    OutputFile(string name): file(name,ios::in|ios::trunc){}

    //Destructor
    ~OutputFile(){file.close();}

    //Member functions
    void datetime(string message=""){
        time_t now=time(0);
        char* dt=ctime(&now);
        string time(dt);
        file<<message<<time;
    }

    void separator(int len=70){
        string sep;
        for(int i=0; i<len; ++i) sep+="-";
        file<<sep<<endl;
    }

    void write(string message){
        file<<message<<endl;
    }


};

//Logfile
class Logfile: public OutputFile{

public:

    //Constructor
    Logfile(string name):OutputFile(name){};

    //Member functions
    void criticalError(string message){
        file<<"Critical error: "<< message<<endl;
        throw(message);
    }

};

#endif //NETMC_OUTPUTFILE_H
