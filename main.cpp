// g++ -o analysis ino_analysis.cpp main.cpp `root-config --glibs --cflags`
//./analysis.out
#include <TTree.h>
#include <TFile.h>
#include <iostream>
#include "ino_analysis.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

using namespace TMVA;
using namespace ::std;
int main()
{
    TString input = "INORUN_20160819_053749.ire"; // input file
    // output will be stored as out.fit.// This name can be changed in ino_analysis.cpp

    std::unique_ptr<TFile> file(TFile::Open(input));
    if (!file || file->IsZombie())
    {
        std::cerr << "Error opening file" << endl;
        exit(-1);
    }

 
    TTree *T = (TTree *)file->Get("evetree");
    int N = T->GetEntries();
    cout << N << " entries in file " << input << endl;

    ino_analysis *a = new ino_analysis(T);
    a->Init(T);
    a->SlaveBegin(T);

    for (int i = 0; i < N; i++)
    {
        a->Process(i);
        if (i % 5000 == 0)
            cout << (float)i/N *100 << "% entries processed..." << endl;
    }

    cout << "Processing completed..." << endl;
    a->SlaveTerminate();

    return 0;
}