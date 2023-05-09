//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Apr 16 09:54:16 2012 by ROOT version 5.32/00-rc2
// from TTree evetree/event tree
// found on file: /home/samuel/DAQDATA/Data/INORUN_150412_203036.ire
//////////////////////////////////////////////////////////

#ifndef main_h
#define main_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

// Header file for the classes stored in the TTree if any.
#include <TTimeStamp.h>
#include <TBits.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TGraphErrors.h>
#include <TBits.h>
#include "TProofOutputFile.h"

#define MAXMULT 3 // Maximum multiplicity reqd. for taking a layer into fitting, Default: 2 // DS: 8
#define MINLAY 4 // Minimum layers reqd for taking an event for fitting. Default: 4, do not decrease less than 3 // DS: 3
#define MAXSHIFT 2 // maximum deviation of the hit point from the fitted point
#define STRIPWIDTH 2.8 // cm
#define STRIPGAP 0.2 // cm
#define LAYERGAP 16 // cm
#define NL 12 // number of layers
#define NS 32 // number of strips
typedef struct {
    int multx[NL];
    int multy[NL];
}MULTIPLICITY;
typedef struct { // structure for timing
    double xtdc[NL];
    double ytdc[NL];
    int xtdcmult[NL];
    int ytdcmult[NL];
    double tdcref;
}TDC;
typedef struct { // structure for timing
    int xlmult;
    int ylmult;
}LAYERMULT;
typedef struct{
    float xgraph[NL]; // number of efficient hits
    float ygraph[NL];
}GRAPHPOINTS;
typedef struct { // structure for linear fitparameters
    float p0,p1,errp0,errp1,ndf,chisq;
}FITPARAM1;
typedef struct { // structure for poly fitparameters
    float p0,p1,p2,errp0,errp1,errp2,ndf,chisq;
}FITPARAM2;
typedef struct{
    float zen;
    float azi;
}ANGLE;

typedef struct{
    int xeff[12]; // number of efficient hits
    int yeff[12];
}EFFICIENCY;
typedef struct{
    int xstrip[12];
    int ystrip[12];
} STRIPNUM; // contains the strip num that was hit if multiplicity =1

// Fixed size dimensions of array or collections stored in the TTree if any.

class ino_analysis : public TSelector {
public :
    TTree          *fChain;   //!pointer to the analyzed TTree or TChain

    // Declaration of leaf types
    TTimeStamp      *Evetime;
    Int_t           ENum;
    TBits           *xstriphitsL0;
    TBits           *ystriphitsL0;
    TBits           *xstriphitsL1;
    TBits           *ystriphitsL1;
    TBits           *xstriphitsL2;
    TBits           *ystriphitsL2;
    TBits           *xstriphitsL3;
    TBits           *ystriphitsL3;
    TBits           *xstriphitsL4;
    TBits           *ystriphitsL4;
    TBits           *xstriphitsL5;
    TBits           *ystriphitsL5;
    TBits           *xstriphitsL6;
    TBits           *ystriphitsL6;
    TBits           *xstriphitsL7;
    TBits           *ystriphitsL7;
    TBits           *xstriphitsL8;
    TBits           *ystriphitsL8;
    TBits           *xstriphitsL9;
    TBits           *ystriphitsL9;
    TBits           *xstriphitsL10;
    TBits           *ystriphitsL10;
    TBits           *xstriphitsL11;
    TBits           *ystriphitsL11;
    ULong64_t       tdcdata[32];
    Int_t           tdcmult[32];
    ULong64_t       tdcref;

    // List of branches
    TBranch        *b_Evetime;   //!
    TBranch        *b_ENum;   //!
    TBranch        *b_xstriphitsL0;   //!
    TBranch        *b_ystriphitsL0;   //!
    TBranch        *b_xstriphitsL1;   //!
    TBranch        *b_ystriphitsL1;   //!
    TBranch        *b_xstriphitsL2;   //!
    TBranch        *b_ystriphitsL2;   //!
    TBranch        *b_xstriphitsL3;   //!
    TBranch        *b_ystriphitsL3;   //!
    TBranch        *b_xstriphitsL4;   //!
    TBranch        *b_ystriphitsL4;   //!
    TBranch        *b_xstriphitsL5;   //!
    TBranch        *b_ystriphitsL5;   //!
    TBranch        *b_xstriphitsL6;   //!
    TBranch        *b_ystriphitsL6;   //!
    TBranch        *b_xstriphitsL7;   //!
    TBranch        *b_ystriphitsL7;   //!
    TBranch        *b_xstriphitsL8;   //!
    TBranch        *b_ystriphitsL8;   //!
    TBranch        *b_xstriphitsL9;   //!
    TBranch        *b_ystriphitsL9;   //!
    TBranch        *b_xstriphitsL10;   //!
    TBranch        *b_ystriphitsL10;   //!
    TBranch        *b_xstriphitsL11;   //!
    TBranch        *b_ystriphitsL11;   //!
    TBranch        *b_TDCdata;   //!
    TBranch        *b_TDCmult;   //!
    TBranch        *b_TDCref;   //!

    ino_analysis(TTree * /*tree*/ =0) : fChain(0) { }
    virtual ~ino_analysis();
    // virtual Int_t   Version() const { return 2; }
    virtual void    Begin(TTree *tree);
    virtual void    SlaveBegin(TTree *tree);
    virtual void    Init(TTree *tree);
    virtual Bool_t  Notify();
    virtual Bool_t  Process(Long64_t entry);
    virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
    virtual void    SetOption(const char *option) { fOption = option; }
    virtual void    SetObject(TObject *obj) { fObject = obj; }
    virtual void    SetInputList(TList *input) { fInput = input; }
    virtual TList  *GetOutputList() const { return fOutput; }
    virtual void    SlaveTerminate();
    virtual void    Terminate();

    /*INO*/
    int  ii,jj,ll,BID,EvtIndx,xxindx,yyindx,fpts,pt,mm,kk;
    float xa,xb,ya,yb,za,zb,hyp,adj,opp;
    double layer;
    double stripNum,strip[MAXMULT];
    double strip0,strip11;

    // set the below flags to true to ignore those layers in fitting (for alignments)
    bool ignoreLayerX[NL]; // ignore layer in fitting
    bool ignoreLayerY[NL]; // ignore layer in fitting
    bool goodlayer[NL];     // true for layer with multiplicity <= 3 and hits in consecutive strips

    //Residuals
    bool useXResiduals;
    bool useYResiduals;
    float Xresidual[NL];
    float Yresidual[NL];
    TGraphErrors* fxgraph; // x hits
    TGraphErrors* fygraph;// y hits
    TGraphErrors* fitgraph; // fit points
    TGraphErrors* fitgraph2; // fit points for second fit (SDG)

    TGraphErrors* GetXGraph();
    TGraphErrors* GetYGraph();

    void SetFitError1(FITPARAM1& f,int errCode);
    void SetFitParameters1(FITPARAM1 &f,TF1* func);
    void FitData1(TGraphErrors *g,FITPARAM1& f, int mult[NL]); //get fit parameters for graph g and save in f using multiplicity
    void SetAngle1(FITPARAM1 &fx, FITPARAM1 &fy);
    void GetEfficiency1(FITPARAM1 &f, TGraphErrors *g,int mult[NL], int eff[NL]);

    void SetFitError2(FITPARAM2& f,int errCode);
    void SetFitParameters2(FITPARAM2 &f,TF1* func);
    void FitData2(TGraphErrors *g,FITPARAM2& f, int mult[NL]); //get fit parameters for graph g and save in f using multiplicity
    void GetEfficiency2(FITPARAM2 &f, TGraphErrors *g,int mult[NL], int eff[NL]);
    bool CheckTrackCoverage(FITPARAM1& f);  // returns true if fitted track is confined within the fiducial volume of Layer 0 and Layer 11

    TF1 *fitfunction;
    MULTIPLICITY multiplicity;
    TDC tdc;
    LAYERMULT layermult;
    FITPARAM1 xfitparam1;
    FITPARAM1 yfitparam1;
    FITPARAM2 xfitparam2;
    FITPARAM2 yfitparam2;
    GRAPHPOINTS graphpts;
    ANGLE angle;
    EFFICIENCY eff;
    STRIPNUM stripnum;
    TBits *xstriphits[NL];
    TBits *ystriphits[NL];
    TString outs;


    //output is stored here
    TProofOutputFile *fProofFile;
    TFile *fFile;
    TTree* fittree;
    TString outFileName;
    //
    // for summary
    TH1F *summEff[12];

    TString option;

};




#endif