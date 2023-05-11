#include "convert_ascii.h"
#include <TH2.h>
#include <TStyle.h>
#include <TH2.h>
#include <TF1.h>
#include <TMath.h>
#include <TGraph.h>
#include <TDirectory.h>
#include <TProof.h>


#define SIGMASTRIP 0.28 // unit of strips  = 2.8/sqrt(12) = 0.8 assuming uniform distribution [0.8/2.8=0.28]
#define SIGMATIME 1.5 //ns
#define SIGMALAYER 0.01 //units of layers

void convert_ascii::Begin(TTree * /*tree*/)
{
    // The Begin() function is called at the start of the query.
    // When running with PROOF Begin() is only called on the client.
    // The tree argument is deprecated (on PROOF 0 is passed).

    TString option = GetOption();

}

void convert_ascii::SlaveBegin(TTree * /*tree*/)
{
    // The SlaveBegin() function is called after the Begin() function.
    // When running with PROOF SlaveBegin() is called on each slave server.
    // The tree argument is deprecated (on PROOF 0 is passed).
    // TFile *fInput=new TFile("INORUN_20161230_161712.ire");
    // option = GetOption();
    fxgraph= new TGraphErrors();
    fygraph= new TGraphErrors();
    fitgraph = new TGraphErrors();
   //residual = xstrip-fitfunc
    useXResiduals=false;
    useYResiduals=false;

     for(int l=0;l<NL;l++){
        ignoreLayerX[l]=false;
        ignoreLayerY[l]=false;
        Xresidual[l]=0.0;
        Yresidual[l]=0.0;
    }

    TString xname,yname;

    for(int ii=0;ii<NL;ii++){
        xstriphits[ii]= new TBits(NS);
        ystriphits[ii]= new TBits(NS);
        xname.Form("xstriphits_%d",ii);
        yname.Form("ystriphits_%d",ii);
        // fittree->Branch(xname,"TBits",&xstriphits[ii],10000,0);
        // fittree->Branch(yname,"TBits",&ystriphits[ii],10000,0);
    }
   
   outfile.open(outFileName);

}

Bool_t convert_ascii::Process(Long64_t entry)
{
    // The Process() function is called for each entry in the tree (or possibly
    // keyed object in the case of PROOF) to be processed. The entry argument
    // specifies which entry in the currently loaded tree is to be processed.
    // It can be passed to either convert_ascii.C::GetEntry() or TBranch::GetEntry()
    // to read either all or the required parts of the data. When processing
    // keyed objects with PROOF, the object is already loaded and is available
    // via the fObject pointer.
    //
    // This function should contain the "body" of the analysis. It can contain
    // simple or elaborate selection criteria, run algorithms on the data
    // of the event and typically fill histograms.
    //
    // The processing can be stopped by calling Abort().
    //
    // Use fStatus to set the return value of TTree::Process().
    //
    // The return value is currently not used.

    // copy values from the original tree
    xstriphits[0] = xstriphitsL0;
    ystriphits[0] = ystriphitsL0;
    xstriphits[1] = xstriphitsL1;
    ystriphits[1] = ystriphitsL1;
    xstriphits[2] = xstriphitsL2;
    ystriphits[2] = ystriphitsL2;
    xstriphits[3] = xstriphitsL3;
    ystriphits[3] = ystriphitsL3;
    xstriphits[4] = xstriphitsL4;
    ystriphits[4] = ystriphitsL4;
    xstriphits[5] = xstriphitsL5;
    ystriphits[5] = ystriphitsL5;
    xstriphits[6] = xstriphitsL6;
    ystriphits[6] = ystriphitsL6;
    xstriphits[7] = xstriphitsL7;
    ystriphits[7] = ystriphitsL7;
    xstriphits[8] = xstriphitsL8;
    ystriphits[8] = ystriphitsL8;
    xstriphits[9] = xstriphitsL9;
    ystriphits[9] = ystriphitsL9;
    xstriphits[10] = xstriphitsL10;
    ystriphits[10] = ystriphitsL10;
    xstriphits[11] = xstriphitsL11;
    ystriphits[11] = ystriphitsL11;
    GetEntry(entry);

    // GetGraphs
    outfile <<ENum<<","<<Evetime->GetSec()<<",";
    GetXGraph();
    GetYGraph();
    outfile <<" "<<std::endl;
    
   
    return kTRUE;

}

TGraphErrors* convert_ascii::GetXGraph(){
    //TGraphErrors* fxgraph= new TGraphErrors();
    ii=0; jj=0;
    fxgraph->Set(0);
    xxindx=0;
    for(ii=0;ii<NL;ii++)
    {
        if(ignoreLayerX[ii]==false){ // set point only if layer is not ignored
            for(jj=0;jj<NS;jj++){
                if(xstriphits[ii]->TestBitNumber(jj)){
                    fxgraph->SetPoint(xxindx,ii,jj-Xresidual[ii]);//jj,BID);
                    fxgraph->SetPointError(xxindx,SIGMALAYER,SIGMASTRIP);
                    xxindx++;
                    outfile << getStripID(ii,jj) <<",";
                }
            }
        }

    }
    fxgraph->Set(xxindx);
    fxgraph->Sort(); // sort along layer
    return fxgraph;
}

TGraphErrors* convert_ascii::GetYGraph(){

    ii=0; jj=0;
    fygraph->Set(0);
    yyindx=0;
    for(ii=0;ii<NL;ii++)
    {
        if(ignoreLayerX[ii]==false){ // set point only if layer is not ignored
            for(jj=0;jj<NS;jj++){
                if(ystriphits[ii]->TestBitNumber(jj)){
                    fygraph->SetPoint(yyindx,ii,jj-Yresidual[ii]);//jj,BID);
                    fygraph->SetPointError(yyindx,SIGMALAYER,SIGMASTRIP);
                    yyindx++;
                    outfile << getStripID(ii,jj,false) <<",";
                }
            }
        }

    }
    fygraph->Set(yyindx);
    fygraph->Sort(); // sort along layer
    return fygraph;
}

void convert_ascii::SlaveTerminate()
{
    // The SlaveTerminate() function is called after all entries or objects
    // have been processed. When running with PROOF SlaveTerminate() is called
    // on each slave server.

outfile.close();
  
}

void convert_ascii::Terminate()
{
    // The Terminate() function is the last function to be called during
    // a query. It always runs on the client, it can be used to present
    // the results graphically or save the results to file.

}
convert_ascii::~convert_ascii()
{

}


void convert_ascii::Init(TTree *tree)
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    // Set object pointer
    Evetime = 0;
    xstriphitsL0 = 0;
    ystriphitsL0 = 0;
    xstriphitsL1 = 0;
    ystriphitsL1 = 0;
    xstriphitsL2 = 0;
    ystriphitsL2 = 0;
    xstriphitsL3 = 0;
    ystriphitsL3 = 0;
    xstriphitsL4 = 0;
    ystriphitsL4 = 0;
    xstriphitsL5 = 0;
    ystriphitsL5 = 0;
    xstriphitsL6 = 0;
    ystriphitsL6 = 0;
    xstriphitsL7 = 0;
    ystriphitsL7 = 0;
    xstriphitsL8 = 0;
    ystriphitsL8 = 0;
    xstriphitsL9 = 0;
    ystriphitsL9 = 0;
    xstriphitsL10 = 0;
    ystriphitsL10 = 0;
    xstriphitsL11 = 0;
    ystriphitsL11 = 0;
    // Set branch addresses and branch pointers
    if (!tree) return;
    fChain = tree;
    fChain->SetMakeClass(1);

    fChain->SetBranchAddress("Evetime", &Evetime, &b_Evetime);
    fChain->SetBranchAddress("ENum", &ENum, &b_ENum);
    fChain->SetBranchAddress("xstriphitsL0", &xstriphitsL0, &b_xstriphitsL0);
    fChain->SetBranchAddress("ystriphitsL0", &ystriphitsL0, &b_ystriphitsL0);
    fChain->SetBranchAddress("xstriphitsL1", &xstriphitsL1, &b_xstriphitsL1);
    fChain->SetBranchAddress("ystriphitsL1", &ystriphitsL1, &b_ystriphitsL1);
    fChain->SetBranchAddress("xstriphitsL2", &xstriphitsL2, &b_xstriphitsL2);
    fChain->SetBranchAddress("ystriphitsL2", &ystriphitsL2, &b_ystriphitsL2);
    fChain->SetBranchAddress("xstriphitsL3", &xstriphitsL3, &b_xstriphitsL3);
    fChain->SetBranchAddress("ystriphitsL3", &ystriphitsL3, &b_ystriphitsL3);
    fChain->SetBranchAddress("xstriphitsL4", &xstriphitsL4, &b_xstriphitsL4);
    fChain->SetBranchAddress("ystriphitsL4", &ystriphitsL4, &b_ystriphitsL4);
    fChain->SetBranchAddress("xstriphitsL5", &xstriphitsL5, &b_xstriphitsL5);
    fChain->SetBranchAddress("ystriphitsL5", &ystriphitsL5, &b_ystriphitsL5);
    fChain->SetBranchAddress("xstriphitsL6", &xstriphitsL6, &b_xstriphitsL6);
    fChain->SetBranchAddress("ystriphitsL6", &ystriphitsL6, &b_ystriphitsL6);
    fChain->SetBranchAddress("xstriphitsL7", &xstriphitsL7, &b_xstriphitsL7);
    fChain->SetBranchAddress("ystriphitsL7", &ystriphitsL7, &b_ystriphitsL7);
    fChain->SetBranchAddress("xstriphitsL8", &xstriphitsL8, &b_xstriphitsL8);
    fChain->SetBranchAddress("ystriphitsL8", &ystriphitsL8, &b_ystriphitsL8);
    fChain->SetBranchAddress("xstriphitsL9", &xstriphitsL9, &b_xstriphitsL9);
    fChain->SetBranchAddress("ystriphitsL9", &ystriphitsL9, &b_ystriphitsL9);
    fChain->SetBranchAddress("xstriphitsL10", &xstriphitsL10, &b_xstriphitsL10);
    fChain->SetBranchAddress("ystriphitsL10", &ystriphitsL10, &b_ystriphitsL10);
    fChain->SetBranchAddress("xstriphitsL11", &xstriphitsL11, &b_xstriphitsL11);
    fChain->SetBranchAddress("ystriphitsL11", &ystriphitsL11, &b_ystriphitsL11);
    fChain->SetBranchAddress("tdcdata", tdcdata, &b_TDCdata);
    fChain->SetBranchAddress("tdcmult", tdcmult, &b_TDCmult);
    fChain->SetBranchAddress("tdcref", &tdcref, &b_TDCref);
}

Bool_t convert_ascii::Notify()
{
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.

    return kTRUE;
}

// a function which returns a unique id for a strip 
// x side layer 0 strip id will be from 0 -31
// x side layer 1 strip id will be from 32 -63 
// max strip id on xside will be 383 for 12 layers and 32 strips
// y side layer 0 strip id will be from 384 -415
// y side layer 1 strip id will be from 416 -447 

int convert_ascii::getStripID(int layer, int strip, bool isXside)
{
    if(isXside)
    return (layer*NS) + strip;
    else
    return (NS*NL) + (layer*NS) + strip;
}