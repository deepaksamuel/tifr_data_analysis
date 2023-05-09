#include "ino_analysis.h"
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

void ino_analysis::Begin(TTree * /*tree*/)
{
    // The Begin() function is called at the start of the query.
    // When running with PROOF Begin() is only called on the client.
    // The tree argument is deprecated (on PROOF 0 is passed).

    TString option = GetOption();

}

void ino_analysis::SlaveBegin(TTree * /*tree*/)
{
    // The SlaveBegin() function is called after the Begin() function.
    // When running with PROOF SlaveBegin() is called on each slave server.
    // The tree argument is deprecated (on PROOF 0 is passed).
    // TFile *fInput=new TFile("INORUN_20161230_161712.ire");
    // option = GetOption();
    fxgraph= new TGraphErrors();
    fygraph= new TGraphErrors();
    fitgraph = new TGraphErrors();


     for(int l=0;l<NL;l++){
        ignoreLayerX[l]=false;
        ignoreLayerY[l]=false;
        Xresidual[l]=0.0;
        Yresidual[l]=0.0;
    }
	//for even layers only: CUK changes
	//  ignoreLayerX[0]=true;
	//  ignoreLayerX[2]=true;
    //  ignoreLayerX[4]=true;
    //  ignoreLayerX[6]=true;
	//  ignoreLayerX[8]=true;
	//  ignoreLayerX[9]=true;
    //  ignoreLayerX[10]=true;
    //  ignoreLayerX[11]=true;
    //  ignoreLayerX[1]=true;
    //  ignoreLayerX[1]=true;
	
    //residual = xstrip-fitfunc
    useXResiduals=false;
    useYResiduals=false;

    // use residuals this way
    //    if(useXResiduals==true){
    //        Xresidual[0]=0.253534;
    //    }
    //    if(useYResiduals==true){
    //        Yresidual[0]=0.185485;
    //    }
    // ignore layers this way
    //    ignoreLayerX[2]=true;

    // Example of adding graphs and histograms in output list
    //histo = new TH1F("Proof","I am Proof",1000,0,60);
    //graph = new TGraph();
    //graph->SetName("graph");
    //GetOutputList()->Add(histo);
    //GetOutputList()->Add(graph);
    
    //fInput="INORUN_20160918_205857.ire"
    // option = GetOption();
    option = "NOMAG";

    if(option==TString("MAG")){ // magnetic field present
        // Warning("SlaveBegin","Option: MAG - Using Second order fitting");
        fitfunction=new TF1("Poly_Fit","pol2", 0,0);
    } else{ // no magnetic field present
        // Warning("SlaveBegin","Option: None - Using First order fitting");
        fitfunction=new TF1("Lin_Fit","pol1", 0,0);
    }

    TFile *outFile = new TFile("out.fit","RECREATE");
    fFile = outFile;
    fittree = new TTree("analysis_results","Analysis Results");  //GetOutputList()->Add(fittree);


    fittree->Branch("ENum",&ENum,"ENum/i");
    fittree->Branch("layer_multiplicity",&layermult,"xlmult/I:ylmult/I");
    fittree->Branch("tdc",&tdc,"xtdc[12]/D:ytdc[12]/D:xtdcmult[12]/i:ytdcmult[12]/i:tdcref/D");
    fittree->Branch("time","TTimeStamp",&Evetime,10000,0);
    if(option==TString("MAG")){
        fittree->Branch("X_Fit",&xfitparam2,"xp0:xp1:xp2:xerrp0:xerrp1:xerrp2:xndf:xchisq");
        fittree->Branch("Y_Fit",&yfitparam2,"yp0:yp1:yp2:yerrp0:yerrp1:yerrp2:yndf:ychisq");
    }
    else{
        fittree->Branch("X_Fit",&xfitparam1,"xp0:xp1:xerrp0:xerrp1:xndf:xchisq");
        fittree->Branch("Y_Fit",&yfitparam1,"yp0:yp1:yerrp0:yerrp1:yndf:ychisq");
    }
    fittree->Branch("angle",&angle,"zen/f:azi/f");
    fittree->Branch("efficiency",&eff,"xeff[12]/i:yeff[12]/i");
    TString xname,yname;

    for(int ii=0;ii<NL;ii++){
        xstriphits[ii]= new TBits(NS);
        ystriphits[ii]= new TBits(NS);
        xname.Form("xstriphits_%d",ii);
        yname.Form("ystriphits_%d",ii);
        // fittree->Branch(xname,"TBits",&xstriphits[ii],10000,0);
        // fittree->Branch(yname,"TBits",&ystriphits[ii],10000,0);
    }


}

Bool_t ino_analysis::Process(Long64_t entry)
{
    // The Process() function is called for each entry in the tree (or possibly
    // keyed object in the case of PROOF) to be processed. The entry argument
    // specifies which entry in the currently loaded tree is to be processed.
    // It can be passed to either ino_analysis.C::GetEntry() or TBranch::GetEntry()
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
    for(ii=0;ii<12;ii++){
        eff.xeff[ii]=0;
        eff.yeff[ii]=0;
    }
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
    layermult.xlmult=0;
    layermult.ylmult=0;
    GetEntry(entry);

    // GetGraphs
    GetXGraph();
    GetYGraph();

    //    Multiplicity and Tdc
    for(ii=0;ii<12;ii++){
        multiplicity.multx[ii]=xstriphits[ii]->CountBits();
        multiplicity.multy[ii]=ystriphits[ii]->CountBits();
        if(multiplicity.multx[ii]>0) layermult.xlmult++;
        if(multiplicity.multy[ii]>0) layermult.ylmult++;
        tdc.xtdc[ii]=tdcdata[ii];
        tdc.ytdc[ii]=tdcdata[ii+16];
        tdc.xtdcmult[ii]=tdcmult[ii];
        tdc.ytdcmult[ii]=tdcmult[ii+16];
    }



    if(option==TString("MAG")){
        FitData2(fxgraph,xfitparam2,multiplicity.multx);
        FitData2(fygraph,yfitparam2,multiplicity.multy);
        GetEfficiency2(xfitparam2,fxgraph,multiplicity.multx,eff.xeff);
        GetEfficiency2(yfitparam2,fygraph,multiplicity.multy,eff.yeff);
    }
    else{
        FitData1(fxgraph,xfitparam1,multiplicity.multx);
        FitData1(fygraph,yfitparam1,multiplicity.multy);
        SetAngle1(xfitparam1,yfitparam1);
        GetEfficiency1(xfitparam1,fxgraph,multiplicity.multx,eff.xeff);
        GetEfficiency1(yfitparam1,fygraph,multiplicity.multy,eff.yeff);
    }
    fittree->Fill();
    // Filling your histograms
    //    if(angle.zen>0)
    //    histo->Fill(angle.zen);
    //    graph->SetPoint(ENum, ENum, angle.zen);
    return kTRUE;

    return kTRUE;
}

TGraphErrors* ino_analysis::GetXGraph(){
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
                }
            }
        }

    }
    fxgraph->Set(xxindx);
    fxgraph->Sort(); // sort along layer
    return fxgraph;
}

TGraphErrors* ino_analysis::GetYGraph(){

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
                }
            }
        }

    }
    fygraph->Set(yyindx);
    fygraph->Sort(); // sort along layer
    return fygraph;
}

void ino_analysis::SetFitError1(FITPARAM1 &f,int errCode){
    f.p0=errCode;
    f.p1=errCode;
    f.errp0=errCode;
    f.errp1=errCode;
    f.chisq=errCode;
    f.ndf=errCode;
}
void ino_analysis::SetFitError2(FITPARAM2 &f,int errCode){
    f.p0=errCode;
    f.p1=errCode;
    f.p2=errCode;
    f.errp0=errCode;
    f.errp1=errCode;
    f.errp2=errCode;
    f.chisq=errCode;
    f.ndf=errCode;
}
void ino_analysis::SetFitParameters1(FITPARAM1 &f,TF1* func){
    f.p0=func->GetParameter(0);
    f.p1=func->GetParameter(1);
    f.errp0=func->GetParError(0);
    f.errp1=func->GetParError(1);
    f.chisq=func->GetChisquare();
    f.ndf=func->GetNDF();
}
void ino_analysis::SetFitParameters2(FITPARAM2 &f,TF1* func){
    f.p0=func->GetParameter(0);
    f.p1=func->GetParameter(1);
    f.p2=func->GetParameter(2);

    f.errp0=func->GetParError(0);
    f.errp1=func->GetParError(1);
    f.errp2=func->GetParError(2);

    f.chisq=func->GetChisquare();
    f.ndf=func->GetNDF();
}

void ino_analysis::SetAngle1(FITPARAM1 &fx, FITPARAM1 &fy){
    if(fx.chisq<=-9990 || fy.chisq<=-9990){
        angle.zen=-9990;
        angle.azi=-9990;
        return;
    }
    xa = fx.p0;
    xb = fx.p1*11+fx.p0;
    ya = fy.p0;
    yb = fy.p1*11+fy.p0;
    za = 0;
    zb = 11;
    hyp=pow((pow((xb-xa)*(STRIPWIDTH+STRIPGAP),2)+pow((yb-ya)*(STRIPWIDTH+STRIPGAP),2)+pow((zb-za)*(LAYERGAP),2)),0.5);
    adj=11*LAYERGAP;
    angle.zen=acos(adj/hyp);
    // azimuth angle calculation in (0,90) (DS)
    //adj=yb*(STRIPWIDTH+STRIPGAP);
    //opp=xb*(STRIPWIDTH+STRIPGAP);
    //angle.azi=atan2(opp,adj);
    // azimuthal angle calculation in (0,360) (SDG)
    adj=yb-ya;
    opp=xa-xb;  // xb-xa
    if(fabs(opp)<0.01)
        opp=0;
    if(fabs(adj)<0.01)
        adj=0;
    if(opp==0 && adj==0)
        angle.azi = 0;
    else
        angle.azi=atan2(adj,opp);
    if(angle.azi>=-3.141593&&angle.azi<0)
        angle.azi=TMath::TwoPi()+angle.azi;
}
/**********************************************/
//Fits a straight line to the data and stores in a FITPARAM structure
//Layers with multiplicity>MAXMULT (2) or ==0 are rejected
//For Layers with more than 1 hit, the average is taken as the hit point
//First fit is rejected (err=-9994) if number of layer hits is < 4
//A second fit is done if first fit is successful but with a Red.chisq>10
//Before the second fitting, the residual wrt the first fit is calculated
//and points away from the fit line by greater than 2 are removed.
//Second fit is rejected (err=-9996) if number of layer hits is < 4
// If the fitting routine itself fails -9990 is returned.
/**********************************************/
void ino_analysis::FitData1(TGraphErrors *g,FITPARAM1& f, int mult[12])
{
    // IMPORTANT !!!
    // graphs should be sorted (along X) prior to calling and multiplicity stored correctly
    // X == Layers and Y == Strips
    // Data is now in an sorted array

    // Init
    fpts=0;pt=0;
    layer=0;
    for(mm=0;mm<MAXMULT;mm++)
    {
        strip[mm]=0;
    }

    fitgraph->Set(0);
    // remove points if multiplicity >MAXMULT
    for(mm=0;mm<12;mm++)
    { // loop over multiplicity
        if(mult[mm]>0 && mult[mm]<=MAXMULT)
        { // if multiplicity < maxmult
            pt= TMath::BinarySearch(g->GetN(),g->GetX(),(double)mm);// find the position where the layer is (Layer Number given by mm)
            ///!!!WARNING IN BINARY SEARCH!!!
            ///If match is found, function returns position of element.
            ///If no match found, function gives nearest element smaller than value.=> means another layer could be returned
            for(kk=0;kk<mult[mm];kk++)
            { // store the strips corresponding to a layer in the array strip[kk]
                g->GetPoint(pt+kk,layer,strip[kk]); // get the ypt (strip)
                if(layer!=mm)
                {
                    // in principle layer returned should be equal to mm if the search returned successfully,
                    //else not correct match: binary search failed
                    break;
                }
            }
            if(layer==mm)
            {///Only if correct match, needed to rectify the bug in binary search (see warning note above)

                // checking for consecutive strips (maximum gap between any two strips is 2 strips)

                for(kk=0;kk<mult[mm];kk++)
                {
                    if(fabs(strip[0]-strip[kk])>MAXSHIFT)
                    {
                        goodlayer[mm]=false;
                        break;
                    }
                    else
                        goodlayer[mm]=true;
                }

                if(goodlayer[mm])
                {
                    fitgraph->SetPoint(fpts,layer,TMath::Mean(mult[mm],strip)); // get the mean and store the point
                    fitgraph->SetPointError(fpts,SIGMALAYER,SIGMASTRIP);
                    fpts++;
                }
            }
        }
    }
    fitgraph->Set(fpts);

    if(fpts<MINLAY)     // Min Number of Layers reqd.
    {
        SetFitError1(f,-9994);
        return;
    }

    if(fitgraph->Fit(fitfunction,"NQO")==0) // First fit is successful==0
    {
        SetFitParameters1(f,fitfunction);
        if(f.chisq/f.ndf >10)   // go for next fit
        {
            for(mm=0;mm<fitgraph->GetN();mm++)
            {
                fitgraph->GetPoint(mm,layer,stripNum); // get the points
                if(fabs(stripNum-(f.p0+f.p1*layer))>(MAXSHIFT-1))   // remove point if difference is greater than 1 strip
                {
                    fitgraph->RemovePoint(mm);
                    mm--; fpts--;
                }
            }
            if(fpts<MINLAY)
            {
                SetFitError1(f,-9996);
                return;
            }
            fitgraph->Set(fpts);

            if(fitgraph->Fit(fitfunction,"NQO")==0) // Second Fitting is successful==0
            {
                SetFitParameters1(f,fitfunction);

                if(!CheckTrackCoverage(f))
                    {
                        SetFitError1(f,-9995);
                    }
            }
            else
            {
                SetFitError1(f,-9990);
            }
        }
        else
        {
            if(!CheckTrackCoverage(f))
            {
                SetFitError1(f,-9995);
            }
        }

    }
    else
    {
        SetFitError1(f,-9990);
    }

}

void ino_analysis::FitData2(TGraphErrors *g,FITPARAM2& f, int mult[12])
{
    // IMPORTANT !!!
    // graphs should be sorted (along X) prior to calling and multiplicity stored correctly
    // X == Layers and Y == Strips
    // Data is now in an sorted array

    // Init
    fpts=0;pt=0;
    layer=0;
    for(mm=0;mm<MAXMULT;mm++){
        strip[mm]=0;
    }
    fitgraph->Set(0);
    // remove points if multiplicity >MAXMULT
    for(mm=0;mm<12;mm++){ // loop over multiplicity
        if(mult[mm]>0 && mult[mm]<=MAXMULT){ // if multiplicity < maxmult
            pt= TMath::BinarySearch(g->GetN(),g->GetX(),(double)mm);// find the position where the layer is (Layer Number given by mm)
            ///!!!WARNING IN BINARY SEARCH!!!
            ///If match is found, function returns position of element.
            ///If no match found, function gives nearest element smaller than value.=> means another layer could be returned
            for(kk=0;kk<mult[mm];kk++){ // store the strips corresponding to a layer in the array strip[kk]
                g->GetPoint(pt+kk,layer,strip[kk]); // get the ypt (strip)
                if(layer!=mm){
                    // in principle layer returned should be equal to mm if the search returned successfully,
                    //else not correct match: binary search failed
                    break;
                }
            }
            if(layer==mm){///Only if correct match, needed to rectify the bug in binary search (see warning note above)
                fitgraph->SetPoint(fpts,layer,TMath::Mean(mult[mm],strip)); // get the mean and store the point
                fitgraph->SetPointError(fpts,SIGMALAYER,SIGMASTRIP);
                fpts++;
            }
        }
    }
    fitgraph->Set(fpts);
    if(fpts<MINLAY){ // Min Number of Layers reqd.
        SetFitError2(f,-9994);
        return;
    }
    if(fitgraph->Fit(fitfunction,"NQO")==0){//==0
        SetFitParameters2(f,fitfunction);
        if(f.chisq/f.ndf >10){ // go for next fit
            for(mm=0;mm<fitgraph->GetN();mm++){
                fitgraph->GetPoint(mm,layer,stripNum); // get the points
                if(fabs(stripNum-(f.p0+f.p1*layer))>2){// remove point if difference is greater
                    fitgraph->RemovePoint(mm);
                    mm--; fpts--;
                }
            }
            if(fpts<4){
                SetFitError2(f,-9996);
                return;
            }
            fitgraph->Set(fpts);
            // Second Fitting
            if(fitgraph->Fit(fitfunction,"NQO")==0){//==0
                SetFitParameters2(f,fitfunction);
            }
            else{
                SetFitError2(f,-9990);
            }
        }

    }
    else{
        SetFitError2(f,-9990);
    }

}
/**********************************************/
//Efficiency bit is set to 1 for each layer in each event provided:
//(a)multiplicity is !=0 and < 4
//(b)Absolute difference between hit and fit is less than 2
// Tracking efficiency is calculated as the ratio of
// Number of '1's in a layer where (chisq>=0)
// by the total number accepted events (chisq>=0)
// chisq>=0 removes all rejected events with error codes in fitting
/**********************************************/
void  ino_analysis::GetEfficiency1(FITPARAM1 &f, TGraphErrors *g,int mult[12], int eff[12])
{
    for(ii=0;ii<12;ii++)
    {
        if(mult[ii]==0 || mult[ii]>4)
        {
            eff[ii]=0;
        }
        else
        {
            pt= TMath::BinarySearch(g->GetN(),g->GetX(),(double)ii);// find the point where the layer is;
            for(jj=0;jj<mult[ii];jj++)
            {
                g->GetPoint(pt+jj,layer,stripNum); // get the ypt (strip)
                if(fabs(stripNum-(f.p0+f.p1*layer))<=(MAXSHIFT-1))
                {
                    eff[ii]=1;
                    break;
                }

            }
        }

    }

}
void  ino_analysis::GetEfficiency2(FITPARAM2 &f, TGraphErrors *g,int mult[12], int eff[12])
{
    for(ii=0;ii<12;ii++){
        if(mult[ii]==0 || mult[ii]>4){
            eff[ii]=0;
        }
        else {
            pt= TMath::BinarySearch(g->GetN(),g->GetX(),(double)ii);// find the point where the layer is;
            for(jj=0;jj<mult[ii];jj++){
                g->GetPoint(pt+jj,layer,stripNum); // get the ypt (strip)
                if(fabs(stripNum-(f.p0+f.p1*layer+f.p2*layer*layer))<2){
                    eff[ii]=1;
                    break;
                }

            }
        }

    }

}

bool ino_analysis::CheckTrackCoverage(FITPARAM1 &f)
{
    // checking if fitted track is confined within the fiducial volume of layer 0 and layer 11
    strip0 = f.p0+f.p1*0;
    strip11 = f.p0+f.p1*11;

    if((strip0<0 || strip0>31) || (strip11<0 || strip11>31))
        return false;
    else
        return true;
}

void ino_analysis::SlaveTerminate()
{
    // The SlaveTerminate() function is called after all entries or objects
    // have been processed. When running with PROOF SlaveTerminate() is called
    // on each slave server.


    if (fFile) {
        Bool_t cleanup = kFALSE;
        TDirectory *savedir = gDirectory;
        if (fittree->GetEntries() > 0) {
            fFile->cd();
            fittree->Write();
            //fProofFile->Print();
            //fOutput->Add(fProofFile);
        }
        else {
            cleanup = kTRUE;
        }
        fittree->SetDirectory(0);
        gDirectory = savedir;
        fFile->Close();
        // TFile::Open(outfileName);// required as a workaround for a ROOT bug

        // Cleanup, if needed
        // if (cleanup) {
        //     TUrl uf(*(fFile->GetEndpointUrl()));
        //     SafeDelete(fFile);
        //     gSystem->Unlink(uf.GetFile());
        //     SafeDelete(fProofFile);
        // }

    }
}

void ino_analysis::Terminate()
{
    // The Terminate() function is the last function to be called during
    // a query. It always runs on the client, it can be used to present
    // the results graphically or save the results to file.

}
ino_analysis::~ino_analysis()
{
    TFile f(outFileName); // this is required as a workaround for a bug in PROOF, likely to be removed in version 5.32
    f.Close();
}


void ino_analysis::Init(TTree *tree)
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

Bool_t ino_analysis::Notify()
{
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.

    return kTRUE;
}