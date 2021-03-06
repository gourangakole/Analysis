
#include "hplusAnalyzer.h"

using namespace std;
void hplusAnalyzer::CutFlowAnalysis(TString url, string myKey, string evtType){
  
  TString outFile("13TeV/outputDir/");
  TString Filename_ = outFile+evtType+"_Anal.root";
  TFile *outFile_ = TFile::Open( Filename_, "RECREATE" );
  outFile_->SetCompressionLevel( 9 );
  
  TString debug_Filename_ = Filename_+"_debug.txt";
  string debug_file(debug_Filename_);
  outfile_.open(debug_file.c_str());

  CutFlowProcessor(url, myKey, "base", outFile_);
  outfile_.close();
  outFile_->Write(); 
  outFile_->Close(); 
}

//---------------------------------------------------//
//Process the cuts, event by event
//---------------------------------------------------//  
void hplusAnalyzer::CutFlowProcessor(TString url,  string myKey, TString cutflowType, TFile *outFile_){
  
  int input_count = 0;
  bool isPFlow = (myKey.find("PFlow") != string::npos) ? true : false;
  string eAlgo("Electrons"), mAlgo("Muons"), jAlgo("Jets"), metAlgo("METs");
  
  //Uncertainty variations, JES, JER, MET unclustered, bTag
  ///int jes = 0, jer = 0, metuc = 0, bscale = 0;
  int jes = 1, jer = 1; 
  if(cutflowType.Contains("JESPlus"))jes = 1;
  else if (cutflowType.Contains("JESMinus"))jes = -1;
  else if (cutflowType.Contains("JERPlus"))jer = 1;
  else if (cutflowType.Contains("JERMinus"))jer = -1;
  //else if (cutflowType.Contains("METUCPlus"))metuc = 1;
  //else if (cutflowType.Contains("METUCMinus"))metuc = -1;
  //else if (cutflowType.Contains("bTagPlus"))bscale = 1;
  //else if (cutflowType.Contains("bTagMinus"))bscale = -1; 

      evR = new Reader();
      TFile *f = TFile::Open(url);
      if(f==0) return ;
      if(f->IsZombie()) { f->Close(); return; }
     
      //---------------------------------------------------//
      //get initial number of events, from ntuples
      //store initial informations, in a txt file
      //---------------------------------------------------//
      double lumiTotal = 34312;
      int nEntries = evR->AssignEventTreeFrom(f);
      if( nEntries == 0) {return; }
      TH1F* inputcf = (TH1F*)(f->Get("allEventsFilter/totalEvents"))->Clone("inputcf");
      TH1F* intimepu = (TH1F*)(f->Get("myMiniTreeProducer/MCINFO/intimepu"))->Clone("intimepu");
      TH1F* outoftimepu = (TH1F*)(f->Get("myMiniTreeProducer/MCINFO/outoftimepu"))->Clone("outoftimepu");
      TH1F* totalpu = (TH1F*)(f->Get("myMiniTreeProducer/MCINFO/totalpu"))->Clone("totalpu");
      TH1F* trueintimepu = (TH1F*)(f->Get("myMiniTreeProducer/MCINFO/trueintimepu"))->Clone("trueintimepu");
      TH1F* trueoutoftimepu = (TH1F*)(f->Get("myMiniTreeProducer/MCINFO/trueoutoftimepu"))->Clone("trueoutoftimepu");
      TH1F* truetotalpu = (TH1F*)(f->Get("myMiniTreeProducer/MCINFO/truetotalpu"))->Clone("truetotalpu");
      double initialEvents = inputcf->GetBinContent(1);
      cout<<"input file: "<<url<<endl;
      outfile_<<"input file: "<<url<<endl;
      outfile_<<"totalEvents: "<<initialEvents<<endl;
      CreateAnalHistos(cutflowType, outFile_);
      for(int nbin=1; nbin<6001; nbin++){
        fillHistoPU("intimepu", cutflowType, nbin, intimepu->GetBinContent(nbin));  
        fillHistoPU("outoftimepu", cutflowType, nbin, outoftimepu->GetBinContent(nbin));  
        fillHistoPU("totalpu", cutflowType, nbin, totalpu->GetBinContent(nbin));  
        fillHistoPU("trueintimepu", cutflowType, nbin, trueintimepu->GetBinContent(nbin));  
        fillHistoPU("trueoutoftimepu", cutflowType, nbin, trueoutoftimepu->GetBinContent(nbin));  
        fillHistoPU("truetotalpu", cutflowType, nbin, truetotalpu->GetBinContent(nbin));  
      }
      fillHisto("totalEvents", cutflowType, initialEvents, 1);
      MyEvent *ev;
      int nTriggEvent = 0, nSelEvents = 0, matchjetcount= 0, threepairjet = 0;
      double nVerticesFailCount = 0.0;
      ///double corr_jet_pair_cs = 0.0, corr_jet_pair_ud = 0.0;
      double matchedJet_q = 0.0, matchedJet_b = 0.0, matched_quark_eta_pt = 0.0, not_matchedJet_q = 0.0;
      double TotalTopPtWeights = 0, TotalLplusJEvents = 0; 
      double TotalSVEff = 0, TotalEvtsWithSV = 0.;
      
      //---------------------------------------------------//
      //loop over each event, of the ntuple
      //---------------------------------------------------//
      double kfCount = 0;
      for(int i=0; i<nEntries; ++i){
        Long64_t ientry = evR->LoadTree(i);
        if (ientry < 0) break;
        ev = evR->GetNewEvent(i);
        if(ev==0) continue;
        if(i%1000==0) cout<<"\033[01;32mEvent number = \033[00m"<< i << endl;
     
        //---------------------------------------------------//
        //apply lumi, k factor and pileup weight
        //---------------------------------------------------//
        double evtWeight = 1.0;
        if(!ev->isData){
          string sampleName = ev->sampleInfo.sampleName;
          //k-factor weight (along with lumi weight) 
          if(sampleName.find("WJetsToLNu") != string::npos || sampleName.find("W1JetsToLNu") != string::npos || sampleName.find("W2JetsToLNu") != string::npos || sampleName.find("W3JetsToLNu") != string::npos || sampleName.find("W4JetsToLNu") != string::npos){
            int hepNUP = ev->sampleInfo.hepNUP;
            double weightK = reweightHEPNUPWJets(hepNUP) * (lumiTotal/1000.0);
            evtWeight *= weightK;  
            fillHisto("hepNUP", cutflowType, hepNUP, evtWeight);
          }
          else if(sampleName.find("DYJetsToLL") != string::npos || sampleName.find("DY1JetsToLL") != string::npos || sampleName.find("DY2JetsToLL") != string::npos || sampleName.find("DY3JetsToLL") != string::npos || sampleName.find("DY4JetsToLL") != string::npos){
            int hepNUP = ev->sampleInfo.hepNUP;
            std::vector<int> hepIDUP = ev->sampleInfo.hepIDUP;
            std::vector<int> hepISTUP = ev->sampleInfo.hepISTUP;
            int countZ = 0;
            for(size_t p=0; p<hepIDUP.size(); p++){
              if(hepIDUP[p]==23 && hepISTUP[p]==2)
                countZ = countZ + 1;
            }
            if(countZ==0) hepNUP = hepNUP+1;
            double weightK = reweightHEPNUPDYJets(hepNUP) * (lumiTotal/1000.0);
            evtWeight *= weightK;  
            fillHisto("hepNUP", cutflowType, hepNUP, evtWeight);
          }
          //lumi weight
          else {
          double sampleWeight(1.0);
          sampleWeight = lumiTotal* xss[sampleName]/evtDBS[sampleName];
          evtWeight *= sampleWeight; 
          if(i < 1){
            outfile_<<"Sample weight for "<<sampleName<<" is= "<< sampleWeight<<endl;
            }
          }
          //pileup weight 
          vector<double>pu = ev->sampleInfo.truepileup;
          if(pu.size() > 0) {
          float npu = pu[0];
          double weightPU = LumiWeights_.weight(npu);
          evtWeight *= weightPU;  
          }
        } 
/*
        //---------------------------------------------------//
        // apply top re-weighting weight
        //---------------------------------------------------//
        double topPtWeights_offline = 1.0;
        if(!ev->isData){
          string sampleName = ev->sampleInfo.sampleName;
          //(we have to update it ==========)
          if(sampleName.find("TTJets") != string::npos || sampleName.find("WH") != string::npos || sampleName.find("HplusM120") != string::npos){
            vector<double>topptweights = ev->sampleInfo.topPtWeights;
            if(topptweights.size() > 0){
              topPtWeights_offline = topptweights[0];
              if(cutflowType.Contains("TopPtPlus"))
                topPtWeights_offline = topPtWeights_offline*topPtWeights_offline;
              else if(cutflowType.Contains("TopPtMinus"))
                topPtWeights_offline = 1.0;
            }
          }
        }
        TotalTopPtWeights += topPtWeights_offline;
        //cout<<"topPtWeights_offline = "<<topPtWeights_offline<<endl;
        TotalLplusJEvents++;
        evtWeight *= topPtWeights_offline; //Multiply to the total weights
  */      
        //---------------------------------------------------//
        //apply muon triggers
        //---------------------------------------------------//
        double nCutPass = 0.0;
        fillHisto("cutflow", cutflowType+"/Iso", nCutPass, evtWeight);
        double nCutPass_NonIso = 0.0;
        fillHisto("cutflow", cutflowType+"/NonIso", nCutPass_NonIso, evtWeight);
        
        bool passTrig = false;
        vector<string> trig = ev->hlt;
        for(size_t it = 0; it < trig.size(); it++){
          if(trig[it].find("Mu") != string::npos) passTrig = true;
        }
        if(!passTrig){
                //cout << "not satisfying trigger" << endl;
          continue;
        }
        nTriggEvent++;
        nCutPass++;
        nCutPass_NonIso++;
        fillHisto("cutflow", cutflowType+"/Iso", nCutPass, evtWeight);
        fillHisto("cutflow", cutflowType+"/NonIso", nCutPass_NonIso, evtWeight);
       
        //---------------------------------------------------//
        //get all objets e.g. leptons, jets, vertices etc.
        //---------------------------------------------------//
        vector<MyVertex> Vertices = ev->PrimaryVtxs;
        if(Vertices.size() <= 0){
          nVerticesFailCount+=evtWeight;
          cout<<" no vertexes , exit"<<endl;
          continue;
        }
        vector<MyMuon> pfMuons = evR->getMuons(ev, mAlgo);
        vector<MyMuon> pfMuons_noiso = evR->getMuons(ev, mAlgo);
        vector<MyElectron> pfElectrons = evR->getElectrons(ev, eAlgo);
        vector<MyJet> pfJets = evR->getJets(ev, jAlgo);
        MyMET met = evR->getMET(ev, metAlgo);

        // preselect objects 
        vector<int> m_init; m_init.clear();
        preSelectMuons(&m_init, pfMuons, Vertices[0], isPFlow);
        vector<int> m_init_noiso; m_init_noiso.clear();
        preSelectMuonsNoIso(&m_init_noiso, pfMuons_noiso, Vertices[0], isPFlow);
        vector<int> e_init; e_init.clear();
        preSelectElectrons(&e_init, pfElectrons, Vertices[0], isPFlow);
        vector<int> j_init; j_init.clear();
        preSelectJets(jAlgo, &j_init, pfJets, jes, jer);
        
        // clean objects //
        vector<int> e_final; e_final.clear();
        ElectronCleaning( pfElectrons, pfMuons, &e_init, &e_final, &m_init, DRMIN_ELE);
        vector<int> j_final; j_final.clear();
        vector<int> t_final; t_final.clear();
        JetCleaning(pfJets, pfMuons, pfElectrons,  &j_init, &j_final, &m_init, &e_final, DRMIN_JET);
        
        //Get MC partons
        vector<MyLorentzVector> bquarks; bquarks.clear();
        vector<MyLorentzVector> lquarks; lquarks.clear();
        MyLorentzVector mcTop, mcAntiTop;
        mcTop.SetPxPyPzE(0., 0., 0., 0.); mcAntiTop.SetPxPyPzE(0., 0., 0., 0.);
        if(!ev->isData){
          vector<MyMCParticle>allMCParticles = ev->mcParticles;
          for(size_t imc=0; imc < allMCParticles.size(); ++imc){
            if(abs(allMCParticles[imc].pid) == 5 && allMCParticles[imc].mother.size() > 0 && (abs(allMCParticles[imc].mother[0])==6) )
              bquarks.push_back(allMCParticles[imc].p4Gen);
            else if(abs(allMCParticles[imc].pid) <= 4 && allMCParticles[imc].mother.size() > 0 && (abs(allMCParticles[imc].mother[0])==24 || abs(allMCParticles[imc].mother[0])==37) )
              lquarks.push_back(allMCParticles[imc].p4Gen); 
            else if(allMCParticles[imc].pid == 6 && allMCParticles[imc].status == 3)
              mcTop = allMCParticles[imc].p4Gen;
            else if(allMCParticles[imc].pid == -6 && allMCParticles[imc].status == 3)
              mcAntiTop = allMCParticles[imc].p4Gen;
          }
        }

        //---------------------------------------------------//
        //get KinFit objects
        //---------------------------------------------------//
        vector<MyLorentzVector> kfJets; kfJets.clear();
        vector<MyLorentzVector> kfLepton; kfLepton.clear();
        vector<MyLorentzVector> kfMet; kfMet.clear();
     
        double chi2OfKinFit=999.;
        double statusOfKinFit=-99;
        double probOfKinFit=-99;
        vector<MyKineFitParticle> allKineFitParticles = ev->KineFitParticles;
	for(size_t imk=0; imk < allKineFitParticles.size(); imk++){
          string labelName = "";
	  if(cutflowType.Contains("JESPlus"))labelName="JESUp";
          else if(cutflowType.Contains("JESMinus"))labelName="JESDown";
          else if(cutflowType.Contains("JERPlus"))labelName="JERUp";
          else if(cutflowType.Contains("JERMinus"))labelName="JERDown";
          if(labelName=="JESUp" || labelName=="JESDown" || labelName=="JERUp" ||labelName=="JERDown"){
            if(allKineFitParticles[imk].labelName.find(labelName) != string::npos ){
              if(allKineFitParticles[imk].partName.find("PartonsHadP") != string::npos ||
                 allKineFitParticles[imk].partName.find("PartonsHadQ") != string::npos ||
                 allKineFitParticles[imk].partName.find("PartonsHadB") != string::npos
                 ){
                kfJets.push_back(allKineFitParticles[imk].p4);
              }
              else if(allKineFitParticles[imk].partName.find("Leptons") != string::npos ){
                kfLepton.push_back(allKineFitParticles[imk].p4);
              }
              else if(allKineFitParticles[imk].partName.find("Neutrinos") != string::npos ){
                kfMet.push_back(allKineFitParticles[imk].p4);
              }
              if(imk<1){
                chi2OfKinFit = allKineFitParticles[imk].chi2OfFit;
                statusOfKinFit = allKineFitParticles[imk].statusOfFit;
                probOfKinFit = allKineFitParticles[imk].probOfFit;
              }
            }
          }
          else{
             if(allKineFitParticles[imk].labelName.find("JESUp") == string::npos && allKineFitParticles[imk].labelName.find("JESDown") == string::npos &&
                allKineFitParticles[imk].labelName.find("JERUp") == string::npos && allKineFitParticles[imk].labelName.find("JERDown") == string::npos)
               {
                 if(allKineFitParticles[imk].partName.find("PartonsHadP") != string::npos ||
                    allKineFitParticles[imk].partName.find("PartonsHadQ") != string::npos ||
                    allKineFitParticles[imk].partName.find("PartonsHadB") != string::npos
                    ){
                   kfJets.push_back(allKineFitParticles[imk].p4);
                 }
                 else if(allKineFitParticles[imk].partName.find("Leptons") != string::npos ){
                   kfLepton.push_back(allKineFitParticles[imk].p4);
                 }
                 else if(allKineFitParticles[imk].partName.find("Neutrinos") != string::npos ){
                   kfMet.push_back(allKineFitParticles[imk].p4);
                 }
                 if(imk<1){
                   chi2OfKinFit = allKineFitParticles[imk].chi2OfFit;
                   statusOfKinFit = allKineFitParticles[imk].statusOfFit;
                   probOfKinFit = allKineFitParticles[imk].probOfFit;
                 }
               }
           }
	}
        //---------------------------------------------------//
        //apply selection cuts on leptons
        //---------------------------------------------------//
        if(m_init_noiso.size() > 0){
          int m_i_noiso = m_init_noiso[0];
          double mRelIso_no_iso = pfMuons_noiso[m_i_noiso].pfRelIso;
          fillHisto("pre_RelIso_mu",cutflowType+"/Iso", mRelIso_no_iso, evtWeight);
          fillHisto("pre_RelIso_mu",cutflowType+"/NonIso", mRelIso_no_iso, evtWeight);
        }
        bool noisofound = false;
        bool isofound = false;
        int nLepton = m_init_noiso.size();
        //double pri_vtxs = Vertices.size();
        double pri_vtxs = Vertices[0].totVtx;
        
        if(nLepton != 1)continue;
        if( looseMuonVeto( m_init_noiso[0],pfMuons_noiso, isPFlow) ) continue;
        if( looseElectronVeto(-1,pfElectrons, isPFlow) ) continue;

        int m_i = m_init_noiso[0];
        double tmp_iso = pfMuons_noiso[m_i].pfRelIso;
        if(tmp_iso > 0.12 && tmp_iso < 0.3){
          noisofound = true;
        }
        fillHisto("RelIso_mu",cutflowType, tmp_iso, evtWeight);
        if(tmp_iso < 0.12) isofound = true;

        //Transverse mass b/w lepton and MET
        double metPt = 0; 
	double   leptonPt(0), deltaPhi(0);
	metPt= met.p4.pt();
        leptonPt = TMath::Abs(pfMuons[m_i].p4.pt());
        deltaPhi = ROOT::Math::VectorUtil::DeltaPhi(pfMuons[m_i].p4, met.p4);
        double mt = sqrt (  2*leptonPt*metPt*(1 - cos(deltaPhi) ) ) ;
        fillHisto("wmt_1mu", cutflowType+"/Iso", mt, evtWeight);
	
	////////////////////////////////////////////////////////////////  
        //  Isolation region
        ////////////////////////////////////////////////////////////////  
        if(isofound){
          //apply muon SF to eventWeights 
           
          //double musfWeight = 1.0;
          //if(fabs(pfMuons[m_i].p4.eta()) < 0.9)musfWeight = muSF["sfEta1"];
          //else if(fabs(pfMuons[m_i].p4.eta()) > 0.9 && fabs(pfMuons[m_i].p4.eta()) < 1.2)musfWeight = muSF["sfEta2"];
          //else musfWeight = muSF["sfEta3"];
          //evtWeight *= musfWeight;
          //cout<<"evtWeight musf = "<<evtWeight<<endl;
             
          nCutPass++; 
          fillHisto("cutflow", cutflowType+"/Iso", nCutPass, evtWeight); // one electron
          int count_muon = m_init_noiso.size();
          ///int muCharge = pfMuons[m_i].charge;

          ///double nCutPass = 0.0;// double nCutPass_plus = 0.0; double nCutPass_minus = 0.0;
          double mRelIso = pfMuons[m_i].pfRelIso;
          fillHisto("pt_mu", cutflowType+"/Iso", pfMuons[m_i].p4.pt(), evtWeight);
          fillHisto("final_RelIso_mu",cutflowType+"/Iso", mRelIso, evtWeight);
          fillHisto("final_multi_mu",cutflowType+"/Iso", count_muon, evtWeight);

          // Fill histogram after trigger and one offline isolated muon and applied 2nd lepton veto
          int nJet = j_final.size();
          fillHisto("eta_mu", cutflowType+"/Iso", pfMuons[m_i].p4.eta(), evtWeight);
          fillHisto("phi_mu", cutflowType+"/Iso", pfMuons[m_i].p4.phi(), evtWeight);
          // vertex just after one lepton selection
          //double pri_vtxs = Vertices.size();
          fillHisto("nvtx", cutflowType+"/Iso", pri_vtxs, evtWeight);
          fillHisto("nvtx_6Kbins", cutflowType+"/Iso", pri_vtxs, evtWeight);
          for( std::size_t n=0; n<Vertices.size(); n++){
            fillHisto("rhoAll", cutflowType+"/Iso", Vertices[n].rhoAll, evtWeight);
            fillHisto("chi2", cutflowType+"/Iso", Vertices[n].chi2, evtWeight);
            fillHisto("ndof", cutflowType+"/Iso", Vertices[n].ndof, evtWeight);
            }
          //---------------------------------------------------//
          // Apply Jet Selection
          //---------------------------------------------------//

          fillHisto("multi_jet", cutflowType+"/Iso", nJet, evtWeight);
          if(nJet < 4)continue;  // this condition implies event should contain at least 4 jets
          for(size_t ijet = 0; ijet < j_final.size(); ijet++){
            int ind_jet = j_final[ijet];
            fillHisto("pt_jet", cutflowType+"/Iso", pfJets[ind_jet].p4.pt(), evtWeight);
            fillHisto("eta_jet", cutflowType+"/Iso", pfJets[ind_jet].p4.eta(), evtWeight);
            fillHisto("phi_jet", cutflowType+"/Iso", pfJets[ind_jet].p4.phi(), evtWeight);
          }
          fillHisto("wmt_1mu_4jet", cutflowType+"/Iso", mt, evtWeight);
          fillHisto("final_multi_jet", cutflowType+"/Iso", nJet, evtWeight);
          nCutPass++;
          fillHisto("cutflow", cutflowType+"/Iso", nCutPass, evtWeight);
          
          //---------------------------------------------------//
          //apply MET selection   
          //---------------------------------------------------//
	  double metPtJESJER =0;
	  metPtJESJER= metWithJESJER(pfJets, &j_final, met, 0, 0);
          fillHisto("pt_metJESJER", cutflowType+"/Iso", metPtJESJER, evtWeight);
	  ///if(!metuc)metPt = metWithJESJER(pfJets, &j_final, met, jes, jer);
          ///else metPt = metWithUncl(pfJets, &j_final, pfMuons, &m_init, pfElectrons, &e_final, met, metuc);
          fillHisto("pt_met", cutflowType+"/Iso", metPt, evtWeight);
          fillHisto("phi_met", cutflowType+"/Iso", met.p4.phi(), evtWeight);
          if(metPt < 20) continue;  // Missing transverse energy cut 30 GeV(CMS) for ATLAS 20 GeV 
          fillHisto("final_pt_met", cutflowType+"/Iso", metPt, evtWeight);
          fillHisto("final_phi_met", cutflowType+"/Iso", met.p4.phi(), evtWeight);
          nCutPass++;
          fillHisto("cutflow", cutflowType+"/Iso", nCutPass, evtWeight);
          fillHisto("wmt", cutflowType+"/Iso", mt, evtWeight);
          
          //if((mt + metPt ) < 30 ) continue;  // extra condition to reduce multijet bkgs
          //nCutPass++;
          //fillHisto("cutflow", cutflowType+"/Iso", nCutPass, evtWeight);
          
          
          //---------------------------------------------------//
          //apply B-tagging, C-tagging
          //---------------------------------------------------//
          vector<int> j_final_nob; j_final_nob.clear();
          double pfCISV = 0.0; //pfCombinedInclusiveSecondaryVertexV2BJetTags
          double pfCMVA = 0.0; //pfCombinedMVAV2BJetTags
          double pfCCvsL = 0.0;//pfCombinedCvsLJetTags
          double pfCCvsB = 0.0; //pfCombinedCvsBJetTags
          double pfCCvsL_0 = 0.0;
          double pfCCvsL_1 = 0.0;
          double pfCCvsB_0 = 0.0;
          double pfCCvsB_1 = 0.0;
          int bscale = 0;

	  //LOOSE BTAG
	  int count_CSVL = 0; int count_CSVL_SF = 0; 
          for(size_t ijet = 0; ijet < j_final.size(); ijet++){
            int ind_jet = j_final[ijet];
            pfCISV = pfJets[ind_jet].bDiscriminator["pfCombinedInclusiveSecondaryVertexV2BJetTags"];
            pfCMVA = pfJets[ind_jet].bDiscriminator["pfCombinedMVAV2BJetTags"];
            pfCCvsL= pfJets[ind_jet].bDiscriminator["pfCombinedCvsLJetTags"];
            pfCCvsB = pfJets[ind_jet].bDiscriminator["pfCombinedCvsBJetTags"];
            fillHisto("pfCISV", cutflowType+"/Iso/BTag", pfCISV , evtWeight); 
            fillHisto("pfCMVA", cutflowType+"/Iso/BTag", pfCMVA , evtWeight); 
            fillHisto("pfCCvsL", cutflowType+"/Iso/BTag", pfCCvsL, evtWeight); 
            fillHisto("pfCCvsB", cutflowType+"/Iso/BTag", pfCCvsB , evtWeight); 
            
	    bool isBtag = getBtagWithSF(pfJets[ind_jet], ev->isData, bscale, true); 
	    if(isBtag)count_CSVL_SF++; 
	    
	    if(pfCISV > 0.5426){
              count_CSVL++;
            }
            else j_final_nob.push_back(ind_jet);  
          }
	  ///cout<<"count_CSVL = "<<count_CSVL<<endl;
	  ///cout<<"count_CSVL_SF = "<<count_CSVL_SF<<endl;

          fillHisto("CSVL_count", cutflowType+"/Iso/BTag", count_CSVL, evtWeight);
          fillHisto("CSVM_count", cutflowType+"/Iso/BTagM", count_CSVL_SF, evtWeight);
          //Get the invariant mass of csbar pair
          if(j_final_nob.size() >= 2){
            int first_index = j_final_nob[0];
            int sec_index = j_final_nob[1];
            pfCCvsL_0 = pfJets[first_index].bDiscriminator["pfCombinedCvsLJetTags"];
            pfCCvsL_1 = pfJets[sec_index].bDiscriminator["pfCombinedCvsLJetTags"];
            pfCCvsB_0 = pfJets[first_index].bDiscriminator["pfCombinedCvsBJetTags"];
            pfCCvsB_1 = pfJets[sec_index].bDiscriminator["pfCombinedCvsBJetTags"];
            fillHisto("pfCCvsL_0", cutflowType+"/Iso/BTag", pfCCvsL_0, evtWeight); 
            fillHisto("pfCCvsL_1", cutflowType+"/Iso/BTag", pfCCvsL_1, evtWeight); 
            fillHisto("pfCCvsB_0", cutflowType+"/Iso/BTag", pfCCvsB_0 , evtWeight); 
            fillHisto("pfCCvsB_1", cutflowType+"/Iso/BTag", pfCCvsB_1 , evtWeight); 
            MyLorentzVector diJet = pfJets[first_index].p4 + pfJets[sec_index].p4;
            fillHisto("mjj", cutflowType+"/Iso/BTag", diJet.M(), evtWeight);
          }
          if(count_CSVL_SF <= 1) continue; // Demanding for 2L b-tagged jets
	  fillHisto("wmt_1mu_4jet_btag", cutflowType+"/Iso", mt, evtWeight);
          nCutPass++;
          fillHisto("cutflow", cutflowType+"/Iso", nCutPass, evtWeight); 
          nSelEvents++;  // this is the counter for counting final number of events  

          //---------------------------------------------------//
          // add set of plots after BTag:
          //---------------------------------------------------//
          fillHisto("pt_mu", cutflowType+"/Iso"+"/BTag", pfMuons[m_i].p4.pt(), evtWeight);
          fillHisto("eta_mu", cutflowType+"/Iso"+"/BTag", pfMuons[m_i].p4.eta(), evtWeight);
          fillHisto("phi_mu", cutflowType+"/Iso"+"/BTag", pfMuons[m_i].p4.phi(), evtWeight);
          fillHisto("final_RelIso_mu",cutflowType+"/Iso"+"/BTag", mRelIso, evtWeight);
          for(size_t ijet = 0; ijet < j_final.size(); ijet++){
            int ind_jet = j_final[ijet];
            double jetPtJESJER = jetPtWithJESJER(pfJets[ind_jet], 0, 0);
	    fillHisto("pt_jetJESJER", cutflowType+"/Iso"+"/BTag", jetPtJESJER, evtWeight);
	    double jetPt = pfJets[ind_jet].p4.pt(); 
            fillHisto("pt_jet", cutflowType+"/Iso"+"/BTag", jetPt, evtWeight);
            fillHisto("eta_jet", cutflowType+"/Iso"+"/BTag", pfJets[ind_jet].p4.eta(), evtWeight);
            fillHisto("phi_jet", cutflowType+"/Iso"+"/BTag", pfJets[ind_jet].p4.phi(), evtWeight);
          }
          fillHisto("final_multi_jet", cutflowType+"/Iso"+"/BTag", nJet, evtWeight);
          fillHisto("final_pt_met", cutflowType+"/Iso"+"/BTag", metPt, evtWeight);
          fillHisto("nvtx", cutflowType+"/Iso"+"/BTag", pri_vtxs, evtWeight);
          fillHisto("nvtx_6Kbins", cutflowType+"/Iso"+"/BTag", pri_vtxs, evtWeight);
          for(std::size_t n=0; n<Vertices.size(); n++){
            fillHisto("rhoAll", cutflowType+"/Iso"+"/BTag", Vertices[n].rhoAll, evtWeight);
            fillHisto("chi2", cutflowType+"/Iso"+"/BTag", Vertices[n].chi2, evtWeight);
            fillHisto("ndof", cutflowType+"/Iso"+"/BTag", Vertices[n].ndof, evtWeight);
            }
          fillHisto("wmt", cutflowType+"/Iso"+"/BTag", mt, evtWeight);
          
          
          input_count++;
          if(input_count%10==0)
          cout << "input count iso: "<< input_count << endl;
          //if(i > 10000) break;
          //---------------------------------------------------//
          // add set of plots after KinFit:
          //---------------------------------------------------//
          bool foundkfMuon = false;
          if(kfLepton.size()>0){
            if(DeltaR(pfMuons_noiso[m_i].p4 , kfLepton[0]) < 0.2)foundkfMuon = true;
          }
          //outfile_ << "chi2OfKinFit just after lepton selection:    " << chi2OfKinFit << endl;
          int count_cjet = 0, count_sjet = 0, count_bjet = 0, count_ujet = 0, count_djet = 0, other_jet = 0;
          if(foundkfMuon){
            //---------------------------------------------------//
            //select maximum b-tag discriminator jet in KF
            //---------------------------------------------------//
            int maxBtagJet = -1;
            double maxBDiscr = -999.;
            for(size_t ik = 0; ik < kfJets.size(); ik++){
              for(size_t ij = 0; ij < j_final.size(); ij++){
                int ind_ij = j_final[ij];
                if(DeltaR(kfJets[ik], pfJets[ind_ij].p4) < 0.2){
                  double discr = pfJets[ind_ij].bDiscriminator["pfCombinedInclusiveSecondaryVertexV2BJetTags"];
                  if(fabs(pfJets[ind_ij].partonFlavour) == 4) count_cjet++;
                  else if(fabs(pfJets[ind_ij].partonFlavour) == 3) count_sjet++;
                  else if(fabs(pfJets[ind_ij].partonFlavour) == 5) count_bjet++;
                  else if(fabs(pfJets[ind_ij].partonFlavour) == 2) count_ujet++;
                  else if(fabs(pfJets[ind_ij].partonFlavour) == 1) count_djet++;
                  else other_jet++;
                  if(discr > maxBDiscr){
                    maxBDiscr = discr;
                    maxBtagJet = ik;
                  }
                }
              }
            }
            vector<MyLorentzVector> kfLightJets; kfLightJets.clear();
            if(kfJets.size() >=3 && maxBtagJet >= 0){
              for(unsigned ik = 0; ik < kfJets.size(); ik++){
                if(ik != maxBtagJet)kfLightJets.push_back(kfJets[ik]);
              }
            }
            
            //---------------------------------------------------//
            //loop over light KF jets
            //---------------------------------------------------//
            if(kfLightJets.size() >= 2){
              kfCount++;
	      MyLorentzVector diJet = kfLightJets[0]+kfLightJets[1];
	      
	      // category depending on C-tagging
	      //pfCCvsL= pfJets[ind_jet].bDiscriminator["pfCombinedCvsLJetTags"];
	      //pfCCvsB = pfJets[ind_jet].bDiscriminator["pfCombinedCvsBJetTags"];

              fillHisto("mjj_kfit", cutflowType+"/Iso/KinFit", diJet.mass(), evtWeight);
              bool match_j1 = false, match_j2 = false;
	      int indexForCTag0 = 0, indexForCTag1 = 0;
	      
              for(size_t ij = 0; ij < j_final.size(); ij++){
                int ind_ij = j_final[ij];
                if(DeltaR(kfLightJets[0], pfJets[ind_ij].p4) < 0.2){
                  match_j1=true;
		  indexForCTag0 = ind_ij;
		  //cout << "ind_ij for 0th Jet: " << ind_ij << endl;
                }
                if(DeltaR(kfLightJets[1], pfJets[ind_ij].p4) < 0.2){
                  match_j2=true;
		  indexForCTag1 = ind_ij;
		  //cout << "ind_ij for 1st Jet: " << ind_ij << endl;
                }
              }
	      
	      // category depending on C-tagging
              double pfCCvsL_ForCTag0 = pfJets[indexForCTag0].bDiscriminator["pfCombinedCvsLJetTags"];
	      double pfCCvsL_ForCTag1 = pfJets[indexForCTag1].bDiscriminator["pfCombinedCvsLJetTags"];

              double pfCCvsB_ForCTag0 = pfJets[indexForCTag0].bDiscriminator["pfCombinedCvsBJetTags"]; 
	      double pfCCvsB_ForCTag1 = pfJets[indexForCTag1].bDiscriminator["pfCombinedCvsBJetTags"];
	      
	      cout << "=========================" << endl;
	      cout << "pfCCvsL_ForCTag0: " << pfCCvsL_ForCTag0 << endl;
	      cout << "pfCCvsL_ForCTag1: " << pfCCvsL_ForCTag1 << endl;
	      
	      cout << "pfCCvsB_ForCTag0: " << pfCCvsB_ForCTag0 << endl;
	      cout << "pfCCvsB_ForCTag1: " << pfCCvsB_ForCTag1 << endl;
	      cout << "=========================" << endl;

	      //MyLorentzVector diJet = kfLightJets[0]+kfLightJets[1];
	      //CharmL-> (CvsL: -0.48 and CvsB: -0.17)
	      if(match_j1 && match_j2){
		MyLorentzVector diJet_tag = kfLightJets[0]+kfLightJets[1];
		if((pfCCvsL_ForCTag0 > -0.48 || pfCCvsL_ForCTag1 > -0.48) && (pfCCvsB_ForCTag0 > -0.17 || pfCCvsB_ForCTag1 > -0.17)){
		  fillHisto("mjj_kfit_CTag", cutflowType+"/Iso/KinFit", diJet_tag.mass(), evtWeight);
		}else {
		  fillHisto("mjj_kfit_noCTag", cutflowType+"/Iso/KinFit", diJet_tag.mass(), evtWeight);
		}
		
	      }




              if(match_j1 && match_j2){
                if(kfLightJets[0].pt() > 25 && kfLightJets[1].pt() > 25 ){
                  if(probOfKinFit > 0.1){
                    fillHisto("pt_mu", cutflowType+"/Iso/KinFit", pfMuons_noiso[m_i].p4.pt(), evtWeight);
                    fillHisto("eta_mu", cutflowType+"/Iso/KinFit", pfMuons_noiso[m_i].p4.eta(), evtWeight);
                    fillHisto("phi_mu", cutflowType+"/Iso/KinFit", pfMuons_noiso[m_i].p4.phi(), evtWeight);
                    for(size_t ijet = 0; ijet < j_final.size(); ijet++){
                      int ind_jet = j_final[ijet];
                      //double jetPt = jetPtWithJESJER(isData, pfJets[ind_jet], jes, jer);
                      double jetPt = pfJets[ind_jet].p4.pt();
                      fillHisto("pt_jet", cutflowType+"/Iso/KinFit", jetPt, evtWeight);
                      fillHisto("eta_jet", cutflowType+"/Iso/KinFit", pfJets[ind_jet].p4.eta(), evtWeight);
                      fillHisto("phi_jet", cutflowType+"/Iso/KinFit", pfJets[ind_jet].p4.phi(), evtWeight);
                    }
                  }
                  if(probOfKinFit > 0.1){
                    for( std::size_t n=0; n<Vertices.size(); n++){
                      nCutPass++;
                      fillHisto("cutflow", cutflowType+"/Iso", nCutPass, evtWeight); 
                      fillHisto("rhoAll", cutflowType+"/Iso/KinFit", Vertices[n].rhoAll, evtWeight);
                      fillHisto("chi2", cutflowType+"/Iso/KinFit", Vertices[n].chi2, evtWeight);
                      fillHisto("ndof", cutflowType+"/Iso/KinFit", Vertices[n].ndof, evtWeight);
                    }
                    fillHisto("final_multi_jet", cutflowType+"/Iso/KinFit", nJet, evtWeight);
                    fillHisto("final_pt_met", cutflowType+"/Iso/KinFit", metPt, evtWeight);
                    fillHisto("nvtx", cutflowType+"/Iso/KinFit", pri_vtxs, evtWeight);
                    //fillHisto("CSVM_count", cutflowType+"/Iso/KinFit", count_CSVM, evtWeight);
                    fillHisto("wmt", cutflowType+"/Iso/KinFit", mt, evtWeight);
                    fillHisto("kfJet1_pt",cutflowType+"/Iso/KinFit", kfLightJets[0].pt(), evtWeight);
                    fillHisto("kfJet2_pt",cutflowType+"/Iso/KinFit", kfLightJets[1].pt(), evtWeight);
                    fillHisto("kfJet1_eta",cutflowType+"/Iso/KinFit", kfLightJets[0].eta(), evtWeight);
                    fillHisto("kfJet2_eta",cutflowType+"/Iso/KinFit", kfLightJets[1].eta(), evtWeight);
                    fillHisto("kfJet1_phi",cutflowType+"/Iso/KinFit", kfLightJets[0].phi(), evtWeight);
                    fillHisto("kfJet2_phi",cutflowType+"/Iso/KinFit", kfLightJets[1].phi(), evtWeight);
                  }
                fillHisto("mjj_kfit_Id",cutflowType+"/Iso/KinFit", diJet.mass(), evtWeight);
                if(probOfKinFit > 0.1) fillHisto("mjj_kfit_Id_probfit1",cutflowType+"/Iso/KinFit", diJet.mass(), evtWeight);
                if(probOfKinFit > 0.2) fillHisto("mjj_kfit_Id_probfit2",cutflowType+"/Iso/KinFit", diJet.mass(), evtWeight);
                } //end-of  if(kfLightJets[0].pt() > 25 && kfLightJets[1].pt() > 25 )
              } // end-of if(match1 & match2) 
            } // end-of if(kfLightJets.size() >= 2)
          } // end-of foundMu
    	}//isomuon
  }//event loop
  
  cout<<"kfCount = "<<kfCount<<endl;
  outfile_ << "Number of times HadP and HadQ matches with jets" << matchjetcount << endl;
  outfile_ << "total number of selected events " << nSelEvents <<endl; 
  outfile_ << "No of times three pair jet matched:    " << threepairjet << endl;
  outfile_ << "No of correct jet pair assign to W:   " << matchedJet_q << endl; 
  outfile_ << "No of not correct jet pair assign to W:   " << not_matchedJet_q << endl;
  outfile_ << "No of correct jet pair not assign to W but |eta| < 2.5 and pt > 25:   " << matched_quark_eta_pt << endl;
  outfile_ << "No of correct jet matched to b (from top) :   " << matchedJet_b << endl;
  outfile_ << "Average weights for top-reweighting : "<<TotalTopPtWeights<<"/"<<TotalLplusJEvents<<" = "<<TotalTopPtWeights/TotalLplusJEvents<<endl; 
  outfile_ << "Efficiency Uncertainty on the SV finding : "<<TotalSVEff<<"/"<<TotalEvtsWithSV<<" = "<<TotalSVEff/TotalEvtsWithSV<<endl;
  fillHisto("AvTopPtWeight", cutflowType, 1.0, TotalTopPtWeights/TotalLplusJEvents);
  fillHisto("SVEffUncert", cutflowType, 2.0, TotalSVEff/TotalEvtsWithSV);
  fillHisto("SVEffUncert", cutflowType, 1.0, TotalEvtsWithSV);
  fillHisto("SVEffUncert", cutflowType, 0.0, TotalSVEff);
  f->Close(); 
  delete f;
}

void hplusAnalyzer::processEvents(){ 
  
  //Data, MC sample from lxplus and T2
  //CutFlowAnalysis("TTJets_ntuple_MuChannel.root", "PF", "TTJets_MuMC_check"); 
  
  //CutFlowAnalysis("root://se01.indiacms.res.in:1094//cms/store/user/rverma/ntuple_MuMC_kfitL_20170620/MuMC_20170620/TTJets_MuMC_20170620/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/TTJets_MuMC_20170620/170620_175939/0000/TTJets_MuMC_20170620_Ntuple_31.root", "PF", "ntupleToHistos");
  
		  //CutFlowAnalysis("root://se01.indiacms.res.in:1094//cms/store/user/rverma/ntuple_MuMC_kfitT_20170610/MuMC_20170608/TTJets_MuMC_20170608/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/TTJets_MuMC_20170608/170608_190531/0000/TTJets_MuMC_20170608_Ntuple_1.root", "PF", "TT");
  //CutFlowAnalysis("root://se01.indiacms.res.in:1094//cms/store/user/rverma/ntuple_MuMC_kfitT_20170610/MuMC_20170608/WZ_MuMC_20170608/WZ_TuneCUETP8M1_13TeV-pythia8/WZ_MuMC_20170608/170608_190850/0000/WZ_MuMC_20170608_Ntuple_3.root", "PF","DY1");
  
  //condor submission
  //CutFlowAnalysis("root://se01.indiacms.res.in:1094/inputFile", "PF", "outputFile");
  CutFlowAnalysis("/hadoop/cms/store/user/gkole/chargedHiggs/13TeV/Ntuple/June2017/TTJets_MuMC_20170620_Ntuple_1.root", "PF", "outputFile");
} 

float hplusAnalyzer::reweightHEPNUPWJets(int hepNUP) {

  int nJets = hepNUP-5;
  if(nJets==0)      return 2.11;
  else if(nJets==1) return 0.23;
  else if(nJets==2) return 0.119;
  else if(nJets==3) return 0.0562;
  else if(nJets>=4) return 0.0671;
  else return 1 ;
}

float hplusAnalyzer::reweightHEPNUPDYJets(int hepNUP){

  int nJets = hepNUP-5;
  if(nJets==0)      return 0.120;
  else if(nJets==1) return 0.0164;
  else if(nJets==2) return 0.0167;
  else if(nJets==3) return 0.0167;
  else if(nJets>=4) return 0.0128;
  else return 1 ;
}

