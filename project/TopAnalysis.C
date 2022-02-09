#define TopAnalysis_cxx
// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.

#include "TopAnalysis.h"
#include "Tophistograms.h"
#include <iostream>
#include <cstring>
#include <string>

#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TMath.h>
#include <TLorentzVector.h>


/*

Int_t           jet_trueflav[10]
Float_t         lep_charge[4]


*/

string name;

void TopAnalysis::Begin(TTree * )
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).
}

void TopAnalysis::SlaveBegin(TTree * )
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();
  printf("Starting analysis with process option: %s \n", option.Data());

  name=option;

  define_histograms();

  FillOutputList();
}

Bool_t TopAnalysis::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // When processing keyed objects with PROOF, the object is already loaded
  // and is available via the fObject pointer.
  //
  // This function should contain the \"body\" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.
  //
  // The processing can be stopped by calling Abort().
  //
  // Use fStatus to set the return value of TTree::Process().
  //
  // The return value is currently not used.

  fChain->GetTree()->GetEntry(entry);
  //  int cut1_mc = 0;

  if(fChain->GetTree()->GetEntries()>0)
  {
    //Do analysis

    Float_t scaleFactor = scaleFactor_ELE*scaleFactor_MUON*scaleFactor_TRIGGER;   //SF
    Float_t eventWeight = mcWeight*scaleFactor_PILEUP*scaleFactor_ZVERTEX;        //EventW
    Double_t weight = scaleFactor*eventWeight;                                    //weight = SF * EventW

    // Make difference between data and MC
    if (weight == 0. && channelNumber != 110090 && channelNumber != 110091) weight = 1.;

    // Missing Et of the event in GeV
    Float_t missingEt = met_et/1000.;


    // PRESELECTION CUT
    if(trigE || trigM)        // Single electron or muon trigger is satisfied
    {
      if(passGRL)             // Event in real data passes the Good Run List
      {
        if(hasGoodVertex)     // good vertex
        {
          //Find the good leptons
          int goodlep_n = 0;
          int goodlep_index[lep_n];
          int lep_index = 0;

          for(int i=0; i<lep_n; i++)
          {
            if(lep_pt[i]>25000. && (lep_ptcone30[i]/lep_pt[i]) < 0.15 && (lep_etcone20[i]/lep_pt[i]) < 0.15 )  // isolation criteria
    	      {
          		// electron selection in fiducial region excluding candidates in the transition region between the barrel and endcap electromagnetic calorimeters
          		if ( lep_type[i]==11 && TMath::Abs(lep_eta[i]) < 2.47 && ( TMath::Abs(lep_eta[i]) < 1.37 || TMath::Abs(lep_eta[i]) > 1.52 ) ) {
          		  goodlep_n++;
                goodlep_index[lep_index] = i;  // NUEVA
                lep_index++;                     // NUEVA
          		}
          		if ( lep_type[i] == 13 && TMath::Abs(lep_eta[i]) < 2.5 ) {   // muon selection
          		  goodlep_n++;
                goodlep_index[lep_index] = i;  // NUEVA
                lep_index++;                     // NUEVA
          		}
    	      }
	        }


          //Zero cut
          if(goodlep_n==2)  // Exactly two good leptons with pT > 25 GeV //CAMBIO A 2
          {

            if(lep_type[goodlep_index[0]] == lep_type[goodlep_index[1]]) {            // same flavour
            if(lep_charge[goodlep_index[0]] * lep_charge[goodlep_index[1]] < 0 ) {    // different charges

            float mZ = 91.18;

            // TLorentzVector definitions
            TLorentzVector Lepton_1  = TLorentzVector();
            TLorentzVector Lepton_2  = TLorentzVector();
            TLorentzVector m_ll  = TLorentzVector();
            // TLorentzVector      MeT  = TLorentzVector();

            Lepton_1.SetPtEtaPhiE(lep_pt[goodlep_index[0]], lep_eta[goodlep_index[0]], lep_phi[goodlep_index[0]],lep_E[goodlep_index[0]]);
            Lepton_2.SetPtEtaPhiE(lep_pt[goodlep_index[1]], lep_eta[goodlep_index[1]], lep_phi[goodlep_index[1]],lep_E[goodlep_index[1]]);
            // MeT.SetPtEtaPhiE(met_et, 0, met_phi , met_et);


            //Calculation of the Invariant Mass using TLorentz vectors (First Lepton + MeT)
            m_ll = Lepton_1 + Lepton_2;
            float InvMass1       = m_ll.M();
            float InvMass1_inGeV = InvMass1/1000.;


            if(TMath::Abs(InvMass1_inGeV - mZ) < 20) {   // mll - mZ

    	      if(jet_n == 0)  {                            // no jets required

  			    if(missingEt >= 0.)  {                       // all missing ET

					double dphi=Lepton_1.Phi() - Lepton_2.Phi();
      				double names_of_global_variable[]={InvMass1_inGeV, missingEt,dphi};
              // double names_of_global_variable[]={InvMass1_inGeV, missingEt, Lepton_1.Eta()};
      				double names_of_leadlep_variable[]={Lepton_1.Pt()/1000., Lepton_1.Eta(), Lepton_1.E()/1000., Lepton_1.Phi(), lep_charge[goodlep_index[0]], (double)lep_type[goodlep_index[0]], lep_ptcone30[goodlep_index[0]], lep_etcone20[goodlep_index[0]], lep_z0[goodlep_index[0]], lep_trackd0pvunbiased[goodlep_index[0]]};
              double names_of_subleadlep_variable[]={Lepton_2.Pt()/1000., Lepton_2.Eta(), Lepton_2.E()/1000., Lepton_2.Phi(), lep_charge[goodlep_index[1]], (double)lep_type[goodlep_index[1]], lep_ptcone30[goodlep_index[1]], lep_etcone20[goodlep_index[1]], lep_z0[goodlep_index[1]], lep_trackd0pvunbiased[goodlep_index[1]]};

      				TString histonames_of_global_variable[]={"hist_vismass","hist_etmiss", "hist_d_phi_ll"};
      				TString histonames_of_leadlep_variable[]={"hist_leadleptpt", "hist_leadlepteta","hist_leadleptE","hist_leadleptphi","hist_leadleptch","hist_leadleptID","hist_leadlept_ptc","hist_leadleptetc","hist_leadlepz0","hist_leadlepd0"};
              TString histonames_of_subleadlep_variable[]={"hist_subleadleptpt", "hist_subleadlepteta","hist_subleadleptE","hist_subleadleptphi","hist_subleadleptch","hist_subleadleptID","hist_subleadlept_ptc","hist_subleadleptetc","hist_subleadlepz0","hist_subleadlepd0"};

      				int length_global = sizeof(names_of_global_variable)/sizeof(names_of_global_variable[0]);
      				int length_leadlep = sizeof(names_of_leadlep_variable)/sizeof(names_of_leadlep_variable[0]);
              int length_subleadlep = sizeof(names_of_subleadlep_variable)/sizeof(names_of_subleadlep_variable[0]);

      				for (int i=0; i<length_global; i++)
      				  {
      				    FillHistogramsGlobal( names_of_global_variable[i], weight, histonames_of_global_variable[i]);
      				  }
      				for (int i=0; i<10; i++)
      				  {
      				    FillHistogramsLeadlept( names_of_leadlep_variable[i], weight, histonames_of_leadlep_variable[i]);
      				  }
              for (int i=0; i<10; i++)
      				  {
      				    FillHistogramsSubLeadlept( names_of_subleadlep_variable[i], weight, histonames_of_subleadlep_variable[i]);
      				  }
			       }
			      }
	         }
	        }
	       }
	      }
	     }
      }
    }}
  //  std::cout<<cut1_mc<<std::endl;
  return kTRUE;
}

void TopAnalysis::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.
}

void TopAnalysis::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  name="output_Top/"+name+".root";

  const char* filename = name.c_str();

  TFile physicsoutput_Top(filename,"recreate");
  WriteHistograms();
  physicsoutput_Top.Close();


}
