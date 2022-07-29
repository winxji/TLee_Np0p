#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include "stdlib.h"
using namespace std;

#include<map>
#include<vector>

#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TLine.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "THStack.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"

#include "TRandom3.h"
#include "TGaxis.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TText.h"
#include "TLatex.h"

#include "TCanvas.h"
#include "TVirtualPad.h"
#include "TPad.h"
#include "TLegend.h"
#include "TString.h"
#include "TColor.h"

//////////////////////////////////////////////////////////////////////////////////

//
// How to run: 
//              root -l cal_CL.cc+
//

void cal_CL()
{
  TString str_file_data = "./fc_files/sub_fit_data.root ";
  TString str_file_distribution = "./fc_files/sub_fit_distribution.root";

  TFile *file_data = new TFile(str_file_data, "read");
  TFile *file_distribution = new TFile(str_file_distribution, "read");

  /////////////////////////////////

  // Declaration of leaf types
  Int_t           grid_Np;
  Int_t           grid_0p;
  Double_t        true_Np;
  Double_t        true_0p;
  vector<int>     *vec_min_status;
  vector<double>  *vec_chi2_var;
  vector<double>  *vec_min_chi2;
  vector<double>  *vec_dchi2;
  vector<double>  *vec_min_fNp_val;
  vector<double>  *vec_min_fNp_err;
  vector<double>  *vec_min_f0p_val;
  vector<double>  *vec_min_f0p_err;

  // List of branches
  TBranch        *b_grid_Np;   //!
  TBranch        *b_grid_0p;   //!
  TBranch        *b_true_Np;   //!
  TBranch        *b_true_0p;   //!
  TBranch        *b_vec_min_status;   //!
  TBranch        *b_vec_chi2_var;   //!
  TBranch        *b_vec_min_chi2;   //!
  TBranch        *b_vec_dchi2;   //!
  TBranch        *b_vec_min_fNp_val;   //!
  TBranch        *b_vec_min_fNp_err;   //!
  TBranch        *b_vec_min_f0p_val;   //!
  TBranch        *b_vec_min_f0p_err;   //!
  
  // Set object pointer
  vec_min_status = 0;
  vec_chi2_var = 0;
  vec_min_chi2 = 0;
  vec_dchi2 = 0;
  vec_min_fNp_val = 0;
  vec_min_fNp_err = 0;
  vec_min_f0p_val = 0;
  vec_min_f0p_err = 0;

  ///////

  TTree *tree_data = (TTree*)file_data->Get("tree");
  tree_data->SetBranchAddress("grid_Np", &grid_Np, &b_grid_Np);
  tree_data->SetBranchAddress("grid_0p", &grid_0p, &b_grid_0p);
  tree_data->SetBranchAddress("true_Np", &true_Np, &b_true_Np);
  tree_data->SetBranchAddress("true_0p", &true_0p, &b_true_0p);
  tree_data->SetBranchAddress("vec_min_status", &vec_min_status, &b_vec_min_status);
  tree_data->SetBranchAddress("vec_chi2_var", &vec_chi2_var, &b_vec_chi2_var);
  tree_data->SetBranchAddress("vec_min_chi2", &vec_min_chi2, &b_vec_min_chi2);
  tree_data->SetBranchAddress("vec_dchi2", &vec_dchi2, &b_vec_dchi2);
  tree_data->SetBranchAddress("vec_min_fNp_val", &vec_min_fNp_val, &b_vec_min_fNp_val);
  tree_data->SetBranchAddress("vec_min_fNp_err", &vec_min_fNp_err, &b_vec_min_fNp_err);
  tree_data->SetBranchAddress("vec_min_f0p_val", &vec_min_f0p_val, &b_vec_min_f0p_val);
  tree_data->SetBranchAddress("vec_min_f0p_err", &vec_min_f0p_err, &b_vec_min_f0p_err);

  TTree *tree_distribution = (TTree*)file_distribution->Get("tree");
  tree_distribution->SetBranchAddress("grid_Np", &grid_Np, &b_grid_Np);
  tree_distribution->SetBranchAddress("grid_0p", &grid_0p, &b_grid_0p);
  tree_distribution->SetBranchAddress("true_Np", &true_Np, &b_true_Np);
  tree_distribution->SetBranchAddress("true_0p", &true_0p, &b_true_0p);
  tree_distribution->SetBranchAddress("vec_min_status", &vec_min_status, &b_vec_min_status);
  tree_distribution->SetBranchAddress("vec_chi2_var", &vec_chi2_var, &b_vec_chi2_var);
  tree_distribution->SetBranchAddress("vec_min_chi2", &vec_min_chi2, &b_vec_min_chi2);
  tree_distribution->SetBranchAddress("vec_dchi2", &vec_dchi2, &b_vec_dchi2);
  tree_distribution->SetBranchAddress("vec_min_fNp_val", &vec_min_fNp_val, &b_vec_min_fNp_val);
  tree_distribution->SetBranchAddress("vec_min_fNp_err", &vec_min_fNp_err, &b_vec_min_fNp_err);
  tree_distribution->SetBranchAddress("vec_min_f0p_val", &vec_min_f0p_val, &b_vec_min_f0p_val);
  tree_distribution->SetBranchAddress("vec_min_f0p_err", &vec_min_f0p_err, &b_vec_min_f0p_err);

  ///////////////////////////////// data
  ///////////////////////////////// (1) find the real minimum chi2
  ///////////////////////////////// (2) map chi2 at each grid

  int entries_data = tree_data->GetEntries();
  cout<<endl<<" ---> entries_data: "<<entries_data<<endl<<endl;
  
  double data_min_chi2 = 1e8;
  double data_min_fNp_val = 0;
  double data_min_f0p_val = 0;
  double data_min_fNp_err = 0;
  double data_min_f0p_err = 0;
  
  map<int, map<int, double> >map_data_chi2_var;
  map<int, map<int, double> >map_data_dchi2;
  
  for(int ientry=0; ientry<entries_data; ientry++) {
    tree_data->GetEntry(ientry);

    int vec_size = vec_min_status->size();
    for(int isize=0; isize<vec_size; isize++) {
      if( vec_min_status->at(isize)==0 && (vec_min_chi2->at(isize) < data_min_chi2) ) {
	data_min_chi2 = vec_min_chi2->at(isize);
	data_min_fNp_val = vec_min_fNp_val->at(isize);
	data_min_f0p_val = vec_min_f0p_val->at(isize);
	data_min_fNp_err = vec_min_fNp_err->at(isize);
	data_min_f0p_err = vec_min_f0p_err->at(isize);	
      }

      map_data_chi2_var[grid_Np][grid_0p] = vec_chi2_var->at(isize);
      
    }// for(int isize=0; isize<vec_size; isize++)        
  }// for(int ientry=0; ientry<entries_data; ientry++)
  
  if( data_min_chi2==1e8 ) {
    cout<<endl<<" ERROR: No valid data min_chi2"<<endl<<endl;
  }
  
  for(auto it_Np=map_data_chi2_var.begin(); it_Np!=map_data_chi2_var.end(); it_Np++) {
    for(auto it_0p=it_Np->second.begin(); it_0p!=it_Np->second.end(); it_0p++) {
      grid_Np = it_Np->first;
      grid_0p = it_0p->first;

      map_data_dchi2[grid_Np][grid_0p] = map_data_chi2_var[grid_Np][grid_0p] - data_min_chi2;
    }
  }

  int number_grid_Np = map_data_dchi2.size();
  int number_grid_0p = map_data_dchi2.begin()->second.size();
  cout<<TString::Format(" ---> The parameters space is divided into %dx%d grids", number_grid_Np, number_grid_0p)<<endl;
  
  //cout<<endl;
  cout<<" ---> fitting: data_min_chi2 "<<data_min_chi2<<endl;
  cout<<" ---> fitting: data_min_fNp  "<<data_min_fNp_val<<" +/- "<<data_min_fNp_err<<" (error by default Minuit2)"<<endl;
  cout<<" ---> fitting: data_min_f0p  "<<data_min_f0p_val<<" +/- "<<data_min_f0p_err<<" (error by default Minuit2)"<<endl;  
  cout<<endl;

  TH2D *h2d_space = new TH2D("h2d_space", "x for index of Np and y for index of 0p", number_grid_Np, 0, number_grid_Np, number_grid_0p, 0, number_grid_0p);
  
  ///////////////////////////////// toy: distribution
  ///////////////////////////////// 
  ///////////////////////////////// 
  
  int entries_distribution = tree_distribution->GetEntries();
  //cout<<endl<<" ---> entries_distribution: "<<entries_distribution<<endl<<endl;

  map<int, map<int, vector<double> > >map_vec_toys_chi2_var;
  map<int, map<int, vector<double> > >map_vec_toys_min_chi2;
  map<int, map<int, vector<double> > >map_vec_toys_dchi2;

  for(int ientry=0; ientry<entries_distribution; ientry++) {
    tree_distribution->GetEntry(ientry);
    
    int vec_size = vec_min_status->size();
    for(int isize=0; isize<vec_size; isize++) {
      if( vec_min_status->at(isize)==0 ) {
	double chi2_val = vec_chi2_var->at(isize);
	double min_chi2 = vec_min_chi2->at(isize);
	double dchi2 = chi2_val - min_chi2;
	if( dchi2<0 ) {
	  cout<<" Invalid dchi2 less than 0: "<<dchi2<<endl;
	  continue;
	}

	map_vec_toys_chi2_var[grid_Np][grid_0p].push_back( chi2_val );
	map_vec_toys_min_chi2[grid_Np][grid_0p].push_back( min_chi2 );
	map_vec_toys_dchi2[grid_Np][grid_0p].push_back( dchi2 );
      }
    }// for(int isize=0; isize<vec_size; isize++)
    
  }// for(int ientry=0; ientry<entries_distribution; ientry++)

  ///////

  int min_valid_toys_at_grid = 10000000;
  int max_valid_toys_at_grid = 0;
  
  for(auto it_Np=map_vec_toys_dchi2.begin(); it_Np!=map_vec_toys_dchi2.end(); it_Np++) {
    for(auto it_0p=it_Np->second.begin(); it_0p!=it_Np->second.end(); it_0p++) {
      int vec_size = it_0p->second.size();
      if( min_valid_toys_at_grid>vec_size ) min_valid_toys_at_grid = vec_size;
      if( max_valid_toys_at_grid<vec_size ) max_valid_toys_at_grid = vec_size;

      grid_Np = it_Np->first;
      grid_0p = it_0p->first;
      
      int line_CL = 0;
      for(int isize=0; isize<vec_size; isize++) {
	double dchi2_toy = map_vec_toys_dchi2[grid_Np][grid_0p].at(isize);
	if( dchi2_toy<map_data_dchi2[grid_Np][grid_0p] ) line_CL++;
      }// for(int isize=0; isize<vec_size; isize++)

      double value_CL = line_CL*1./vec_size;
      
      h2d_space->SetBinContent(grid_Np, grid_0p, value_CL);
      
    }// for(auto it_0p=it_Np->second.begin(); it_0p!=it_Np->second.end(); it_0p++)
  }// for(auto it_Np=map_vec_toys_dchi2.begin(); it_Np!=map_vec_toys_dchi2.end(); it_Np++)

  cout<<" ---> minimum/maximum valid toys at grid: "<<min_valid_toys_at_grid<<"\t"<<max_valid_toys_at_grid<<endl<<endl;

  h2d_space->SetStats(0);
  h2d_space->Draw("colz");
  h2d_space->SaveAs("file_map_CL.root");
  
}
