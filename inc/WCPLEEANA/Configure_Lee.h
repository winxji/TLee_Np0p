namespace config_Lee
{
  ////////// input files for spectra and covariance matrixes

  TString spectra_file = "./data_inputs/merge.root";
  TString flux_Xs_directory = "./data_inputs/hist_rootfiles/XsFlux/";
  TString detector_directory = "./data_inputs/hist_rootfiles/DetVar/";
  TString mc_directory = "./data_inputs/mc_stat/";

  //////////
  
  int syst_cov_flux_Xs_begin = 1;// files in flux_Xs_directory above
  int syst_cov_flux_Xs_end   = 17;
 
  int syst_cov_mc_stat_begin = 0;// files in mc_directory above
  //int syst_cov_mc_stat_end   = 99;
  int syst_cov_mc_stat_end   = 0; // this version should be for no LEE

  
  int channels_observation = 0;// data channels (=hdata_obsch_# in spectra_file above)
                               // which is equal to the channels after collapse
                               // NOTE: This value is not used in the lastest version

  ///////////////////////////////
  
  int array_LEE_ch[1] = {0}; // for 1d strength. Pay attention to the setting.
                             // Confict will happen if both the 1d and 2d(defined below) work.
                             // Three options: 1d works, 2d works, or neither works

  //int array_LEE_ch[12] = {2, 5, 8, 11, 14, 17, 3, 6, 9, 12, 15, 18};
  
  /////////////////////////////// separated Np and 0p strengthes: 2d strengthes

  int array_LEE_Np_ch[6] = {2, 5, 8, 11, 14, 17};// element value "0" will not be set to the LEE_ch
  int array_LEE_0p_ch[6] = {3, 6, 9, 12, 15, 18};// element value "0" will not be set to the LEE_ch

  //int array_LEE_Np_ch[1] = {0};
  //int array_LEE_0p_ch[1] = {0};
  
  ///////////////////////////////
  
  /// some places may need to be changed when use different file-formats
  /// void TLee::Set_Spectra_MatrixCov()
  /// (*) map_input_spectrum_ch_str      -----> prediction channels before collapse
  /// (*) map_Lee_ch                     -----> tag LEE channels
  /// (*) map_detectorfile_str           -----> detector files
   
  ////////// display graphics flag

  bool flag_display_graphics = 0;
  
  ////////// systematics flag
  
  bool flag_syst_flux_Xs    = 1;
  bool flag_syst_detector   = 1;
  bool flag_syst_additional = 1;
  bool flag_syst_mc_stat    = 1;

  //double Lee_strength_for_outputfile_covariance_matrix = 1;
  double Lee_strength_for_outputfile_covariance_matrix = 0;

  double Lee_Np_strength_for_outputfile_covariance_matrix = 0;
  double Lee_0p_strength_for_outputfile_covariance_matrix = 0;
  
  bool flag_plotting_systematics   = 0;
  
  ////////// goodness of fit
  
  double Lee_strength_for_GoF         = 0;
  
  double Lee_Np_strength_for_GoF      = 0;
  double Lee_0p_strength_for_GoF      = 0;

  bool flag_GoF_output2file_default_0 = 0;

  bool flag_lookelsewhere             = 0;
  
  bool flag_both_numuCC            = 0;// 1
  bool flag_CCpi0_FC_by_numuCC     = 0;// 2
  bool flag_CCpi0_PC_by_numuCC     = 0;// 3
  bool flag_NCpi0_by_numuCC        = 0;// 4
  bool flag_nueCC_PC_by_numuCC_pi0 = 0;// 5
  bool flag_nueCC_HghE_FC_by_numuCC_pi0_nueFC = 0;// 6, HghE>800 MeV
  bool flag_nueCC_LowE_FC_by_all   = 0;// 7
  bool flag_nueCC_FC_by_all        = 0;// 8

  ////////// Lee strength fitting -- data

  bool flag_Lee_strength_data = 0;
  bool flag_Lee_scan_data     = 0;

  bool flag_GOF = 0;
}
