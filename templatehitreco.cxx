//! \file template_code6d.cpp
//!
//! Template hit reconstruction algorithms, add FPix to templates, add double pixels to templates, 
//! change dp calling seq (Add PSI46v2 response function and allow for zero-suppressed ROC output)
//! Change to Root historgraams for long-term compatibility
//! Add angle vs resolution for templates and "Standard Algorithms"
//! Tune CMSSW simulation for template 1 reconstruction
//! Change standard algorithm to always use edge method for y-reconstruction
//! Add Estar template number 4
//! Do cosmics
//! Add clustering  algorithm
//! Change response function, parameters to that used for real analysis


#include <stdlib.h> 
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include "boost/multi_array.hpp"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <sys/time.h>
#include "SiPixelTemplate.cc"
static int theVerboseLevel = {2};
#include "SiPixelTemplateReco.cc"
#include "VVIObjF.cc"

using namespace std;

#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TObject.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPostScript.h"

//--- FastSim resolution histograms:
#include "PixelResolutionHistograms.h"


// Global definitions 
const float cmtomicron = 10000.0;
const int   debug      = 0;

// Main program  

int main(int argc, char *argv[])
{
  // *****************************************
  //  Global config for this macro

  // New electronic response (phase 1), all layers except layer 1.  Not used
  // if we selected linear response below.
  const double gain = 3.19;
  const double ped = 16.46;
  const double p0 = 0.01218;
  const double p1 = 0.711;
  const double p2 = 203.;
  const double p3 = 148.;	
  static double vcal = 47.;	
  static double vcaloffst = 60.;

  // Layer 1 in phase 1 needs special handling.
  int layer = 2; 
  // If this is Phase 1, Layer 1, change Vcal settings.
  if (layer == 1) {
    vcal = 50.;	
    vcaloffst = 670;
  } 


  // Local variables 
  int Ntot = 0;
  printf("number 3 %s",argv[3]);
  Ntot = strtol(argv[3], NULL, 10);;
  if (argv[3]==NULL) {
    Ntot = 10000;           // exit after this many events   *** CHANGE TO RUN OVER ALL EVENTS                                                                                 
  }
  // FIXME: Read total number of lines e.g. cat $TEMPLATE_EVENTS_NFILE | wc -l


  // *****************************************
  
  int i_nlines_debug = 0;
  float max_cotalpha = -1.0;
  float max_cotbeta  = -1.0;
  float min_cotalpha = 99.0;
  float min_cotbeta  = 99.0;

  std::vector<float> pvec(6), wgauss(TYSIZE), vgauss(TYSIZE), xgauss(TYSIZE), ygauss(TYSIZE), zgauss(TYSIZE);
  float pixin[TXSIZE][TYSIZE];
  bool ydouble[TYSIZE], xdouble[TXSIZE];
  float thick, xsize, ysize, noise, zcen, gain_frac, q100_frac, common_frac, readout_noise, qscale, qperbit, pixmax;
  float xhit, yhit, xrec, yrec, sigmax, sigmay, probx, proby, signal, cotalpha, cotbeta, qclust, locBz, locBx, probxy, probQ;  
  int nfile, neh, nevent, ID, non_linear, nbad, icol, ndcol, iqbin;
  vector<int> nbin(5,0);
  int i, j, k, ierr, qbin, nypix(0), nxpix(0), etabin, jmin, jmax, imin, imax, numadd, idcol;
  double dx, dy, eta, log10probxy, log10probQ, tote, bade, weight, alpha, adc;
  int iyd, ixd, speed, nbits;
  float q100, q101, q50, q10, qmax;
  int mrow = TXSIZE, mcol = TYSIZE;
  float sigtmp, qin;
  char infile[80], header[80], c, dtitle[80], outfile[80], outfile0[80], outfile1[80], outfile2[80];
  int triplg(std::vector<float>&);
  float cluster[TXSIZE][TYSIZE], clust[TXSIZE][TYSIZE], rclust[TXSIZE][TYSIZE];
  bool bclust[TXSIZE][TYSIZE];
  std::pair<int, int> pixel, max;
  const bool deadpix = false;
  std::vector<std::pair<int, int> > zeropix;
  unsigned int detType(0);

  // File to open, C-style
  FILE *ifp;
  
  struct timeval now0, now1;
  struct timezone timz;
  long deltas, deltaus;
  double deltat;
  
  
  // Basic controls are via phase2.txt file.  The best is to have a copy for each application,
  // and then use ln -sf to make a soft link to phase2.txt.
  //
  // Below, we will read which data and inputs to use.  The file
  // phase2.txt contains the configuration of electronics.
  //
  // ~% cat phase2_forward_50x50.txt 
  // 49291 289 150. 1000. 1000. 0.000 0.000 0.00 000. 0 6
  //
  // #1: run number of the PIXELAV file, in this case template_events_d49291.out.
  //     It must match the template!  For instance, this is what's in this file:
  // zavist~/FastSim_Template_Code] head -1 template_events_d49291.out
  //     Dot1_50x50_prd2010lk:2010lk@-150V,3.8T@00d,263K,rhe=1.10,rhh=0.7   <-- this has 50x50 pixels!
  //
  // #2: template number (here 289)
  // zavist~/FastSim_Template_Code] head -1 template_summary_zp0289.out
  //     Dot1_50x50_prd2010lk:2010lk@-150V,3.8T@00d,263K,rhe=1.10,rhh=0.7
  
  // #3: 150 = preamp noise (applied to threshold)
  // #4 and #5: 1000 and 1000 = readout thresholds for single hit, and multiple hit per doublecolumn
  //    (probably will not apply to phase 2)
  // next 4 arguments (zeros here): threshold variations (noise applied after threshold)
  // #9:  0 = linear or not (0 is linear)
  // #10: 1 = phase1 realistic                 *** WILL NEED TO STUDY
  // #11: 6 = number of bits                   *** WILL NEED TO STUDY

  
  // (use c file i/o which is much faster than c++ i/o)
  const std::string infilename  = argv[1];
  ifp = fopen(infilename.c_str(), "r");
  if (ifp==NULL) {
    printf("no pixel initialization file found/n");
    return 0;
  }
  printf("Opening %s\n",infilename.c_str());
  
  fscanf(ifp,"%d %d %f %f %f %f %f %f %f %d %d", &nfile, &ID, &noise, &q100, &q101, &q100_frac, &common_frac, &gain_frac, &readout_noise, &non_linear, &nbits);
  //nbits = 4;
  fclose(ifp);
  printf("mc file %d noise = %f, threshold0 = %f, threshold1 = %f, rms threshold frac = %f, common_frac = %f, gain fraction = %f, readout noise = %f, linear response = %d, number of adc bits = %d \n", nfile, noise, q100, q101, q100_frac, common_frac, gain_frac, readout_noise, non_linear, nbits);
  
  //  Calculate 50% of threshold in q units and enc noise in adc units
  q50=0.5*q100;
  q10=0.2*q50;
  

  //  Create an input filename for this run 
  //
  if(nfile < 10000) {
    sprintf(infile,"./template_events_d%4.4d.out",nfile);
  } else {
    sprintf(infile,"./template_events_d%5.5d.out",nfile);
  }
  printf(infile,"%s \n");

  
  //  Open input file and read header info 
  //
  ifp = fopen(infile, "r");
  if (ifp==NULL) {
    printf("no pixel data file found\n");
    return 0;
  }
  
  // Read-in a header string first and print it    
  //
  for (i=0; (c=getc(ifp)) != '\n'; ++i) {
    if(i < 79) {header[i] = c;}
  }
  if(i > 78) {i=78;}
  header[i+1] ='\0';
  printf("%s\n", header);
  
  double  halfxs=300.;
  int nx=120;	
  gStyle->SetOptFit(101);
  gStyle->SetHistLineWidth(2);
  static vector<TH1F*> hp(30);
  hp[0] = new TH1F("h101","dy_temp (all sig); #Deltay (#mum)",nx,-halfxs,halfxs);
  hp[1] = new TH1F("h102","dy_temp (signal @> 1.5mn); #Deltay (#mum)",nx,-halfxs,halfxs);      
  hp[2] = new TH1F("h103","dy_temp (1.5mn @> signal @> 1.0mn); #Deltay (#mum)",nx,-halfxs,halfxs);      
  hp[3] = new TH1F("h104","dy_temp (1.0mn @> signal @> 0.85mn); #Deltay (#mum)",nx,-halfxs,halfxs);     
  hp[4] = new TH1F("h105","dy_temp (0.85mn @> signal); #Deltay (#mum)",nx,-halfxs,halfxs);      
  hp[5] = new TH1F("h201","log10(probxy) (all sig)",nx,-12.,0.);
  hp[6] = new TH1F("h202","log10(probxy) (signal @> 1.5mn)",nx,-12.,0.);      
  hp[7] = new TH1F("h203","log10(probxy) (1.5mn @> signal @> 1.0mn)",nx,-12.,0.);      
  hp[8] = new TH1F("h204","log10(probxy) (1.0mn @> signal @> 0.85mn)",nx,-12.,0.);     
  hp[9] = new TH1F("h205","log10(probxy) (0.85mn @> signal)",nx,-12.,0.);      
  hp[10] = new TH1F("h106","dx_temp (all sig); #Deltax (#mum)",nx,-halfxs,halfxs);
  hp[11] = new TH1F("h107","dx_temp (signal @> 1.5mn); #Deltax (#mum)",nx,-halfxs,halfxs);      
  hp[12] = new TH1F("h108","dx_temp (1.5mn @> signal @> 1.0mn); #Deltax (#mum)",nx,-halfxs,halfxs);      
  hp[13] = new TH1F("h109","dx_temp (1.0mn @> signal @> 0.85mn); #Deltax (#mum)",nx,-halfxs,halfxs);      
  hp[14] = new TH1F("h110","dx_temp (0.85mn @> signal); #Deltax (#mum)",nx,-halfxs,halfxs);    
  hp[15] = new TH1F("h206","log10(probQ) (all sig)",nx,-12.,0.);
  hp[16] = new TH1F("h207","log10(probQ) (signal @> 1.5mn)",nx,-12.,0.);      
  hp[17] = new TH1F("h208","log10(probQ) (1.5mn @> signal @> 1.0mn)",nx,-12.,0.);      
  hp[18] = new TH1F("h209","log10(probQ) (1.0mn @> signal @> 0.85mn)",nx,-12.,0.);     
  hp[19] = new TH1F("h210","log10(probQ) (0.85mn @> signal)",nx,-12.,0.);      
  hp[20] = new TH1F("h300","cotbeta",nx,0.1,5.7);
  hp[21] = new TH1F("h301","cotalpha",nx,0.1,0.8);
  //hp[20] = new TH1F("h300","cotbeta",nx,0.1,0.8);
  //hp[21] = new TH1F("h301","cotalpha",nx,0.1,0.8);
  //hp[20] = new TH1F("h300","cotbeta",nx,0.,8);
  //hp[21] = new TH1F("h301","cotalpha",nx,-0.42,0.42);
  hp[22] = new TH1F("h111","dy_temp (one pixel clust); #Deltay (#mum)",nx,-halfxs,halfxs);
  hp[23] = new TH1F("h112","dx_temp (one pixel clust); #Deltax (#mum)",nx,-halfxs,halfxs);
  hp[24] = new TH1F("h501","xy Probability; probxy",101,0,1.01);
  hp[25] = new TH1F("h502","Q Probability; probQ",101,0,1.01);
  hp[26] = new TH1F("h113","dy_temp (>one pixel clust); #Deltay (#mum)",nx,-halfxs,halfxs);
  hp[27] = new TH1F("h114","dx_temp (>one pixel clust); #Deltax (#mum)",nx,-halfxs,halfxs);
  
  // extra
  hp[28] = new TH1F("h302","cotbeta",nx,-5,5);
  hp[29] = new TH1F("h303","cotalpha",nx,-5,5);

  // Set style for the the histograms	
  
  for(i=0; i<30; ++i) {
    hp[i]->SetLineColor(2);
    hp[i]->SetFillColor(38);
  }
  
  


  fscanf(ifp,"%f  %f  %f", &ysize, &xsize, &thick);
  zcen = thick/2.;
  printf("xsize/ysize/thick = %f/%f/%f \n", xsize, ysize, thick);
  
  nevent=0;
  nbad = 0;
  tote = 0.;
  bade = 0.;
  
  static vector<float> sx(5,0.), sx2(5,0.), scx(5,0.), sy(5,0.), sy2(5,0.), scy(5,0.); 
  static vector<float> sxp(5,0.), sxp2(5,0.), scxp(5,0.), syp(5,0.), syp2(5,0.), scyp(5,0.); 
  static vector<float> sxc(5,0.), scxc(5,0.), sxc2(5,0.), syc(5,0.), scyc(5,0.), syc2(5,0.), nt(16,0.), nb(16,0.);
  std::vector<std::pair<int, int> > pixlst;
  
  // Decide if this file corresponds to a CMSSW run
  //
  /// printf("Enter ID of Template \n");
  /// scanf("%d", &ID);
  printf("Template ID = %d, barrel \n", ID);
  // if barrel: detType = 1;
  detType = 1;
  	
  // Decide if this file corresponds to a CMSSW run 


  // Initialize template store 
  
  std::vector< SiPixelTemplateStore > thePixelTemp_;
  SiPixelTemplate templ(thePixelTemp_);
  
  // Initialize template store, Pixelav 100V/300V simulation, +20C as thePixelTemp[6] 
	
  templ.pushfile(ID,thePixelTemp_);
  templ.interpolate(ID, 0.f, 0.f, -1.f);
  qscale = templ.qscale();
  pixmax = templ.pixmax(); // any charge > pixmax gets truncated
  // 2*pixmax as max scale of adc bit
  // other charge information does not provide more
  qperbit = pixmax/(pow(2.,(double)(nbits)-1.)); 
  // all of the bits up to pixmax
  // 4 or 5 bits for new chip
  // 2 bits above pedestal
  // 4 or 5 bits --> res
  
  // Ask for speed info
  
  //printf("enter algorithm speed (-2->5)\n");
  //scanf("%d", &speed);
  speed = 0;
  
  printf("speed = %d\n", speed);
  
  // Ask for external input into which pixels to join 
  //
  //printf("enter the first y- and x-columns to double (-1 -1 = none)\n");
  //scanf("%d %d", &iyd, &ixd);
  iyd = -1; ixd = -1;
  
  printf("y/x double = %d/%d\n", iyd,ixd);
	
  int totale=0;  int goodx=0; int goody=0;
   
  ndcol = TYSIZE/2;
  std::vector<int> ndhit(ndcol, 0);
  
  //  Determine current time
  gettimeofday(&now0, &timz);


  // Resolution binning
  // max_cotalpha=0.775855  max_cotbeta=0.775861min_cotalpha=0.10345  min_cotbeta=0.15001 Forward Phase1
  // max_cotalpha=0.4  max_cotbeta=9.47363min_cotalpha=0  min_cotbeta=6.00094e-05 Barrel Phase1
  // max_cotalpha=0.133333  max_cotbeta=0.666667min_cotalpha=0  min_cotbeta=1.00055e-06
  // barrel phase 2
  //max_cotalpha=0.416667  max_cotbeta=8.00002min_cotalpha=0  min_cotbeta=2.05427e-06

  const double	cotbetaBinWidth = 0.25; // 0.25 forward2 //0.21 //0.15 // 3.3. barrel, 0.21 forward1, 0.2 forward0
  const double	cotbetaLowEdge	= 0.0 ; // 0.0 forward2 //0.15 forward1, 0.25 forward0
  const int	cotbetaBins	= 3; //3
  const double    cotalphaBinWidth = 0.10; // 0.10 forward2 //0.38//0.02 //0.2 barrel 0.34 forward1, 0.28 forward0
  const double    cotalphaLowEdge = 0.0; // -0.04, 0.10 forward0
  const int       cotalphaBins    = 2; //2
  // barrel
  // cotbeta -8 8, cotalph 0.3 0.45
  //const double	cotalphaBinWidth = 0.02;
  //const double	cotalphaLowEdge	= -0.04;
  ///const int	qbinWidth	= 1;
  ///const int	qbins		= 4;
  
  // Output filename
  sprintf(outfile,"pixel_histos%5.5d_%3.3d_%i.root",nfile,ID,nbits);

  // Descriptive title
  // e.g. Forward 50x50x150 flat disk pixel resolution histograms
  sprintf(dtitle,"%s",argv[2]);
  printf("dtitle %s\n",argv[2]);
  if (argv[2]==NULL) {
    sprintf(dtitle,"Forward pixel resolution histograms");
  }

  //--- Create the FastSim ResolutionHistograms storage object.
  //    Internally it will book all resolution histograms.
  //    It also creates a dummy 2D histogram to store the binning
  //    along cotbeta (x-axis) and cotalpha (y-axis).  (The qBin always
  //    has the same binning.)
  //
  PixelResolutionHistograms
    fastSimResHistoStore( outfile,                                // File name for histograms
			  "",                                     // No subdirectory
			  dtitle,                                 // Descriptive title	     
			  detType, // unsigned int detType,             // Do we need this?
			  cotbetaBinWidth, cotbetaLowEdge, cotbetaBins,    // Binning along cot\beta
			  cotalphaBinWidth, cotalphaLowEdge, cotalphaBins); // ... along cot\alpha
    //                    qbinWidth, qbins );                              // ... for qBin

  
  //
  // ==============================================================
  // Loop over PIXELAV hits (until end of input file) 
  //
  //
  while(fscanf(ifp,"%f %f %f %f %f %f %d", &pvec[0], &pvec[1], &pvec[2], &pvec[3], &pvec[4], &pvec[5], &neh) != EOF) {
    
    i_nlines_debug++;
    if (i_nlines_debug > Ntot)
      break;
    
    // Read the input PIXELAV-simulated hits
    //
    for (k=0; k < TXSIZE; ++k) {
      fscanf(ifp,
	     "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f", 
	     &pixin[k][0],&pixin[k][1],&pixin[k][2],&pixin[k][3],&pixin[k][4],&pixin[k][5],&pixin[k][6],&pixin[k][7],&pixin[k][8],&pixin[k][9],
	     &pixin[k][10],&pixin[k][11],&pixin[k][12],&pixin[k][13],&pixin[k][14],&pixin[k][15],&pixin[k][16],&pixin[k][17],&pixin[k][18], 
	     &pixin[k][19],&pixin[k][20]);
    }
    ++nevent;
    if (nevent % 1000 == 0) {
      std::cout << "Processing event: " << nevent << std::endl;
    }
    if (debug > 1)
      cout << "\n--------------------------- \nNew event" << endl;

    
    // Add noise and analog response to cluster, reformat for flipped barrel coordinate system
    // Nose: before and after threshold, some of it is correlated.  (Preamp noise is not very relevant.
    // some comes from pixel-to-pixel variations not calibrated in CMSSW.)
    //
    triplg(vgauss);
    pixlst.resize(0);
    for(i=0; i<ndcol; ++i) {ndhit[i] = 0;}
    icol = 0;
    if(vgauss[1] < 0.) {icol = 1;}
    for(j=0; j<TXSIZE; ++j) {
      triplg(wgauss);
      triplg(xgauss);
      triplg(ygauss);
      triplg(zgauss);
      for(i=0; i<TYSIZE; ++i) {
	bclust[j][i] = false;
	qin = (10.*pixin[j][i] + xgauss[i]*noise);
	rclust[TXSIZE-1-j][TYSIZE-1-i] = qin;
	if(qin < q100*(1.+wgauss[i]*q100_frac)) {
	  clust[TXSIZE-1-j][TYSIZE-1-i] = 0.;
	} else {
	  idcol = (TYSIZE-1-i+icol)/2;
	  ++ndhit[idcol];
	  // hyptg response function
	  if(non_linear == 0) {
	    qin *= (1.+gain_frac*ygauss[i]);
	    signal = (qin + zgauss[i]*readout_noise)/qscale;
	    int ibits = (int)((qin + zgauss[i]*readout_noise)/qscale/qperbit); //charge/chargeperbit + some noise 
	    if (ibits > pow(2,double(nbits)) - 1){ ibits = pow(2,double(nbits)) -1;}
	    signal = qperbit*((float)ibits + 0.5); // don't truncate down to 0 -- saturate max charge
	  } else {
	    adc = (double)((int)(p3+p2*tanh(p0*(qin + vcaloffst)/(7.0*vcal) - p1)));
	    signal = ((float)((1.+gain_frac*ygauss[i])*(vcal*gain*(adc-ped))) - vcaloffst + zgauss[i]*readout_noise)/qscale;
	  }	 
	  clust[TXSIZE-1-j][TYSIZE-1-i] = (1.+vgauss[0]*common_frac)*signal;
	}
      }
    }
    
    // Single pixel per double column need to have a higher readout
    // threshold, so we simulate it here.  (Need more charge in single pixel
    // to get noticed by the token bit manager.)
    //
    for(j=0; j<TXSIZE; ++j) {
      for(i=0; i<TYSIZE; ++i) {
	if(clust[j][i] > 0.) {
	  idcol = (i+icol)/2;
	  if(ndhit[idcol] == 1) {
	    // Apply higher threshold on single hits in double columns
	    if(rclust[j][i] < q101*(1.+wgauss[i]*q100_frac)) {
	      clust[j][i] = 0.;
	    }
	  }
	}
      }
    }
    

    // Beginning of clustering (collecting pixel charges into a cluster).
    // This reproduces the algorithm we run in CMSSW.

    // Step 1: find the seed pixel.
    //
    qmax = 0.;
    for(j=0; j<TXSIZE; ++j) {
      for(i=0; i<TYSIZE; ++i) {
	if(clust[j][i] > qmax) {
	  qmax = clust[j][i];
	  max.first = j; max.second = i;
 	}
      }
    }
    if(qmax < 1500.) {
      if (debug > 5) 
	cout << "Qmax=" << qmax << " is below threshold, continuing to next cluster." << endl;
      continue;   // This PIXELAV hit would fail clustering -- skip to next one.
      //          // This was 5000 in phase0 and phase1, but changed for phase 2 because 
    }
    
    
    // OK, now we have the seed pixel.
    // Step 2: Simulate clustering around maximum pixel (seed), whose location is
    // kept in "max" (which is pair of ints).
    //
    pixlst.push_back(max);
    j=max.first; i=max.second;
    bclust[j][i] = true;

    if (debug > 5) 
      cout << "Seed pixel j=" << max.first << " i=" << max.second << endl;
    
    
    std::vector<std::pair<int, int> >::const_iterator pixIter, pixEnd;

    // "rescn" is label, we'll jump back here if we don't exit the loop.
    // (NB: this could have been written with a while() or repeat-until...)
  rescn: pixIter = pixlst.begin();
    pixEnd = pixlst.end();
    numadd = 0;
    for ( ; pixIter != pixEnd; ++pixIter ) {
      jmin = pixIter->first-1; 
      jmax = jmin+3;
      if(jmin < 0) {jmin = 0;}
      if(jmax > TXSIZE) {jmax = TXSIZE;}
      imin = pixIter->second-1;
      imax = imin+3;
      if(imin < 0) {imin = 0;}
      if(imax > TYSIZE) {imax = TYSIZE;}
      for(j=jmin; j<jmax; ++j) {
	for(i=imin; i<imax; ++i) {
	  if(clust[j][i] > 0.) {
	    if(!bclust[j][i]) {
	      bclust[j][i] = true;
	      pixel.first = j; pixel.second = i;
	      pixlst.push_back(pixel);
	      ++numadd;
	    }
	  }
	}
      }
    }
    if (numadd > 0) goto rescn;   // If have pixels to add, jump back up.

    //
    for(j=0; j<TXSIZE; ++j) {
      for(i=0; i<TYSIZE; ++i) {
	cluster[j][i] = 0.;
      }
    }
    pixIter = pixlst.begin();
    pixEnd = pixlst.end();
    qclust=0.;
    for ( ; pixIter != pixEnd; ++pixIter ) {
      j = pixIter->first; 
      i = pixIter->second;
      cluster[j][i] = clust[j][i];
      qclust += clust[j][i];
    }
    
    // Calculate the hit coordinates in the flipped coordinate system 
    
    yhit = -(pvec[0] + (zcen-pvec[2])*pvec[3]/pvec[5]);
    xhit = -(pvec[1] + (zcen-pvec[2])*pvec[4]/pvec[5]);
    
    // Calculate local particle angles, cotalpha, and cotbeta
    
    cotalpha = pvec[4]/pvec[5];
    alpha = atan((double)cotalpha);
    cotbeta = pvec[3]/pvec[5];
    //		 if(cotbeta < 0.) continue;
    eta = fabs(-log((double)(-cotbeta+sqrt((double)(1.+cotbeta*cotbeta)))));
    weight = 1./cosh((double)eta);
    tote += weight;
    etabin = (int)(eta/0.25);
    if(etabin > 15) {etabin = 15;}
    ++nt[etabin];

    // Keep track of the maximal abs cotalpha and cotbeta
    if (fabs(cotalpha) > max_cotalpha)
      max_cotalpha = fabs(cotalpha);
    if (fabs(cotbeta) > max_cotbeta)
      max_cotbeta = fabs(cotbeta);
    if (fabs(cotalpha) < min_cotalpha)
      min_cotalpha = fabs(cotalpha);
    if (fabs(cotbeta) < min_cotbeta)
      min_cotbeta = fabs(cotbeta);
    
    //if (debug > 5) {
    //cout << "cotalpha=" << cotalpha << "  cotbeta=" << cotbeta << endl;
    //cout << "xhit=" << xhit << "    yhit=" << yhit << endl;
    //}


#if 0    
    // Combine two single pixels into a double pixel. (Just to simulate
    // it... should be disabled code... ?)
    //
    // Along Y:
    //
    for(i=0; i<TYSIZE; ++i) {
      ydouble[i] = false;           // turning it off...
    }
    if(iyd >= 0 && iyd < TYSIZE) {
      ydouble[iyd] = true;
      for(j=0; j<TXSIZE; ++j) {
	sigtmp = cluster[j][iyd]+cluster[j][iyd+1];
	cluster[j][iyd] = sigtmp;
      }
      for(i=iyd+1; i<TYSIZE-1; ++i) {
	for(j=0; j<TXSIZE; ++j) {
	  cluster[j][i] = cluster[j][i+1];
	}
      }    
      for(j=0; j<TXSIZE; ++j) {
	cluster[j][TYSIZE-1] = 0.;
      }
    }
    //
    // Along X:
    //
    for(j=0; j<TXSIZE; ++j) {
      xdouble[j] = false;           // turning it off...
    }
    if(ixd >= 0 && ixd < TXSIZE-1) {
      xdouble[ixd] = true;
      for(i=0; i<TYSIZE; ++i) {
	sigtmp = cluster[ixd][i]+cluster[ixd+1][i];
	cluster[ixd][i] = sigtmp;
      }
      for(j=ixd+1; j<TXSIZE-1; ++j) {
	for(i=0; i<TYSIZE; ++i) {
	  cluster[j][i] = cluster[j+1][i];
	}
      }    
      for(i=0; i<TYSIZE; ++i) {
	cluster[TXSIZE-1][i] = 0.;
      }
    }
#endif

    
    // We have simulated the clustering, and now we have the cluster.
    // We can now run the template reconstruction on the cluster.
    //
    SiPixelTemplateReco::ClusMatrix clusterPayload{&cluster[0][0], xdouble, ydouble, mrow,mcol};

    // Flip coordinate systems so that we increase template/res histo statistics,
    // by only looking at cotbeta>0 and cotalpha>0.
    //
    locBx = 1.;
    if(cotbeta < 0.) locBx = -1.;
    locBz = locBx;
    if(cotalpha < 0.) locBz = -locBx;
    //
    // Template reconstruction:
    //
    ierr = PixelTempReco2D(ID, cotalpha, cotbeta, locBz, locBx, clusterPayload, templ, yrec, sigmay, proby, xrec, sigmax, probx, qbin, speed, deadpix, zeropix, probQ, nypix, nxpix);		 
    
    if(ierr != 0) {
      // Template reconstruction has failed.
      ++nbad; ++nb[etabin]; bade +=weight;
      printf("reconstruction failed with error %d \n", ierr);
    }
    else {
      // Template reconstruction succeeded.
      
      // Check resolution and weights

      // First, figure out which qBin is appropriate for this cluster.
      if (debug > 1)
	cout << "qbin=" << qbin << endl;

      if(qbin < 0) qbin = 0;
      iqbin = qbin;
      if(qbin > 3) {
	qbin = 4;
	iqbin = 3;
	//printf(" qbin = %d \n", qbin);
      }
      ++nbin[qbin];

      // Calculate dy shifts between the reconstructed hit and the true hit.  This involves
      // a change in convention where the hit is, which is different in PIXELAV and CMSSW.
      //
      if(iyd != 0) {
	dy = yrec - (TYSIZE/2)*ysize - yhit;
	if (debug > 10)
	  cout << "TYSIZE/2*ysize=" << (TYSIZE/2)*ysize << endl;
      } else {
	dy = yrec - ((TYSIZE/2)-0.5)*ysize - yhit;
      }
      ++totale;
      probxy = probx*proby*(1.-log((double)(probx*proby)));
      if(probxy > 1.e-2) {++goodx;}
      if(probxy > 1.e-2) {++goody;}
      hp[20]->Fill((double)cotbeta);
      hp[21]->Fill((double)cotalpha);

      hp[28]->Fill((double)cotbeta);
      hp[29]->Fill((double)cotalpha);

      log10probxy = log10((double)probxy); log10probQ = log10((double)probQ);
      sy[qbin] += dy; sy2[qbin] += dy*dy; scy[qbin] += dy*dy/(sigmay*sigmay);
      if(probxy > 1.e-2) {syp[qbin] += dy; syp2[qbin] += dy*dy; scyp[qbin] += dy*dy/(sigmay*sigmay);}
      hp[0]->Fill(dy);
      hp[1+iqbin]->Fill(dy);
      hp[5]->Fill(log10probxy);
      hp[6+iqbin]->Fill(log10probxy);

      // Calculate dx shifts between the reconstructed hit and the true hit.  This involves
      // a change in convention where the hit is, which is different in PIXELAV and CMSSW.
      //
      if(ixd != 0) {
	dx = xrec - (TXSIZE/2)*xsize - xhit;
	if (debug > 10)
	  cout << "TXSIZE/2*xsize=" << (TXSIZE/2)*xsize << endl;
      } else {
	dx = xrec - ((TXSIZE/2)-0.5)*xsize - xhit;
      }
      sx[qbin] += dx; sx2[qbin] += dx*dx; scx[qbin] += dx*dx/(sigmax*sigmax);
      if(probx > 1.e-2) {sxp[qbin] += dx; sxp2[qbin] += dx*dx; scxp[qbin] += dx*dx/(sigmax*sigmax);}

      if (debug > 1) {
	cout << "xrec=" << xrec << "  xhit=" << xhit << "     "
	     << "yrec=" << xrec << "  yhit=" << xhit << endl;
	cout << "Shifted positions: dx=" << dx << "   dy=" << dy << endl;
      }

      hp[10]->Fill(dx);
      hp[11+iqbin]->Fill(dx);
      hp[15]->Fill(log10probQ);
      hp[16+iqbin]->Fill(log10probQ);
      hp[24]->Fill((double)probxy);
      hp[25]->Fill((double)probQ);
      
      if(nypix == 1) {hp[22]->Fill(dy);} else {hp[26]->Fill(dy);}
      if(nxpix == 1) {hp[23]->Fill(dx);} else {hp[27]->Fill(dx);}

      // 
      // Fill the Fastsim histograms
      //
      //cout << "nxpix " << nxpix << endl;
      (void) fastSimResHistoStore.Fill( dx, dy, cotalpha, cotbeta, qbin, nxpix, nypix );

    }
  }  // end of the loop over hits i.e. end of ... while(fscanf(...))
  //
  // ========================================================================
  

  //  Determine current time (end of the job)
  gettimeofday(&now1, &timz);

  //  Timing printout
  deltas = now1.tv_sec - now0.tv_sec;
  deltaus = now1.tv_usec - now0.tv_usec;
  deltat = ((double)deltaus)/1000000.;
  deltat += (double)deltas;
  printf("ellapsed time = %f seconds \n", deltat);

  // Other printout
  printf(" total events = %d, probx > 10^{-3} = %d, proby > 10^{-3} = %d \n", totale, goodx, goody);
  printf(" low q failures = %d, malformed clusters = %d \n", nbin[4], nbad);
	   
  for(j=0; j<5; ++j) {
    sy[j] /= (float)nbin[j]; sy2[j] /= (float)nbin[j]; scy[j] /= (float)nbin[j];
    sy2[j] = sqrt((double)(sy2[j] - sy[j]*sy[j]));
    printf(" avg y residual[%1d] = %f +- %f, avg y chi^2 = %f \n", j, sy[j], sy2[j], scy[j]);       
    sx[j] /= (float)nbin[j]; sx2[j] /= (float)nbin[j]; scx[j] /= (float)nbin[j];
    sx2[j] = sqrt((double)(sx2[j] - sx[j]*sx[j]));
    printf(" avg x residual[%1d] = %f +- %f, avg x chi^2 = %f \n", j, sx[j], sx2[j], scx[j]); 
  }
  printf(" After 10^{-2}(x) 10^{-3}(y) probability cuts: \n");
  
  for(j=0; j<5; ++j) {
    syp[j] /= (float)nbin[j]; syp2[j] /= (float)nbin[j]; scyp[j] /= (float)nbin[j];
    syp2[j] = sqrt((double)(syp2[j] - syp[j]*syp[j]));
    printf(" avg y residual[%1d] = %f +- %f, avg y chi^2 = %f \n", j, syp[j], syp2[j], scyp[j]);       
    sxp[j] /= (float)nbin[j]; sxp2[j] /= (float)nbin[j]; scxp[j] /= (float)nbin[j];
    sxp2[j] = sqrt((double)(sxp2[j] - sxp[j]*sxp[j]));
    printf(" avg x residual[%1d] = %f +- %f, avg x chi^2 = %f \n", j, sxp[j], sxp2[j], scxp[j]); 
  }
  
  
  // 
  //  Histograms plotting
  //
  for(i=0; i<5; ++i) {hp[i]->Fit("gaus"); hp[i+10]->Fit("gaus");}
  for(i=22; i<24; ++i) {hp[i]->Fit("gaus");}
  for(i=26; i<28; ++i) {hp[i]->Fit("gaus");}
  
  //  Create an output filename for this run.  Need to run pstopdf to conver to pdf
  //  (was getting corrupt pdf files when saving to pdf).
  
  sprintf(outfile0,"pixel_histos%5.5d_%3.3d_%i.pdf[",nfile,ID,nbits);
  sprintf(outfile1,"pixel_histos%5.5d_%3.3d_%i.pdf",nfile,ID,nbits);
  sprintf(outfile2,"pixel_histos%5.5d_%3.3d_%i.pdf]",nfile,ID,nbits);
  TCanvas c1("c1", header);
  c1.SetFillStyle(4000);
  c1.Print(outfile0);
  for(i=0; i<30; ++i) {
    hp[i]->Draw();
    if(i != 29){ c1.Print(outfile1);}
    else{ c1.Print(outfile2);}
  }
  // for( int ii=0; ii<cotbetaBins; ii++)
  //   for( int jj=0; jj<cotalphaBins; jj++)   {
  //     qbinForward[ii][jj]->Draw();
  //     c1.Print(outfile1);
  //   }
  // c1.Print(outfile2);

  // Report maximal abs cotalpha and abs cotbeta seen in this file:
  cout << "max_cotalpha=" << max_cotalpha
       << "  max_cotbeta=" << max_cotbeta 
       << "min_cotalpha=" << min_cotalpha
       << "  min_cotbeta=" << min_cotbeta
       << endl;
  
  
  return 0;
} // MAIN__
// At this point fastSimResHistoStore object goes out of scope,
// and in the destructor the histograms are written out into the file.




// ***************************************************************** 
//! Calculate 21 gaussianly-distributed random numbers.
//! \param x - a vector holding 21 random numbers generated with rms = 1.            
// ***************************************************************** 
int triplg(std::vector<float>& x)
{
  // Initialized data 
  
  static int fcall = -1;
  
  // Local variables 
  static float r1, r2;
  static int i__;
  static float r__;
  static int ibase;
  static std::vector<float> rbuff(TYTEN);
  static float twopi;
  double arg, phi;
  
  // Function Body 

  //  Initalize the parameters 
  
  if (fcall) {
    twopi = 2.*acos((double)-1.);
    ibase = TYTEN;
    fcall = 0;
  }
  
  //  If all random numbers used up, generate 210 more 
  
  if (ibase == TYTEN) {
    for (i__ = 0; i__ < TYTEN-1; i__ += 2) {
      r1 = ((float)random())/((float)RAND_MAX);
      r2 = ((float)random())/((float)RAND_MAX);
      arg = (double)(1. - r1);
      if (arg < 1.e-30) {arg = 1.e-30;}
      r__ = sqrt(log(arg) * (-2.));
      phi = twopi * r2;
      rbuff[i__] = r__ * cos(phi);
      rbuff[i__+1] = r__ * sin(phi);
    }
    ibase = 0;
  }
  for (i__ = 0; i__ < TYSIZE; ++i__) {
    x[i__] = rbuff[ibase + i__];
  }
  ibase += TYSIZE;
  return 0;
} // triplg 

