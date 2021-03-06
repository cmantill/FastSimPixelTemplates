//
//  The implementation of the PixelResolutinHistograms.cc class.  Please
//  look at PixelResolutionHistograms.h header file for the interface.
//
//------------------------------------------------------------------------------

// The switch, undefined in CMSSW release, and defined by standalone compilation script:


#ifdef SI_PIXEL_TEMPLATE_STANDALONE
//
//--- Stand-alone: Include a the header file from the local directory, as well as
//    a dummy implementation of SimpleHistogramGenerator.
//
#include "PixelResolutionHistograms.h"
//
class TH1F;
class TH2F;
class TF1;
class SimpleHistogramGenerator {
public:
  SimpleHistogramGenerator(TH1F * hist) : hist_(hist) {};
private:
  TH1F * hist_;
};
//
#else
//--- We're inside a CMSSW release: Include the real thing.
//
#include "FastSimulation/TrackingRecHitProducer/interface/PixelResolutionHistograms.h"
#include "FastSimulation/Utilities/interface/RandomEngineAndDistribution.h"
#include "FastSimulation/Utilities/interface/SimpleHistogramGenerator.h"
//
#endif




// Generic C stuff
#include <math.h>       /* floor */
#include <iostream>
#include <string>

// ROOT
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>

// Global definitions 
const float cmtomicron = 10000.0;


//------------------------------------------------------------------------------
//  Constructor: Books the FastSim Histograms, given the input parameters
//  which are provided as arguments. These variables are then const inside
//  the class. (That is, once we make the histograms, we can't change the
//  definition of the binning.)
//------------------------------------------------------------------------------
PixelResolutionHistograms::
PixelResolutionHistograms( const char * filename,   // ROOT file for histograms
			   const char * rootdir,    // Subdirectory in the file, "" if none
			   const char * descTitle,  // Descriptive title	     
			   unsigned int detType,    // Where we are... (&&& do we need this?)
			   double 	cotbetaBinWidth,
			   double 	cotbetaLowEdge,
			   int	cotbetaBins,
			   double	cotalphaBinWidth,
			   double	cotalphaLowEdge,
			   int	cotalphaBins) //,
			   //int	qbinWidth,
			   //int	qbins )
  : weOwnHistograms_(true),                          // we'll be making some histos
    detType_          ( detType ),
    cotbetaBinWidth_  ( cotbetaBinWidth ),
    cotbetaLowEdge_   ( cotbetaLowEdge ), 
    cotbetaBins_      ( cotbetaBins ),
    cotalphaBinWidth_ ( cotalphaBinWidth  ),
    cotalphaLowEdge_  ( cotalphaLowEdge ),
    cotalphaBins_     ( cotalphaBins ),
    qbinWidth_        ( 1 ),
    qbins_            ( 4 ),
    binningHisto_(nullptr),
    resMultiPixelXHist_(), resSinglePixelXHist_(),  // all to nullptr
    resMultiPixelYHist_(), resSinglePixelYHist_(),  // all to nullptr
    qbinHist_(),                                    // all to nullptr
    file_(nullptr),
    status_(0),
    resMultiPixelXGen_(), resSinglePixelXGen_(), 
    resMultiPixelYGen_(), resSinglePixelYGen_(),
    qbinGen_()
{

  file_ = new TFile( filename, "recreate" );
  //Resolution binning
  // const double 	cotbetaBinWidth = 1.0;
  // const double 	cotbetaLowEdge	= -11.5 ;
  // const int	cotbetaBins	= 23;
  // const double	cotalphaBinWidth	= 0.08 ;
  // const double	cotalphaLowEdge	= -0.36 ;
  // const int	cotalphaBins	= 9;
  // const int	qbinWidth	= 1;
  // const int	qbins		= 4;

  // Dummy 2D histogram to store binning:
  binningHisto_ = new TH2F( "ResHistoBinning", descTitle,
			    cotbetaBins,  cotbetaLowEdge,  cotbetaLowEdge  + cotbetaBins*cotbetaBinWidth,
			    cotalphaBins, cotalphaLowEdge, cotalphaLowEdge + cotalphaBins*cotalphaBinWidth );

  // Store detType in the underflow bin
  binningHisto_->SetBinContent(0, 0, detType_);
  
  // All other histograms:
  Char_t histo[200];
  Char_t title[200];
  //
  for( int ii=0; ii<cotbetaBins_; ii++ ) {
    for( int jj=0; jj<cotalphaBins_; jj++ ) {
      for( int kk=0; kk<qbins_; kk++ ) {
	//
	sprintf( histo, "hx%d1%02d%d%d", detType_, ii+1, jj+1, kk+1 );  //information of bits of histogram names
	//--- First bit 1/0 barrel/forward, second 1/0 multi/single, cotbeta, cotalpha, qbins
	sprintf( title, "cotbeta %.1f-%.1f cotalpha %.2f-%.2f qbin %d npixel>1 X",
		 cotbetaLowEdge_ + ii*cotbetaBinWidth_ , cotbetaLowEdge_ + (ii+1)*cotbetaBinWidth_,
		 cotalphaLowEdge_ +jj*cotalphaBinWidth_, cotalphaLowEdge_ +(jj+1)*cotalphaBinWidth_,
		 kk+1 );
	//
	resMultiPixelXHist_ [ ii ][ jj ][ kk ] = new TH1F(histo, title, 1000, -0.05, 0.05);

	sprintf( histo, "hy%d1%02d%d%d", detType_, ii+1, jj+1, kk+1 );
	sprintf( title, "cotbeta %.1f-%.1f cotalpha %.2f-%.2f qbin %d npixel>1 Y",
		 cotbetaLowEdge_ + ii*cotbetaBinWidth_ , cotbetaLowEdge_ + (ii+1)*cotbetaBinWidth_,
		 cotalphaLowEdge_ +jj*cotalphaBinWidth_, cotalphaLowEdge_ +(jj+1)*cotalphaBinWidth_,
		 kk+1 );
	//
	resMultiPixelYHist_ [ ii ][ jj ][ kk ] = new TH1F(histo, title, 1000, -0.05, 0.05);
      }
    }
  }


  for( int ii=0; ii<cotbetaBins_; ii++) {
    for( int jj=0; jj<cotalphaBins_; jj++) {

      sprintf( histo, "hx%d0%02d%d", detType_, ii+1, jj+1 );  //information of bits of histogram names
      //first bit 1/0 barrel/forward, second 1/0 multi/single, cotbeta, cotalpha
      sprintf( title, "cotbeta %.1f-%.1f cotalpha %.2f-%.2f npixel=1 X",
		 cotbetaLowEdge_ + ii*cotbetaBinWidth_ , cotbetaLowEdge_ + (ii+1)*cotbetaBinWidth_,
		 cotalphaLowEdge_ +jj*cotalphaBinWidth_, cotalphaLowEdge_ +(jj+1)*cotalphaBinWidth_ );
      //
      resSinglePixelXHist_ [ ii ][ jj ] = new TH1F(histo, title, 1000, -0.05, 0.05);

      sprintf( histo, "hy%d0%02d%d", detType_, ii+1, jj+1 );
      sprintf( title, "cotbeta %.1f-%.1f cotalpha %.2f-%.2f npixel=1 Y",
		 cotbetaLowEdge_ + ii*cotbetaBinWidth_ , cotbetaLowEdge_ + (ii+1)*cotbetaBinWidth_,
	       cotalphaLowEdge_ +jj*cotalphaBinWidth_, cotalphaLowEdge_ +(jj+1)*cotalphaBinWidth_ );
      //
      resSinglePixelYHist_ [ ii ][ jj ] = new TH1F(histo, title, 1000, -0.05, 0.05);

      sprintf( histo, "hqbin%d%02d%d", detType_, ii+1, jj+1 );
      sprintf( title, "cotbeta %.1f-%.1f cotalpha %.2f-%.2f qbin",
	       cotbetaLowEdge_ + ii*cotbetaBinWidth_ , cotbetaLowEdge_ + (ii+1)*cotbetaBinWidth_,
	       cotalphaLowEdge_ +jj*cotalphaBinWidth_, cotalphaLowEdge_ +(jj+1)*cotalphaBinWidth_ );
      //
      qbinHist_ [ ii ][ jj ] = new TH1F(histo, title, 4, -0.49, 3.51);
      
    }
  }

}



//------------------------------------------------------------------------------
//  Another constructor: load the histograms from one file.
//     filename = full path to filename
//     rootdir  = ROOT directory inside the file
//
//  The other parameters are the same (needed later) and must correspond
//  to the histograms we are loading from the file.
//------------------------------------------------------------------------------
PixelResolutionHistograms::
PixelResolutionHistograms( const char * filename, 
			   const char * rootdir )
  : weOwnHistograms_(false),                          // we'll be making some histos
    detType_          (-1),
    cotbetaBinWidth_  (0),
    cotbetaLowEdge_   (0), 
    cotbetaBins_      (0),
    cotalphaBinWidth_ (0),
    cotalphaLowEdge_  (0),
    cotalphaBins_     (0),
    qbinWidth_        (1),       // &&& Is this ever going to change?
    qbins_            (4),       // &&& Is this ever going to change?
    binningHisto_(nullptr),
    resMultiPixelXHist_(), resSinglePixelXHist_(),  // all to nullptr
    resMultiPixelYHist_(), resSinglePixelYHist_(),  // all to nullptr
    qbinHist_(),                                    // all to nullptr
    file_(nullptr),
    status_(0),
    resMultiPixelXGen_(), resSinglePixelXGen_(), 
    resMultiPixelYGen_(), resSinglePixelYGen_(),
    qbinGen_()
{
  Char_t histo[200];       // the name of the histogram
  Char_t title[200];       // histo title, for debugging and sanity checking (compare inside file)
  TH1F * tmphist = nullptr;  // cache for histo pointer

  //--- Open the file for reading.
  file_ = new TFile( filename  ,"READ");
  if ( !file_ ) {
    status_ = 1;
    return;
  }

  //--- The dummy 2D histogram with the binning of cot\beta and cot\alpha:
  binningHisto_ = (TH2F*) file_->Get( Form( "%s%s" , rootdir, "ResHistoBinning" ) );
  if ( !binningHisto_ ) {
    status_ = 11;
    return;
  }

  //--- Fish out detType from the underflow bin:
  detType_ = binningHisto_->GetBinContent(0, 0);

  //--- Now we fill the binning variables:
  cotbetaAxis_      = binningHisto_->GetXaxis();
  cotbetaBinWidth_  = binningHisto_->GetXaxis()->GetBinWidth(1);  // assume all same width
  cotbetaLowEdge_   = binningHisto_->GetXaxis()->GetXmin();       // low edge of the first bin
  cotbetaBins_      = binningHisto_->GetXaxis()->GetNbins();
  cotalphaAxis_     = binningHisto_->GetYaxis();
  cotalphaBinWidth_ = binningHisto_->GetYaxis()->GetBinWidth(1);  // assume all same width;
  cotalphaLowEdge_  = binningHisto_->GetYaxis()->GetXmin();       // low edge of the first bin;
  cotalphaBins_     = binningHisto_->GetYaxis()->GetNbins();
  
  
  for( int ii=0; ii<cotbetaBins_; ii++ ) {
    for( int jj=0; jj<cotalphaBins_; jj++ ) {
      for( int kk=0; kk<qbins_; kk++ ) {
	//
	sprintf( histo, "hx%d1%02d%d%d", detType_, ii+1, jj+1, kk+1 );  //information of bits of histogram names
	//--- First bit 1/0 barrel/forward, second 1/0 multi/single, cotbeta, cotalpha, qbins
	sprintf( title, "cotbeta %.1f-%.1f cotalpha %.2f-%.2f qbin %d npixel>1 X",
		 cotbetaLowEdge_ + ii*cotbetaBinWidth_ , cotbetaLowEdge_ + (ii+1)*cotbetaBinWidth_,
		 cotalphaLowEdge_ +jj*cotalphaBinWidth_, cotalphaLowEdge_ +(jj+1)*cotalphaBinWidth_,
		 kk+1 );
	std::cout << "Debug: histo=" << std::string(histo) << " title=" << std::string(title) 
		  << std::endl;
	//
	tmphist = (TH1F*) file_->Get( Form( "%s%s" , rootdir, histo ) );
	if ( !tmphist ) {
	  status_ = 2;
	  return;
	}
	std::cout << "Debug: found histo with title = " << std::string( tmphist->GetTitle() ) 
		  << std::endl;
	resMultiPixelXHist_ [ ii ][ jj ][ kk ] = tmphist;
	resMultiPixelXGen_  [ ii ][ jj ][ kk ] = new SimpleHistogramGenerator( tmphist );


	sprintf( histo, "hy%d1%02d%d%d", detType_, ii+1, jj+1, kk+1 );
	sprintf( title, "cotbeta %.1f-%.1f cotalpha %.2f-%.2f qbin %d npixel>1 Y",
		 cotbetaLowEdge_ + ii*cotbetaBinWidth_ , cotbetaLowEdge_ + (ii+1)*cotbetaBinWidth_,
		 cotalphaLowEdge_ +jj*cotalphaBinWidth_, cotalphaLowEdge_ +(jj+1)*cotalphaBinWidth_,
		 kk+1 );
	std::cout << "Debug: histo=" << std::string(histo) << " title=" << std::string(title)  
		  << std::endl;
	//
	tmphist = (TH1F*) file_->Get( Form( "%s%s" , rootdir, histo ) );
	if ( !tmphist ) {
	  status_ = 3;
	  return;
	}
	std::cout << "Debug: found histo with title = " << std::string( tmphist->GetTitle() ) 
		  << std::endl;
	resMultiPixelYHist_ [ ii ][ jj ][ kk ] = tmphist;
	resMultiPixelYGen_  [ ii ][ jj ][ kk ] = new SimpleHistogramGenerator( tmphist );
      }
    }
  }


  for( int ii=0; ii<cotbetaBins_; ii++) {
    for( int jj=0; jj<cotalphaBins_; jj++) {

      sprintf( histo, "hx%d0%02d%d", detType_, ii+1, jj+1 );  //information of bits of histogram names
      //--- First bit 1/0 barrel/forward, second 1/0 multi/single, cotbeta, cotalpha
      sprintf( title, "cotbeta %.1f-%.1f cotalpha %.2f-%.2f npixel=1 X",
	       cotbetaLowEdge_ + ii*cotbetaBinWidth_ , cotbetaLowEdge_ + (ii+1)*cotbetaBinWidth_,
	       cotalphaLowEdge_ +jj*cotalphaBinWidth_, cotalphaLowEdge_ +(jj+1)*cotalphaBinWidth_ );
      std::cout << "Debug: histo=" << std::string(histo) << " title=" << std::string(title)  
		<< std::endl;
      //
      tmphist = (TH1F*) file_->Get( Form( "%s%s" , rootdir, histo ) );
      if ( !tmphist ) {
	  status_ = 4;
	  return;
      }
      std::cout << "Debug: found histo with title = " << std::string( tmphist->GetTitle() ) 
		<< std::endl;
      resSinglePixelXHist_ [ ii ][ jj ] = tmphist;
      resSinglePixelXGen_  [ ii ][ jj ] = new SimpleHistogramGenerator( tmphist );


      sprintf( histo, "hy%d0%02d%d", detType_, ii+1, jj+1 );
      sprintf( title, "cotbeta %.1f-%.1f cotalpha %.2f-%.2f npixel=1 Y",
	       cotbetaLowEdge_ + ii*cotbetaBinWidth_ , cotbetaLowEdge_ + (ii+1)*cotbetaBinWidth_,
	       cotalphaLowEdge_ +jj*cotalphaBinWidth_, cotalphaLowEdge_ +(jj+1)*cotalphaBinWidth_ );
      std::cout << "Debug: histo=" << std::string(histo) << " title=" << std::string(title)  
		<< std::endl;
      //
      tmphist = (TH1F*) file_->Get( Form( "%s%s" , rootdir, histo ) );
      if ( !tmphist ) {
	  status_ = 5;
	  return;
      }
      std::cout << "Debug: found histo with title = " << std::string( tmphist->GetTitle() ) 
		<< std::endl;
      resSinglePixelYHist_ [ ii ][ jj ] = tmphist;
      resSinglePixelYGen_  [ ii ][ jj ] = new SimpleHistogramGenerator( tmphist );

      sprintf( histo, "hqbin%d%02d%d", detType_, ii+1, jj+1 );
      sprintf( title, "cotbeta %.1f-%.1f cotalpha %.2f-%.2f qbin",
	       cotbetaLowEdge_ + ii*cotbetaBinWidth_ , cotbetaLowEdge_ + (ii+1)*cotbetaBinWidth_,
	       cotalphaLowEdge_ +jj*cotalphaBinWidth_, cotalphaLowEdge_ +(jj+1)*cotalphaBinWidth_ );
      std::cout << "Debug: histo=" << std::string(histo) << " title=" << std::string(title)  
		<< std::endl;
      //
      tmphist = (TH1F*) file_->Get( Form( "%s%s" , rootdir, histo ) );
      if ( !tmphist ) {
	  status_ = 6;
	  return;
      }
      std::cout << "Debug: found histo with title = " << std::string( tmphist->GetTitle() ) 
		<< std::endl;
      qbinHist_ [ ii ][ jj ] = tmphist;
      qbinGen_  [ ii ][ jj ] = new SimpleHistogramGenerator( tmphist );
    }
  }
}




//------------------------------------------------------------------------------
//  Destructor.  Use file_ pointer to tell whether we loaded the histograms
//  from a file (and do not own them), or we built them ourselves and thus need
//  to delete them.
//------------------------------------------------------------------------------
PixelResolutionHistograms::~PixelResolutionHistograms()
{
  //--- Delete histograms, but only if we own them. If 
  //--- they came from a file, let them be.
  //
  if ( ! weOwnHistograms_ ) {
    //--- Read the histograms from the TFile, the file will take care of them.
    file_->Close();
    delete file_ ;
    file_ = 0;
  }
  else {
    //--- Make sure that histograms are filled
    TF1 *mygaus = new TF1("mygaus","TMath::Gaus(x,0,.0005)",-0.005,0.005); 
    for( int ii=0; ii<cotbetaBins_; ii++ ) {                                                                                                                                      
      for( int jj=0; jj<cotalphaBins_; jj++ ) {                                                                                                                                   
	for( int kk=0; kk<qbins_; kk++ ) {                                                                                                                                      
	  if(resMultiPixelXHist_ [ ii ][ jj ][ kk ]->GetEntries() <= 5) resMultiPixelXHist_[ ii ][ jj ][ kk ]->FillRandom("mygaus",1000);
          if(resMultiPixelYHist_ [ ii ][ jj ][ kk ]->GetEntries() <= 5) resMultiPixelYHist_[ ii ][ jj ][ kk ]->FillRandom("mygaus",1000);
	}
      }
    }
    for( int ii=0; ii<cotbetaBins_; ii++ ) {                                                                                                                                      
      for( int jj=0; jj<cotalphaBins_; jj++ ) {                                                                                                                                
	if(resSinglePixelXHist_ [ ii ][ jj ]->GetEntries() <= 5) resSinglePixelXHist_[ ii ][ jj ]->FillRandom("mygaus",1000);
	if(resSinglePixelYHist_ [ ii ][ jj ]->GetEntries() <= 5) resSinglePixelYHist_[ ii ][ jj ]->FillRandom("mygaus",1000);
      }
    }

    //--- We made the histograms, so first write them inthe output ROOT file and close it.
    std::cout
      << "PixelResHistoStore: Writing the histograms to the output file. " // << std::string( filename )
      << std::endl;
    file_->Write();
    file_->Close();
    delete file_ ;
    file_ = 0;

    // Well, apparently, ROOT file has the ownership, and once the file is closed,
    // all of these guys are deleted.  So, we don't need to do anything.
    // //--- Then, delete them all.
    // delete binningHisto_ ;

    // for( int ii=0; ii<cotbetaBins_; ii++ ) {
    //   for( int jj=0; jj<cotalphaBins_; jj++ ) {
    // 	for( int kk=0; kk<qbins_; kk++ ) {
    // 	  delete resMultiPixelXHist_ [ ii ][ jj ][ kk ];
    // 	  delete resMultiPixelYHist_ [ ii ][ jj ][ kk ];
    // 	}
    //   }
    // }
    // for( int ii=0; ii<cotbetaBins_; ii++ ) {
    //   for( int jj=0; jj<cotalphaBins_; jj++ ) {
    // 	delete resSinglePixelXHist_ [ ii ][ jj ];
    // 	delete resSinglePixelYHist_ [ ii ][ jj ]; 
    // 	delete qbinHist_ [ ii ][ jj ];
    //   }
    // }
  } // else



  //--- Delete FastSim generators. (It's safe to delete a nullptr.)
  for( int ii=0; ii<cotbetaBins_; ii++ ) {
    for( int jj=0; jj<cotalphaBins_; jj++ ) {
      for( int kk=0; kk<qbins_; kk++ ) {
	delete resMultiPixelXGen_ [ ii ][ jj ][ kk ];
	delete resMultiPixelYGen_ [ ii ][ jj ][ kk ];
      }
    }
  }
  for( int ii=0; ii<cotbetaBins_; ii++ ) {
    for( int jj=0; jj<cotalphaBins_; jj++ ) {
      delete resSinglePixelXGen_ [ ii ][ jj ];
      delete resSinglePixelYGen_ [ ii ][ jj ]; 
      delete qbinGen_ [ ii ][ jj ];
    }
  }

}


  



//------------------------------------------------------------------------------
//  Fills the appropriate FastSim histograms.
//  Returns 0 if the relevant histogram(s) were found and filled, 1 if not.
//------------------------------------------------------------------------------
int
PixelResolutionHistograms::
Fill( double dx, double dy, double cotalpha, double cotbeta, 
      int qbin, int nxpix, int nypix ) 
{
  int icotalpha, icotbeta, iqbin ;
  icotalpha = (int)floor( (cotalpha - cotalphaLowEdge_)/cotalphaBinWidth_ ) ;
  icotbeta  = (int)floor( (cotbeta - cotbetaLowEdge_) /cotbetaBinWidth_ ) ;
  iqbin = qbin > 2 ? 3 : qbin;
  ///std::cout << icotalpha << " < " << cotalphaBins_ << " " << icotbeta << " < " << cotbetaBins_ << std::endl;
  if( icotalpha >= 0 && icotalpha < cotalphaBins_ && icotbeta >= 0 && icotbeta < cotbetaBins_ ) {
    //std::cout << "accepted " << std::endl;
    qbinHist_[icotbeta][icotalpha]->Fill((double)iqbin);
    if( nxpix == 1 )
      resSinglePixelXHist_[icotbeta][icotalpha]->Fill(dx/cmtomicron);
    else
      resMultiPixelXHist_[icotbeta][icotalpha][iqbin]->Fill(dx/cmtomicron);
    if( nypix == 1 )
      resSinglePixelYHist_[icotbeta][icotalpha]->Fill(dy/cmtomicron);
    else
      resMultiPixelYHist_[icotbeta][icotalpha][iqbin]->Fill(dy/cmtomicron);
  }
  //else { 
    //if( nxpix == 1 || nypix == 1) std::cout << "discarded nxpix " << nxpix << " nypix " << nypix << " cotalpha " << cotalpha << " "<<  icotalpha << " < " << cotalphaBins_ << " cotbeta " << cotbeta << " "  << icotbeta << " < " << cotbetaBins_ << std::endl;
  //}
  
  return 0;  
}



  //--- Get generators, for resolution in X and Y.  Use in FastSim.
const SimpleHistogramGenerator * 
PixelResolutionHistograms::
getGeneratorX( double cotalpha, double cotbeta, int qbin, bool single )
{
  int icotalpha, icotbeta, iqbin ;
  icotalpha = (int)floor( (cotalpha - cotalphaLowEdge_)/cotalphaBinWidth_ ) ;
  icotbeta  = (int)floor( (cotbeta - cotbetaLowEdge_) /cotbetaBinWidth_ ) ;
  iqbin = qbin > 2 ? 3 : qbin;
  std::cout << "Debug (getGeneratorX): "<<icotalpha<<" "<<icotbeta<<" "<<iqbin<< std::endl;
  //
  if( icotalpha >= 0 && icotalpha < cotalphaBins_ && icotbeta >= 0 && icotbeta < cotbetaBins_ ) {
    if( single )
      return resSinglePixelXGen_[icotbeta][icotalpha];
    else
      return resMultiPixelXGen_[icotbeta][icotalpha][iqbin];
  }
  else
    return nullptr;
}

const SimpleHistogramGenerator * 
PixelResolutionHistograms::
getGeneratorY( double cotalpha, double cotbeta, int qbin, bool single )
{
  int icotalpha, icotbeta, iqbin ;
  icotalpha = (int)floor( (cotalpha - cotalphaLowEdge_)/cotalphaBinWidth_ ) ;
  icotbeta  = (int)floor( (cotbeta - cotbetaLowEdge_) /cotbetaBinWidth_ ) ;
  iqbin = qbin > 2 ? 3 : qbin;
  std::cout << "Debug (getGeneratorY): "<<icotalpha<<" "<<icotbeta<<" "<<iqbin<< std::endl;
  //
  if( icotalpha >= 0 && icotalpha < cotalphaBins_ && icotbeta >= 0 && icotbeta < cotbetaBins_ ) {
    if( single )
      return resSinglePixelYGen_[icotbeta][icotalpha];
    else
      return resMultiPixelYGen_[icotbeta][icotalpha][iqbin];
  }
  else
    return nullptr;
}
