// all users - please change the name of this file to lhcbStyle.C
// Commits to lhcbdocs svn of .C files are not allowed
//

//
#include <lhcbStyle.h>



void lhcb::lhcbStyle()
{
//	if (option == "LHCb") {
		// Use times new roman, precision 2
	  Int_t kLHCbFont        = 132;  // Old LHCb style: 62;
	  // Line thickness
	  Double_t lhcbWidth    = 2.00; // Old LHCb style: 3.00;
	  // Text size
	  Double_t lhcbTSize    = 0.06;

	  // use plain black on white colors
	  gROOT->SetStyle("Plain");
	  TStyle *lhcbStyle= new TStyle("lhcbStyle","LHCb plots style");

	  //lhcbStyle->SetErrorX(0); //  don't suppress the error bar along X

	  lhcbStyle->SetFillColor(1);
	  lhcbStyle->SetFillStyle(1001);   // solid
	  lhcbStyle->SetFrameFillColor(0);
	  lhcbStyle->SetFrameBorderMode(0);
	  lhcbStyle->SetPadBorderMode(0);
	  lhcbStyle->SetPadColor(0);
	  lhcbStyle->SetCanvasBorderMode(0);
	  lhcbStyle->SetCanvasColor(0);
	  lhcbStyle->SetStatColor(0);
	  lhcbStyle->SetLegendBorderSize(0);

	  // If you want the usual gradient palette (blue -> red)
	  lhcbStyle->SetPalette(1);
	  // If you want colors that correspond to gray scale in black and white:
	  int colors[8] = {0,5,7,3,6,2,4,1};
//	  lhcbStyle->SetPalette(8,colors);
	  lhcbStyle->SetPalette(51);

	  // set the paper & margin sizes
	  lhcbStyle->SetPaperSize(20,26);
	  lhcbStyle->SetPadTopMargin(0.09);
//	  lhcbStyle->SetPadRightMargin(0.05); // increase for colz plots
    lhcbStyle->SetPadRightMargin(0.15); // increase for colz plots
	  lhcbStyle->SetPadBottomMargin(0.16);
	  lhcbStyle->SetPadLeftMargin(0.14);

	  // use large fonts
	  lhcbStyle->SetTextFont(kLHCbFont);
	  lhcbStyle->SetTextSize(lhcbTSize);
	  lhcbStyle->SetLabelFont(kLHCbFont,"x");
	  lhcbStyle->SetLabelFont(kLHCbFont,"y");
	  lhcbStyle->SetLabelFont(kLHCbFont,"z");
	  lhcbStyle->SetLabelSize(lhcbTSize,"x");
	  lhcbStyle->SetLabelSize(lhcbTSize,"y");
	  lhcbStyle->SetLabelSize(lhcbTSize,"z");
	  lhcbStyle->SetTitleFont(kLHCbFont);
	  lhcbStyle->SetTitleFont(kLHCbFont,"x");
	  lhcbStyle->SetTitleFont(kLHCbFont,"y");
	  lhcbStyle->SetTitleFont(kLHCbFont,"z");
	  lhcbStyle->SetTitleSize(1.2*lhcbTSize,"x");
	  lhcbStyle->SetTitleSize(1.2*lhcbTSize,"y");
	  lhcbStyle->SetTitleSize(1.2*lhcbTSize,"z");

	  // use medium bold lines and thick markers
	  lhcbStyle->SetLineWidth(lhcbWidth);
	  lhcbStyle->SetFrameLineWidth(lhcbWidth);
	  lhcbStyle->SetHistLineWidth(lhcbWidth);
	  lhcbStyle->SetFuncWidth(lhcbWidth);
	  lhcbStyle->SetGridWidth(lhcbWidth);
	  lhcbStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
	  lhcbStyle->SetMarkerStyle(20);
	  lhcbStyle->SetMarkerSize(1.0);

	  // label offsets
	  lhcbStyle->SetLabelOffset(0.010,"X");
	  lhcbStyle->SetLabelOffset(0.010,"Y");

	  // by default, do not display histogram decorations:
	  lhcbStyle->SetOptStat(0);
	  //lhcbStyle->SetOptStat("emr");  // show only nent -e , mean - m , rms -r
	  // full opts at http://root.cern.ch/root/html/TStyle.html#TStyle:SetOptStat
	  lhcbStyle->SetStatFormat("6.3g"); // specified as c printf options

	  lhcbStyle->SetOptTitle(0);


	  lhcbStyle->SetOptFit(0);
	  //lhcbStyle->SetOptFit(1011); // order is probability, Chi2, errors, parameters
	  //titles
	  lhcbStyle->SetTitleOffset(0.95,"X");
	  lhcbStyle->SetTitleOffset(0.95,"Y");
	  lhcbStyle->SetTitleOffset(1.2,"Z");
	  lhcbStyle->SetTitleFillColor(0);
	  lhcbStyle->SetTitleStyle(0);
	  lhcbStyle->SetTitleBorderSize(0);
	  lhcbStyle->SetTitleFont(kLHCbFont,"title");
	  lhcbStyle->SetTitleX(0.0);
	  lhcbStyle->SetTitleY(1.0);
	  lhcbStyle->SetTitleW(1.0);
	  lhcbStyle->SetTitleH(0.05);

	  // look of the statistics box:
	  lhcbStyle->SetStatBorderSize(0);
	  lhcbStyle->SetStatFont(kLHCbFont);
	  lhcbStyle->SetStatFontSize(0.05);
	  lhcbStyle->SetStatX(0.9);
	  lhcbStyle->SetStatY(0.9);
	  lhcbStyle->SetStatW(0.25);
	  lhcbStyle->SetStatH(0.15);

	  // put tick marks on top and RHS of plots
	  lhcbStyle->SetPadTickX(1);
	  lhcbStyle->SetPadTickY(1);

	  // histogram divisions: only 5 in x to avoid label overlaps
	  lhcbStyle->SetNdivisions(505,"x");
	  lhcbStyle->SetNdivisions(510,"y");

	  gROOT->SetStyle("lhcbStyle");
	  gROOT->ForceStyle();
//	}

//	else {
//  	gROOT->SetStyle("Plain");
//
//		// text font
//		int font = 132;
//		gStyle->SetTitleFont(font, "xyz");	// set the all 3 axes title font
//		gStyle->SetTitleFont(font, " ");	// set the pad title font
//		gStyle->SetLabelFont(font, "xyz");
//		gStyle->SetTextFont(font);
//		gStyle->SetStatFont(font);
//
//		// no default boxes
//		gStyle->SetOptStat(0);
//		gStyle->SetOptTitle(0);
//		gStyle->SetOptFit(0);
//
//		// rainbow colors in 2d plots
//		gStyle->SetPalette(1);
//
//		// use plain black on white colors
//		gStyle->SetFrameBorderMode(0);
//		gStyle->SetFrameFillColor(0);
//		gStyle->SetCanvasBorderMode(0);
//		gStyle->SetPadBorderMode(0);
//		gStyle->SetPadColor(0);
//		gStyle->SetCanvasColor(0);
//		gStyle->SetStatColor(0);
//		gStyle->SetTitleFillColor(0);
//
//  	/// \todo Put some more documentation in here.
//
//		// canvas default size
//		gStyle->SetCanvasDefH(480);
//		gStyle->SetCanvasDefW(640);
//
//		// text
//		gStyle->SetTextSize(0.07);
//		gStyle->SetTextAlign(13);
//
//		// set the paper & margin sizes Changed the values to look a bit nicer TOBI 20120823 does not work properly yet
//		//gStyle->SetPadTopMargin(0.08);
//		//gStyle->SetPadRightMargin(0.10);
//		//gStyle->SetPadLeftMargin(0.17);
//		//gStyle->SetPadBottomMargin(0.15);
//		//if ( option=="2d" ) gStyle->SetPadBottomMargin(0.15);
//		//gStyle->SetPadLeftMargin(0.14);
//		//if ( option=="2d" ) gStyle->SetPadRightMargin(0.18);
//
//		//old margins
//		gStyle->SetPadTopMargin(0.05);
//		gStyle->SetPadRightMargin(0.05);
//		gStyle->SetPadBottomMargin(0.14);
//		if ( option=="2d" ) gStyle->SetPadBottomMargin(0.15);
//		gStyle->SetPadLeftMargin(0.14);
//		if ( option=="2d" ) gStyle->SetPadRightMargin(0.18);
//
//		// use bold lines and markers
//		gStyle->SetMarkerStyle(8);
//		gStyle->SetHistLineWidth(2);
//		gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
//
//		// get rid of X error bars
//		gStyle->SetErrorX(0.001);
//
//		// axes
//		gStyle->SetTitleSize(0.065,	"xyz");
//		gStyle->SetLabelSize(0.06,	"xyz");
//		gStyle->SetTitleOffset(1.07,	"y");
//		gStyle->SetLabelOffset(0.015,	"x");
//		gStyle->SetNdivisions(507,	"x");
//
//		// put tick marks on top and right hand side of plots
//		//if ( option!="2d" )
//		//{
//    gStyle->SetPadTickX(1);
//    gStyle->SetPadTickY(1);
//		//}
//	}

}

