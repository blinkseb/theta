#include "interface/plugin.hpp"
#include "interface/histogram-function.hpp"
#include "root/root_histogram.hpp"

#include "TH1.h"
#include "TFile.h"

using namespace theta;
using namespace theta::plugin;
using namespace std;

root_histogram::root_histogram(const Configuration & ctx){
    string filename = ctx.replace_theta_dir(ctx.setting["filename"]);
    string histoname = ctx.setting["histoname"];
    int rebin = 1;
    double range_low = -999;
    double range_high = -999;
    if(ctx.setting.exists("rebin")){
       rebin = ctx.setting["rebin"];
    }
    if(ctx.setting.exists("range")){
       range_low = ctx.setting["range"][0];
       range_high = ctx.setting["range"][1];
    }
    bool use_errors = false;
    if(ctx.setting.exists("use_errors")){
         use_errors = ctx.setting["use_errors"];
    }
    
    TFile file(filename.c_str(), "read");
    if(file.IsZombie()){
        stringstream s;
        s << "Could not open file '" << filename << "'";
       throw ConfigurationException(s.str());
    }
    TH1* histo = dynamic_cast<TH1*>(file.Get(histoname.c_str()));
    if(not histo){
       stringstream s;
       s << "Did not find TH1 '" << histoname << "' in file '" << filename << "'";
       throw ConfigurationException(s.str());
    }
    if(ctx.setting.exists("normalize_to")){
       double norm = HistogramFunctionUtils::read_normalize_to(ctx.setting);
       double histo_integral = histo->Integral();
       if(histo_integral==0){
           if(norm!=0){
               throw ConfigurationException("specified non-zero 'normalize_to' setting but original Histogram's integral is zero");
           }
       }
       else{
           histo->Scale(norm / histo_integral);
       }
    }

    Histogram h, h_error;
    //take care of 1D histograms:
    if(histo->GetDimension() == 1){
       histo->Rebin(rebin);
       int bin_low = 1;
       int bin_high = histo->GetNbinsX();
       if(range_low!=-999){
          bin_low = histo->GetXaxis()->FindBin(range_low);
          if(bin_low > 0 && range_low!=histo->GetXaxis()->GetBinLowEdge(bin_low)){
              throw ConfigurationException("'range' setting incompatible with bin borders.");
          }
       }
       if(range_high!=-999){
          if(range_high > histo->GetXaxis()->GetXmax()) bin_high = histo->GetNbinsX()+1;
          else{
             bin_high = histo->GetXaxis()->FindBin(range_high);
             --bin_high;
             if(range_high!=histo->GetXaxis()->GetBinUpEdge(bin_high)){
                 throw ConfigurationException("'range' setting incompatible with bin borders.");
             }
          }
       }
       int nbins = bin_high - bin_low + 1;
    
       double xmin = histo->GetXaxis()->GetXmin();
       double xmax = histo->GetXaxis()->GetXmax();
       if(bin_low > 0)
          xmin = histo->GetXaxis()->GetBinLowEdge(bin_low);
       if(bin_high <= nbins)
          xmax = histo->GetXaxis()->GetBinUpEdge(bin_high);
    
    
       h.reset(nbins, xmin, xmax);
       h_error.reset(nbins, xmin, xmax);
       for(int i = bin_low; i <= bin_high; i++){
          double content = histo->GetBinContent(i);
          h.set(i - bin_low + 1, content);
          //h_error contains the relative errors:
          if(use_errors && content > 0.0){
             h_error.set(i - bin_low + 1, histo->GetBinError(i) / content);
          }
       }
    }
    //take care of 2D/3D histograms:
    else{
        if(rebin!=1 || range_low!=-999 || range_high!=-999) throw ConfigurationException("rebin and range setting not allowed for multidimensional Histograms!");
       size_t nbins_x = histo->GetNbinsX();
       size_t nbins_y = histo->GetNbinsY();
       size_t nbins_z = histo->GetNbinsZ();
       h.reset(nbins_x * nbins_y * max<size_t>(nbins_z, 1), 0, nbins_x * nbins_y * max<size_t>(nbins_z, 1));
       h_error = h;
       for(size_t i=1; i <= nbins_x; ++i){
          for(size_t j=1; j<=nbins_y; ++j){
              // to treat both 2D and 3D histos, use k==0 as z-index for 2D histos,
              // but skip k==0 for 3D histos, as then, this is an underflow bin
              for(size_t k=0; k<=nbins_z; ++k){
                  if(k==0 && nbins_z > 0)continue;
                  double content = histo->GetBinContent(i, j, k);
                  size_t ibin = (i-1) * nbins_y * max<size_t>(nbins_z, 1) + (j-1) * max<size_t>(nbins_z, 1) + k;
                  if(nbins_z==0) ++ibin;
                  h.set(ibin, content);
                  if(use_errors && content > 0.0){
                      h_error.set(ibin, histo->GetBinError(i, j, k) / content);
                  }
              }
          }
       }
        
    }
    
    //apply zerobin_fillfactor:
    double zerobin_fillfactor = 0.0;
    if(ctx.setting.exists("zerobin_fillfactor")){
        zerobin_fillfactor = ctx.setting["zerobin_fillfactor"];
        if(zerobin_fillfactor < 0){
           throw ConfigurationException("zerobin_fillfactor must be >= 0.0!");
        }
    }
    double integral = h.get_sum_of_bincontents();
    for(size_t i=1; i<=h.get_nbins(); ++i){
       h.set(i, max(h.get(i), integral * zerobin_fillfactor));
    }
    
    set_histos(h, h_error);
}

REGISTER_PLUGIN(root_histogram)

