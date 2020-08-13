#include<iostream>
#include<cstdlib>
#include<fstream>
#include<sstream>
#include<vector>
#include<algorithm>
#include<cmath>
#include<string>
#include<unistd.h>
#include<iomanip>
#include<set>

#include<EvenSampledSignal.hpp>
#include<MariaDB.hpp>
#include<GMTPlotSignal.hpp>
#include<ShellExec.hpp>
#include<Float2String.hpp>

using namespace std;

// INPUTs -----------------------------


const string modelingTable="gen2CA_D.ModelingResult_Decon";
const string binTable="gen2CA_D.Bins";
const string propertyTable="gen2CA_D.Properties";

const double dCQThreshold=0.25;
const double dRhoMin=9,dRhoMax=11;

// ------------------------------------


struct BinResult{
    int binN;
    double lon,lat;

    // Best fit model.
    string BestFitModel;
    double h=0,dvs=0,drho=0;

    // compare quality.
    double cq=0,premCQ=0,dcq=0;

    // compare quality matrix.
    vector<vector<double>> cqs;

    BinResult (int num,double la,double lo){
        binN=num;
        lat=la;
        lon=lo;
    }
};

int main(){

    // For each bin, get the location and bin number.
    auto binInfo=MariaDB::Select("bin,lon,lat from "+binTable);
    vector<BinResult> Data;
    for (size_t i=0;i<binInfo.NRow();++i)
        Data.push_back(BinResult(binInfo.GetInt("bin")[i],binInfo.GetDouble("lat")[i],binInfo.GetDouble("lon")[i]));


    // Get model properties.
    auto modelInfo=MariaDB::Select("modelName, thickness, dvs, drho from "+propertyTable+" order by modelName");


    // For each bin, get the cq of all models.
    for (size_t i=0;i<Data.size();++i) {

        // ModelName is in the form: "ModelType_2015xxx"
        auto res=MariaDB::Select("modelName, cq from "+modelingTable+" where bin="+to_string(Data[i].binN));

        Data[i].cq=0;

        for (size_t j=0; j<res.NRow(); ++j) {

            // Find best fit model.
            if (Data[i].cq<res.GetDouble("cq")[j]){
                Data[i].cq=res.GetDouble("cq")[j];
                Data[i].BestFitModel=res.GetString("modelName")[j];
            }

            // Find prem cq.
            if (res.GetString("modelName")[j]=="PREM_201500000000")
                Data[i].premCQ=res.GetDouble("cq")[j];

            
            // get model property.
            auto it=lower_bound(modelInfo.GetString("modelName").begin(), modelInfo.GetString("modelName").end(),res.GetString("modelName")[j]);
            if (it==modelInfo.GetString("modelName").end() || *it!=res.GetString("modelName")[j])
                cerr << "Can't find model property for " << res.GetString("modelName")[j];
            else {
                size_t index=distance(modelInfo.GetString("modelName").begin(), it);
                if (dRhoMin<=modelInfo.GetDouble("drho")[index] && modelInfo.GetDouble("drho")[index]<=dRhoMax)
                    Data[i].cqs.push_back({modelInfo.GetDouble("dvs")[index],modelInfo.GetDouble("thickness")[index],res.GetDouble("cq")[j]});
                if (Data[i].BestFitModel==res.GetString("modelName")[j]) {
                    Data[i].dvs=modelInfo.GetDouble("dvs")[index];
                    Data[i].h=modelInfo.GetDouble("thickness")[index];
                    Data[i].drho=modelInfo.GetDouble("drho")[index];
                }
            }
        }

        // Fix the laziniess for dvs=0, drho=0 and thickness=1-30.
        if (dRhoMax<2) {
            for (double h=1; h<=30.5; h+=1) 
                Data[i].cqs.push_back({0,h,Data[i].premCQ});
        }

        Data[i].dcq=Data[i].cq-Data[i].premCQ;
        if (Data[i].dcq<=0) cout << "Bin " << Data[i].binN << " prem fit best ... Interesting ..." << endl;
    }

    // Some ways to sort the bins.

    // Seperate high velocity from low velocity. Put high velocity group in front.
    // In each group, sort by bin numbering increasing.
//     auto order_HighDvsGroup_bin=[](const BinResult &b1, const BinResult &b2){
//         if (b1.dvs>0 && b2.dvs>0) return b1.binN<b2.binN;
//         else if (b1.dvs==0 && b2.dvs==0) return b1.binN<b2.binN;
//         else if (b1.dvs<0 && b2.dvs<0) return b1.binN<b2.binN;
//         else if (b1.dvs==0) return b2.dvs<0;
//         else if (b2.dvs==0) return b1.dvs>0;
//         else return b1.dvs>b2.dvs;
//     };

    // Seperate high velocity from low velocity. Put low velocity group in front.
    // In each group, sort by bin numbering increasing.
//     auto order_LowDvsGroup_bin=[](const BinResult &b1, const BinResult &b2){
//         if (b1.dvs>0 && b2.dvs>0) return b1.binN<b2.binN;
//         else if (b1.dvs==0 && b2.dvs==0) return b1.binN<b2.binN;
//         else if (b1.dvs<0 && b2.dvs<0) return b1.binN<b2.binN;
//         else if (b1.dvs==0) return b2.dvs>0;
//         else if (b2.dvs==0) return b1.dvs<0;
//         else return b1.dvs<b2.dvs;
//     };

    // Seperate high velocity from low velocity. Put high velocity group in front.
    // In each group, sort by dcq decreasing.
    auto order_HighDvsGroup_dcq=[](const BinResult &b1, const BinResult &b2){
        if (b1.dvs>0 && b2.dvs>0) return b1.dcq>b2.dcq;
        else if (b1.dvs==0 && b2.dvs==0) return b1.dcq>b2.dcq;
        else if (b1.dvs<0 && b2.dvs<0) return b1.dcq>b2.dcq;
        else if (b1.dvs==0) return b2.dvs<0;
        else if (b2.dvs==0) return b1.dvs>0;
        else return b1.dvs>b2.dvs;
    };

    // Seperate high velocity from low velocity. Put low velocity group in front.
    // In each group, sort by dcq decreasing.
    auto order_LowDvsGroup_dcq=[](const BinResult &b1, const BinResult &b2){
        if (b1.dvs>0 && b2.dvs>0) return b1.dcq>b2.dcq;
        else if (b1.dvs==0 && b2.dvs==0) return b1.dcq>b2.dcq;
        else if (b1.dvs<0 && b2.dvs<0) return b1.dcq>b2.dcq;
        else if (b1.dvs==0) return b2.dvs>0;
        else if (b2.dvs==0) return b1.dvs<0;
        else return b1.dvs<b2.dvs;
    };

    // By dcq ascending
//     auto order_cq=[](const BinResult &b1, const BinResult &b2){
//         return b1.dcq<b2.dcq;
//     };


    // Count how many significant bins for each group.
    size_t Count_Significant_Low=0,Count_Significant_High=0;
    for (const auto &item:Data)
        if (item.dcq>dCQThreshold && item.dvs<=0)
            ++Count_Significant_Low;
        else if (item.dcq>dCQThreshold && item.dvs>0)
            ++Count_Significant_High;


    // Plot
    double height=1.7;
    double XSIZE=12,YSIZE=2+height*max(Count_Significant_Low,Count_Significant_High);
    string outfile=GMT::BeginEasyPlot(XSIZE,YSIZE,ShellExec("pwd",true)+"/"+string(__FILE__));
    GMT::set("COLOR_NAN cyan");
    GMT::makecpt("-Cgray -T0/1/0.1 -Z -I > tmp.cpt");


    // Plot low velocity Stacks.
    sort(Data.begin(),Data.end(),order_LowDvsGroup_dcq);
    GMT::MoveReferencePoint(outfile,"-Xf1i -Yf"+to_string(YSIZE-1)+"i");
    for (auto &item:Data) {


        if (item.dcq<dCQThreshold) continue;
        if (item.dvs>0.1) break;

        // plot cqs.
        GMT::MoveReferencePoint(outfile,"-Y-"+to_string(height)+"i");

        GMT::grdimage(outfile,item.cqs,1,1,"-JX4i/1i -R-30.5/20.5/0.5/30.5 -nn -Ctmp.cpt -O -K");

        vector<GMT::Text> texts;
        texts.push_back(GMT::Text(-30.5,34,"Bin "+to_string(item.binN),12,"RB"));
        GMT::pstext(outfile,texts,"-J -R -N -O -K");

        GMT::psbasemap(outfile,"-J -R -Bxa5 -Bya5 -Bx+ldvs(%) -By+lH(km.) -BWS -O -K");
    }

    GMT::psscale(outfile,"-Ctmp.cpt -D2i/-1i/4i/0.1ih -O -K -Bxa0.1 -Bx+l\"CQ\"");


    // Plot high velocity Stacks.
    sort(Data.begin(),Data.end(),order_HighDvsGroup_dcq);
    GMT::MoveReferencePoint(outfile,"-Xf6i -Yf"+to_string(YSIZE-1)+"i");
    for (auto &item:Data) {

        

        if (item.dcq<dCQThreshold) continue;
        if (item.dvs<0.1) break;

        // plot cqs.
        GMT::MoveReferencePoint(outfile,"-Y-"+to_string(height)+"i");

        GMT::grdimage(outfile,item.cqs,1,1,"-JX4i/1i -R-30.5/20.5/0.5/30.5 -nn -Ctmp.cpt -O -K");

        vector<GMT::Text> texts;
        texts.push_back(GMT::Text(-30.5,34,"Bin "+to_string(item.binN),12,"RB"));
        GMT::pstext(outfile,texts,"-J -R -N -O -K");

        GMT::psbasemap(outfile,"-J -R -Bxa5 -Bya5 -Bx+ldvs(%) -By+lThickness(km.) -BWS -O -K");
    }

    // Make PDF.
    GMT::SealPlot(outfile);
    string pdffile=__FILE__;
    ShellExec("ps2pdf "+outfile+" "+pdffile.substr(0,pdffile.find_last_of("."))+".pdf");
    remove("tmp.cpt");
    remove(outfile.c_str());

    return 0;
}
