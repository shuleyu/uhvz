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

// Inputs. -----------------------------

const double dCQThreshold = -0.1;
const bool plotData = false; // plot dataStack or modelStack.
const bool plotDifferentRanks = false;
const vector<double> plotRanks{0.15, 0.20, dCQThreshold};

const double dRhoMin=0,dRhoMax=10;
const double maxThickness=50;

const string modelingTable="gen2CA_D.ModelingResult_Decon";
const string binTable="gen2CA_D.Bins";
const string propertyTable="gen2CA_D.Properties";

const double plotBinSize=1.5; // diameter of the bin size on the figure. (in degree.)
const double plotHeight=0.7;

// ------------------------------------


struct BinResult{
    int binN, traceCnt;
    double lon,lat;

    // Best fit model.
    string BestFitModel;
    double h=0,dvs=0,drho=0;

    // compare quality.
    double cq=0,premCQ=0,dcq=0;

    // Waveforms.
    EvenSampledSignal dataStack, dataStackStd, premStack, modelStack, modelStackStd, dataAlteredPremStack, modelAlteredPremStack, dataFR, modelFR;

    BinResult (int num,double la,double lo){
        binN=num;
        lat=la;
        lon=lo;
    }
};

ostream &operator<<(ostream &os, const vector<BinResult> &bins){
    os << "BinNumber, BestFitModel, CompareQuality(CQ), DifferenceCompareQuality(dCQ), Thickness(km), dVs(%), drho(%), AboveCMB(km)\n";
    for (const auto bin:bins)
        os << bin.binN << " " << bin.BestFitModel << " " << bin.cq << " " << bin.premCQ << " " << bin.h << " " << bin.dvs << " " << bin.drho << '\n';
    return os;
}

int main(){


    // For each bin, get the location and bin number.
    auto binInfo=MariaDB::Select("bin,lon,lat from "+binTable);
    vector<BinResult> Data;
    for (size_t i=0;i<binInfo.NRow();++i)
        Data.push_back(BinResult(binInfo.GetInt("bin")[i],binInfo.GetDouble("lat")[i],binInfo.GetDouble("lon")[i]));


    // Get model properties.
    auto modelInfo=MariaDB::Select("modelName, thickness, dvs, drho from "+propertyTable+" order by modelName");


    // For each bin, get the best fit model.
    for (size_t i=0;i<Data.size();++i) {

        // ModelName is in the form: "ModelType_2015xxx"
        auto res=MariaDB::Select("modelName, cq, concat(dirPrefix,'/',dataStack) as ds, concat(dirPrefix,'/',dataStackStd) as dss, concat(dirPrefix,'/',premStack) as ps, concat(dirPrefix,'/',modelStack) as ms, concat(dirPrefix,'/',modelStackStd) as mss, concat(dirPrefix,'/',dataAlteredPremStack) as daps, concat(dirPrefix,'/',modelAlteredPremStack) as maps, stackTraceCnt from "+modelingTable+" where bin="+to_string(Data[i].binN)+" order by cq desc");



        // Find prem result.
        for (size_t j=0; j<res.NRow(); ++j) {
            if (res.GetString("modelName")[j]=="PREM_201500000000"){
                Data[i].premCQ=res.GetDouble("cq")[j];
                break;
            }
        }


        // Find best fit model.
        for (size_t j=0; j<res.NRow(); ++j) {
            
            // get model property.
            auto it=lower_bound(modelInfo.GetString("modelName").begin(), modelInfo.GetString("modelName").end(),res.GetString("modelName")[j]);
            if (it==modelInfo.GetString("modelName").end() || *it!=res.GetString("modelName")[j])
                cerr << "Can't find model property for " << res.GetString("modelName")[j];
            else {
                size_t index=distance(modelInfo.GetString("modelName").begin(), it);
                if (dRhoMin<=modelInfo.GetDouble("drho")[index] && modelInfo.GetDouble("drho")[index]<=dRhoMax && modelInfo.GetDouble("thickness")[index]<=maxThickness) {
                    Data[i].h=modelInfo.GetDouble("thickness")[index];
                    Data[i].dvs=modelInfo.GetDouble("dvs")[index];
                    Data[i].drho=modelInfo.GetDouble("drho")[index];
                    Data[i].cq=res.GetDouble("cq")[j];
                    Data[i].dcq=Data[i].cq-Data[i].premCQ;


                    // Get the waveform data for this bin.
                    Data[i].dataStack=EvenSampledSignal(res.GetString("ds")[j]);
                    Data[i].dataStackStd=EvenSampledSignal(res.GetString("dss")[j]);
                    Data[i].premStack=EvenSampledSignal(res.GetString("ps")[j]);
                    Data[i].modelStack=EvenSampledSignal(res.GetString("ms")[j]);
                    Data[i].modelStackStd=EvenSampledSignal(res.GetString("mss")[j]);
                    Data[i].dataAlteredPremStack=EvenSampledSignal(res.GetString("daps")[j]);
                    Data[i].modelAlteredPremStack=EvenSampledSignal(res.GetString("maps")[j]);
                    Data[i].traceCnt=res.GetInt("stackTraceCnt")[j];

                    //
                    Data[i].dataStackStd.CheckAndCutToWindow(-29,29);
                    Data[i].modelStackStd.CheckAndCutToWindow(-29,29);

                    Data[i].dataFR=Data[i].dataStack-Data[i].dataAlteredPremStack;
                    Data[i].dataFR.Mask(0,30);
                    Data[i].dataFR.FlipReverseSum(0);

                    Data[i].modelFR=Data[i].modelStack-Data[i].modelAlteredPremStack;
                    Data[i].modelFR.Mask(0,30);
                    Data[i].modelFR.FlipReverseSum(0);

                    break;
                }
            }
        }

        if (Data[i].dcq<=0) cout << "Bin " << Data[i].binN << " prem fit best ... Interesting ..." << endl;

    }

    cout << Data;

    // Some ways to sort the bins.

    // Seperate high velocity from low velocity. Put high velocity group in front.
    // In each group, sort by bin numbering increasing.
    auto order_HighDvsGroup_bin=[](const BinResult &b1, const BinResult &b2){
        if (b1.dvs>0 && b2.dvs>0) return b1.binN<b2.binN;
        else if (b1.dvs==0 && b2.dvs==0) return b1.binN<b2.binN;
        else if (b1.dvs<0 && b2.dvs<0) return b1.binN<b2.binN;
        else if (b1.dvs==0) return b2.dvs<0;
        else if (b2.dvs==0) return b1.dvs>0;
        else return b1.dvs>b2.dvs;
    };

    // Seperate high velocity from low velocity. Put low velocity group in front.
    // In each group, sort by bin numbering increasing.
    auto order_LowDvsGroup_bin=[](const BinResult &b1, const BinResult &b2){
        if (b1.dvs>0 && b2.dvs>0) return b1.binN<b2.binN;
        else if (b1.dvs==0 && b2.dvs==0) return b1.binN<b2.binN;
        else if (b1.dvs<0 && b2.dvs<0) return b1.binN<b2.binN;
        else if (b1.dvs==0) return b2.dvs>0;
        else if (b2.dvs==0) return b1.dvs<0;
        else return b1.dvs<b2.dvs;
    };

    // Seperate high velocity from low velocity. Put high velocity group in front.
    // In each group, sort by dcq decreasing.
//     auto order_HighDvsGroup_dcq=[](const BinResult &b1, const BinResult &b2){
//         if (b1.dvs>0 && b2.dvs>0) return b1.dcq>b2.dcq;
//         else if (b1.dvs==0 && b2.dvs==0) return b1.dcq>b2.dcq;
//         else if (b1.dvs<0 && b2.dvs<0) return b1.dcq>b2.dcq;
//         else if (b1.dvs==0) return b2.dvs<0;
//         else if (b2.dvs==0) return b1.dvs>0;
//         else return b1.dvs>b2.dvs;
//     };

    // Seperate high velocity from low velocity. Put low velocity group in front.
    // In each group, sort by dcq decreasing.
//     auto order_LowDvsGroup_dcq=[](const BinResult &b1, const BinResult &b2){
//         if (b1.dvs>0 && b2.dvs>0) return b1.dcq>b2.dcq;
//         else if (b1.dvs==0 && b2.dvs==0) return b1.dcq>b2.dcq;
//         else if (b1.dvs<0 && b2.dvs<0) return b1.dcq>b2.dcq;
//         else if (b1.dvs==0) return b2.dvs>0;
//         else if (b2.dvs==0) return b1.dvs<0;
//         else return b1.dvs<b2.dvs;
//     };

    // By dcq ascending
    auto order_cq=[](const BinResult &b1, const BinResult &b2){
        return b1.dcq<b2.dcq;
    };


    // Count how many significant bins for each group.
    size_t Count_Significant_Low=0,Count_Significant_High=0;
    for (const auto &item:Data)
        if (item.dcq>dCQThreshold && item.dvs<=0)
            ++Count_Significant_Low;
        else if (item.dcq>dCQThreshold && item.dvs>0)
            ++Count_Significant_High;


    // Plot
    double scale=0.25,XSIZE=scale*35+12,YSIZE=1+max(scale*35 + 9 , 1.5+plotHeight*max(Count_Significant_Low,Count_Significant_High));
    string outfile=GMT::BeginEasyPlot(XSIZE,YSIZE,ShellExec("pwd",true)+"/"+string(__FILE__));


    // Plot low velocity Stacks.
    sort(Data.begin(),Data.end(),order_LowDvsGroup_bin);
    GMT::MoveReferencePoint(outfile,"-Xf1i -Yf"+to_string(YSIZE-1)+"i");
    for (auto &item:Data) {

        if (item.dcq<dCQThreshold) continue;
        if (item.dvs>0.1) break;

        // plot stacks.
        GMT::MoveReferencePoint(outfile,"-Y-"+to_string(plotHeight)+"i");
        GMT::psxy(outfile,vector<double> {-15,15},vector<double> {0,0},"-JX2i/1.3i -R-15/15/-1/1 -W0.5p,black,. -O -K");
        GMT::psxy(outfile,vector<double> {-10,-5,0,5,10},vector<double> {0,0,0,0,0},"-J -R -Sy0.03i -W0.5p,red -O -K");
        
        if (plotData){
            GMT::psxy(outfile,item.dataStack,"-JX -R -W0.5p,black -O -K");
            GMT::psxy(outfile,item.dataStack+item.dataStackStd,"-JX -R -W0.5p,black,- -O -K");
            GMT::psxy(outfile,item.dataStack-item.dataStackStd,"-JX -R -W0.5p,black,- -O -K");
            GMT::psxy(outfile,item.premStack,"-JX -R -W0.5p,darkgreen -O -K");
            GMT::psxy(outfile,item.dataAlteredPremStack,"-JX -R -W0.5p,orange -O -K");
        }
        else {
            GMT::psxy(outfile,item.modelStack,"-JX -R -W0.5p,red -O -K");
            GMT::psxy(outfile,item.modelStack+item.modelStackStd,"-JX -R -W0.5p,black,- -O -K");
            GMT::psxy(outfile,item.modelStack-item.modelStackStd,"-JX -R -W0.5p,black,- -O -K");
            GMT::psxy(outfile,item.premStack,"-JX -R -W0.5p,darkgreen -O -K");
            GMT::psxy(outfile,item.modelAlteredPremStack,"-JX -R -W0.5p,orange -O -K");
        }

        vector<GMT::Text> texts;
        texts.push_back(GMT::Text(-15.5,0.0,"Bin "+to_string(item.binN),12,"RB"));
        texts.push_back(GMT::Text(-15.5,0.0,"# "+to_string(item.traceCnt),10,"RT"));
        GMT::pstext(outfile,texts,"-J -R -N -O -K");
    }
    GMT::MoveReferencePoint(outfile,"-Y-1i");
    GMT::psbasemap(outfile,"-J -R -Bxa5f2.5 -Bya0.5 -Bx+lsec. -By+lAmp. -BWS -O -K");

    GMT::MoveReferencePoint(outfile,"-Xf3.5i -Yf"+to_string(YSIZE-1)+"i");
    for (auto &item:Data) {

        if (item.dcq<dCQThreshold) continue;
        if (item.dvs>0.1) break;

        // plot stacks.
        GMT::MoveReferencePoint(outfile,"-Y-"+to_string(plotHeight)+"i");
        GMT::psxy(outfile,vector<double> {0,15},vector<double> {0,0},"-JX1i/1.3i -R0/15/-0.5/0.5 -W0.5p,black,. -O -K");
        GMT::psxy(outfile,vector<double> {5,10},vector<double> {0,0},"-J -R -Sy0.03i -W0.5p,red -O -K");

        GMT::psxy(outfile,item.dataFR,"-JX -R -W1p,black -O -K");
        GMT::psxy(outfile,item.modelFR,"-JX -R -W1p,red -O -K");

        vector<GMT::Text> texts;
        texts.push_back(GMT::Text(14.9,0.2,"CQ "+Float2String(item.cq,2)+", dCQ "+Float2String(item.dcq,2),6,"RB"));
        texts.push_back(GMT::Text(14.9,0.1,Float2String(item.h,0)+"km,"+Float2String(item.dvs,0)+"%dVs,"+Float2String(item.drho,1)+"%d@~\162@~",6,"RB"));
        GMT::pstext(outfile,texts,"-J -R -N -O -K");
    }
    GMT::MoveReferencePoint(outfile,"-Y-1i");
    GMT::psbasemap(outfile,"-J -R -Bxa5f2.5 -Bya0.5 -Bx+lsec. -By+lAmp. -BWS -O -K");


    // Plot high velocity Stacks.
    sort(Data.begin(),Data.end(),order_HighDvsGroup_bin);
    GMT::MoveReferencePoint(outfile,"-Xf6i -Yf"+to_string(YSIZE-1)+"i");
    for (auto &item:Data) {

        if (item.dcq<dCQThreshold) continue;
        if (item.dvs<0.1) break;

        // plot stacks.
        GMT::MoveReferencePoint(outfile,"-Y-"+to_string(plotHeight)+"i");
        GMT::psxy(outfile,vector<double> {-15,15},vector<double> {0,0},"-JX2i/1.3i -R-15/15/-1/1 -W0.5p,black,. -O -K");
        GMT::psxy(outfile,vector<double> {-10,-5,0,5,10},vector<double> {0,0,0,0,0},"-J -R -Sy0.03i -W0.5p,red -O -K");

        if (plotData) {
            GMT::psxy(outfile,item.dataStack,"-JX -R -W0.5p,black -O -K");
            GMT::psxy(outfile,item.dataStack+item.dataStackStd,"-JX -R -W0.5p,black,- -O -K");
            GMT::psxy(outfile,item.dataStack-item.dataStackStd,"-JX -R -W0.5p,black,- -O -K");
            GMT::psxy(outfile,item.premStack,"-JX -R -W0.5p,darkgreen -O -K");
            GMT::psxy(outfile,item.dataAlteredPremStack,"-JX -R -W0.5p,orange -O -K");
        }
        else {
            GMT::psxy(outfile,item.modelStack,"-JX -R -W0.5p,blue -O -K");
            GMT::psxy(outfile,item.modelStack+item.modelStackStd,"-JX -R -W0.5p,black,- -O -K");
            GMT::psxy(outfile,item.modelStack-item.modelStackStd,"-JX -R -W0.5p,black,- -O -K");
            GMT::psxy(outfile,item.premStack,"-JX -R -W0.5p,darkgreen -O -K");
            GMT::psxy(outfile,item.modelAlteredPremStack,"-JX -R -W0.5p,orange -O -K");
        }

        vector<GMT::Text> texts;
        texts.push_back(GMT::Text(-15.5,0,"Bin "+to_string(item.binN),12,"RB"));
        texts.push_back(GMT::Text(-15.5,0,"# "+to_string(item.traceCnt),10,"RT"));
        GMT::pstext(outfile,texts,"-J -R -N -O -K");
    }

    // Plot high velocity Stacks.
    GMT::MoveReferencePoint(outfile,"-Xf8.5i -Yf"+to_string(YSIZE-1)+"i");
    for (auto &item:Data) {

        if (item.dcq<dCQThreshold) continue;
        if (item.dvs<0.1) break;

        // plot stacks.
        GMT::MoveReferencePoint(outfile,"-Y-"+to_string(plotHeight)+"i");
        GMT::psxy(outfile,vector<double> {0,15},vector<double> {0,0},"-JX1i/1.3i -R0/15/-0.5/0.5 -W0.5p,black,. -O -K");
        GMT::psxy(outfile,vector<double> {5,10},vector<double> {0,0},"-J -R -Sy0.03i -W0.5p,red -O -K");

        GMT::psxy(outfile,item.dataFR,"-JX -R -W1p,black -O -K");
        GMT::psxy(outfile,item.modelFR,"-JX -R -W1p,blue -O -K");

        vector<GMT::Text> texts;
        texts.push_back(GMT::Text(14.9,0.2,"CQ "+Float2String(item.cq,2)+", dCQ "+Float2String(item.dcq,2),6,"RB"));
        texts.push_back(GMT::Text(14.9,0.1,Float2String(item.h,0)+"km,"+Float2String(item.dvs,0)+"%dVs,"+Float2String(item.drho,1)+"%d@~\162@~",6,"RB"));
        GMT::pstext(outfile,texts,"-J -R -N -O -K");
    }


    // Plot the map (bin locations).

    GMT::MoveReferencePoint(outfile,"-Xf11i -Yf"+to_string(YSIZE-scale*35)+"i");
    GMT::psbasemap(outfile,"-Jx"+to_string(scale)+"id/"+to_string(scale)+"id -R-99/-64/-5/27 -Bxa10f5 -Bya10f5 -Bx+lLongitude -By+lLatitude -BWSne -O -K");
    GMT::pscoast(outfile,"-J -R -O -A10000 -W0.5p,black -K");

    // Plot the bins.
    sort(Data.begin(),Data.end(),order_cq);


    // color scheme for significant bins.
    ofstream fpout("tmp.cpt");
    fpout << 
"-0.1   180/215/180    0.1    180/215/180\n"
"B  red\n"
"F  blue\n"
"N  127.5";
    fpout.close();

    for (size_t i=0;i<Data.size();++i) {

        // choose the bin boundary pen color.
        // using the category (high -> blue, low -> red).
        string PenColor="black";
        if (Data[i].dvs<0.1) PenColor="red";
        else if (Data[i].dvs>0.1) PenColor="blue";


        // Plot the bin as a circle with blue/red/cyan/fill to indicate the significance.
        double sizeWeight=1;
        if (plotDifferentRanks) {
            if (Data[i].dcq <= plotRanks[0]) {
                sizeWeight=0.25;
            }
            else if (Data[i].dcq <= plotRanks[1]) {
                sizeWeight=0.5;
            }
            else if (Data[i].dcq <= plotRanks[2]) {
                sizeWeight=0.75;
            }
            else {
                sizeWeight=1;
            }
        }
        else if (Data[i].dcq<dCQThreshold) {
            sizeWeight=0.5;
        }
        vector<vector<double>> plot_bin{{Data[i].lon, Data[i].lat, (fabs(Data[i].dcq)<dCQThreshold?0.0:(Data[i].dvs<0.1?-1.0:1.0)), scale * plotBinSize * sizeWeight}};

        if (plotDifferentRanks && Data[i].dcq <= plotRanks[0]) {
            GMT::psxy(outfile,plot_bin,"-J -R -Sx -Ctmp.cpt -W0.5p,"+PenColor+" -O -K");
        }
        else {
            GMT::psxy(outfile,plot_bin,"-J -R -Sc -Ctmp.cpt -W0.5p,"+PenColor+" -O -K");
        }


        // Add bin number to significant bins.
        if (Data[i].dcq>dCQThreshold){
            vector<GMT::Text> texts;
            texts.push_back(GMT::Text(Data[i].lon,Data[i].lat,"@;white;"+to_string(Data[i].binN)+"@;;",14,"CM","Helvetica-Bold"));
            GMT::pstext(outfile,texts,"-J -R -N -O -K");
        }
    }


    // Plot best fit model properties for the significant bins.
    vector<GMT::Text> texts;
    for (size_t i=0;i<Data.size();++i)
        if (Data[i].dcq>dCQThreshold)
            texts.push_back(GMT::Text{Data[i].dvs,Data[i].h,to_string(Data[i].binN),12,"CM"});


    GMT::MoveReferencePoint(outfile,"-Y-9i");
    GMT::psbasemap(outfile,"-Jx0.14i -R-35/25/0/55 -Bxa5g5f1 -Bya5g5f1 -Bx+l\"dVs (%)\" -By+l\"Thickness (km)\" -BWS -O -K");
    GMT::pstext(outfile,texts,"-J -R -O -K");


    // Make PDF.
    GMT::SealPlot(outfile);
    string pdffile=__FILE__;
    ShellExec("ps2pdf "+outfile+" "+pdffile.substr(0,pdffile.find_last_of("."))+".pdf");
    remove("tmp.cpt");
    remove(outfile.c_str());

    return 0;
}
