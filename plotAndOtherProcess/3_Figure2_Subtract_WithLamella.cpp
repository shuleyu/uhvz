#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include <unistd.h>
#include <iomanip>
#include <set>

#include <EvenSampledSignal.hpp>
#include <MariaDB.hpp>
#include <GMTPlotSignal.hpp>
#include <ShellExec.hpp>
#include <Float2String.hpp>

using namespace std;

// Inputs. -----------------------------

const double dCQThreshold = 0.45;

const double maxThickness = 50;
const double dRhoMin = -10,dRhoMax = 20;

const double plotBinSize = 1.5;
const double plotHeight = 0.7;

const string modelingTable = "gen2CA_D.ModelingResult_Subtract";
const string binTable = "gen2CA_D.Bins";
const string propertyTable = "gen2CA_D.Properties";

// ------------------------------------


struct BinResult {

    int binN, traceCnt;

    double lon, lat;

    // Best fit model.
    double ulvzH=-1, ulvzdVs=0, ulvzdRho=0;

    double uhvzH=-1, uhvzdVs=0, uhvzdRho=0;

    double lamellaH = -1, lamelladVs = 0, lamelladRho = 0, lamellaAway = -1;

    double h=-1, dvs=0, drho=0, away = -1;

    // compare quality.
    double ulvzCQ=0, premCQ=0, uhvzCQ=0, lamellaCQ = 0, cq=0, dcq=0;

    // Best fit waveforms.
    EvenSampledSignal dataStack, dataStackStd, modelStack, modelStackStd;

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
    for (size_t i=0;i<binInfo.NRow();++i) {

        Data.push_back(BinResult(binInfo.GetInt("bin")[i], binInfo.GetDouble("lat")[i], binInfo.GetDouble("lon")[i]));
    }


    // Get model properties.
    auto modelInfo = MariaDB::Select("modelName, thickness, dvs, drho, awayFromCMB from " + propertyTable + " order by modelName");


    // For each bin, get the best fit model.
    for (size_t i = 0; i < Data.size(); ++i) {

        // ModelName is in the form: "ModelType_2015xxx"
        auto res=MariaDB::Select("bin, modelName, cq, concat(dirPrefix,'/', dataScSStack) as dataFile, concat(dirPrefix,'/', modelScSStack) as modelFile , concat(dirPrefix,'/', dataScSStackStd) as dataFileStd, concat(dirPrefix,'/', modelScSStackStd) as modelFileStd , stackTraceCnt from " + modelingTable + " where bin=" + to_string(Data[i].binN) + " order by cq desc");


        // Find prem result.
        for (size_t j=0; j<res.NRow(); ++j) {
            if (res.GetString("modelName")[j]=="PREM_201500000000"){

                Data[i].premCQ=res.GetDouble("cq")[j];
                break;
            }
        }

        // Find best fit model.
        for (size_t j=0; j<res.NRow(); ++j) {


            // If all found, break;
            if (Data[i].ulvzH >= 0 && Data[i].uhvzH >= 0 && Data[i].lamellaH >= 0) {

                break;
            }

            string modelName = res.GetString("modelName")[j];

            if ( (modelName.find("UHVZ") != string::npos && Data[i].uhvzH >= 0)       ||
                 (modelName.find("ULVZ") != string::npos && Data[i].ulvzH >= 0)       ||
                 (modelName.find("Lamella") != string::npos && Data[i].lamellaH >=0 ) ){

                continue;
            }



            // get model property.
            auto it=lower_bound(modelInfo.GetString("modelName").begin(), modelInfo.GetString("modelName").end(), modelName);
            if (it==modelInfo.GetString("modelName").end() || *it!=res.GetString("modelName")[j]) {
                cerr << "Can't find model property for " << res.GetString("modelName")[j];
            }
            else {

                size_t index=distance(modelInfo.GetString("modelName").begin(), it);

                if (dRhoMin <= modelInfo.GetDouble("drho")[index] && modelInfo.GetDouble("drho")[index] <= dRhoMax && modelInfo.GetDouble("thickness")[index] <= maxThickness) {


                    if (modelName.find("ULVZ") != string::npos) {

                        Data[i].ulvzCQ = res.GetDouble("cq")[j];
                        Data[i].ulvzH = modelInfo.GetDouble("thickness")[index];
                        Data[i].ulvzdVs = modelInfo.GetDouble("dvs")[index];
                        Data[i].ulvzdRho = modelInfo.GetDouble("drho")[index];
                    }
                    else if (modelName.find("UHVZ") != string::npos) {

                        Data[i].uhvzCQ = res.GetDouble("cq")[j];
                        Data[i].uhvzH = modelInfo.GetDouble("thickness")[index];
                        Data[i].uhvzdVs = modelInfo.GetDouble("dvs")[index];
                        Data[i].uhvzdRho = modelInfo.GetDouble("drho")[index];
                    }
                    else if (modelName.find("Lamella") != string::npos) {

                        Data[i].lamellaCQ = res.GetDouble("cq")[j];
                        Data[i].lamellaH = modelInfo.GetDouble("thickness")[index];
                        Data[i].lamelladVs = modelInfo.GetDouble("dvs")[index];
                        Data[i].lamelladRho = modelInfo.GetDouble("drho")[index];
                        Data[i].lamellaAway = modelInfo.GetDouble("awayFromCMB")[index];
                    }


//if (modelName.find("Lamella") == string::npos && modelName.find("PREM") == string::npos) {
if (modelName.find("PREM") == string::npos) {


                    double oldCQ = Data[i].dcq;
                    Data[i].dcq = max(max(Data[i].ulvzCQ, Data[i].uhvzCQ), Data[i].lamellaCQ) - Data[i].premCQ;

                    if (oldCQ < Data[i].dcq) {

                        Data[i].cq = res.GetDouble("cq")[j];
                        Data[i].h = modelInfo.GetDouble("thickness")[index];
                        Data[i].dvs = modelInfo.GetDouble("dvs")[index];
                        Data[i].drho = modelInfo.GetDouble("drho")[index];
                        Data[i].away = modelInfo.GetDouble("awayFromCMB")[index];

                        // Get the waveform data for this bin.
                        Data[i].dataStack=EvenSampledSignal(res.GetString("dataFile")[j]);
                        Data[i].dataStackStd=EvenSampledSignal(res.GetString("dataFileStd")[j]);
                        Data[i].modelStack=EvenSampledSignal(res.GetString("modelFile")[j]);
                        Data[i].modelStackStd=EvenSampledSignal(res.GetString("modelFileStd")[j]);
                        Data[i].traceCnt=res.GetInt("stackTraceCnt")[j];
                        Data[i].dataStackStd.CheckAndCutToWindow(-29,29);
                        Data[i].modelStackStd.CheckAndCutToWindow(-29,29);
                    }
}

                }
            }
        }
    }

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
    // sort(Data.begin(),Data.end(),order_LowDvsGroup_dcq);
    GMT::MoveReferencePoint(outfile,"-Xf1i -Yf"+to_string(YSIZE-1)+"i");
    for (auto &item:Data) {

        if (item.dcq<dCQThreshold) continue;
        if (item.dvs>0.1) break;

        // plot stacks.
        GMT::MoveReferencePoint(outfile,"-Y-"+to_string(plotHeight)+"i");
        GMT::psxy(outfile,vector<double> {-15,15},vector<double> {0,0},"-JX2i/1.3i -R-15/15/-1/1 -W0.5p,black,. -O -K");
        GMT::psxy(outfile,vector<double> {-10,-5,0,5,10},vector<double> {0,0,0,0,0},"-J -R -Sy0.03i -W0.5p,red -O -K");


        GMT::psxy(outfile,item.dataStack,"-JX -R -W0.5p,black -O -K");
        GMT::psxy(outfile,item.dataStack+item.dataStackStd,"-JX -R -W0.5p,black,- -O -K");
        GMT::psxy(outfile,item.dataStack-item.dataStackStd,"-JX -R -W0.5p,black,- -O -K");

        GMT::psxy(outfile,item.modelStack,"-JX -R -W0.5p,red -O -K");
        GMT::psxy(outfile,item.modelStack+item.modelStackStd,"-JX -R -W0.5p,black,- -O -K");
        GMT::psxy(outfile,item.modelStack-item.modelStackStd,"-JX -R -W0.5p,black,- -O -K");

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
//         GMT::psxy(outfile,vector<double> {0,15},vector<double> {0,0},"-JX1i/1.3i -R0/15/-0.25/0.25 -W0.5p,black,. -O -K");
//         GMT::psxy(outfile,vector<double> {5,10},vector<double> {0,0},"-J -R -Sy0.03i -W0.5p,red -O -K");
//         GMT::psxy(outfile,item.dataStack,"-JX -R -W0.5p,black -O -K");
//         GMT::psxy(outfile,item.modelStack,"-JX -R -W0.5p,red -O -K");





        // signal to plot.
        vector<double> dataAmp, modelAmp;

        item.dataStack.CheckAndCutToWindow(0, 15);
        item.modelStack.CheckAndCutToWindow(0, 15);

        if (item.dataStack.GetAmp().size() != item.modelStack.GetAmp().size()) {
            throw runtime_error("Hmmm....");
        }   
        size_t npts = item.dataStack.GetAmp().size();

        for (size_t i = 0; i < npts; ++i) {
            dataAmp.push_back(item.dataStack.GetAmp()[npts - 1 - i] * -1);
            modelAmp.push_back(item.modelStack.GetAmp()[npts - 1 - i] * -1);
        }   

        EvenSampledSignal dataForPlot(dataAmp, item.dataStack.GetDelta(), -15);
        EvenSampledSignal modelForPlot(modelAmp, item.modelStack.GetDelta(), -15);

        GMT::psxy(outfile,vector<double> {-15, 0},vector<double> {0,0},"-JX0.543i/1.4118i -R-15/0/-0.5/0.5 -W0.5p,black -O -K");
//         GMT::psxy(outfile,vector<double> {-10, -5},vector<double> {0,0},"-JX0.543i/1.4118i -R-15/0/-0.5/0.5 -Sy0.03i -W0.5p,red -O -K");
        GMT::psxy(outfile,dataForPlot,"-JX0.543i/1.4118i -R-15/0/-0.5/0.5 -W0.5p,black -O -K");
        GMT::psxy(outfile,modelForPlot,"-JX -R -W0.5p,blue -O -K");











        vector<GMT::Text> texts;
        if (item.away > 0.5) {
            texts.push_back(GMT::Text(14.9,0.16,"CQ "+Float2String(item.cq,2)+", dCQ "+Float2String(item.dcq,2),6,"RB"));
            texts.push_back(GMT::Text(14.9,0.13,Float2String(item.h,0)+"km,"+Float2String(item.dvs,0)+"%dVs,"+Float2String(item.drho,1)+"%d@~\162@~",6,"RB"));
            texts.push_back(GMT::Text(14.9,0.1, "Dist to CMB: " + Float2String(item.away,1)+" km.",6,"RB"));
        }
        else {
            texts.push_back(GMT::Text(14.9,0.16,"CQ "+Float2String(item.cq,2)+", dCQ "+Float2String(item.dcq,2),6,"RB"));
            texts.push_back(GMT::Text(14.9,0.13,Float2String(item.h,0)+"km,"+Float2String(item.dvs,0)+"%dVs,"+Float2String(item.drho,1)+"%d@~\162@~",6,"RB"));
        }
        GMT::pstext(outfile,texts,"-J -R -N -O -K");
    }
    GMT::MoveReferencePoint(outfile,"-Y-1i");
    GMT::psbasemap(outfile,"-J -R -Bxa5f2.5 -Bya0.25 -Bx+lsec. -By+lAmp. -BWS -O -K");


    // Plot high velocity Stacks.
    sort(Data.begin(),Data.end(),order_HighDvsGroup_bin);
    // sort(Data.begin(),Data.end(),order_HighDvsGroup_dcq);
    GMT::MoveReferencePoint(outfile,"-Xf6i -Yf"+to_string(YSIZE-1)+"i");
    for (auto &item:Data) {

        if (item.dcq<dCQThreshold) continue;
        if (item.dvs<0.1) break;

        // plot stacks.
        GMT::MoveReferencePoint(outfile,"-Y-"+to_string(plotHeight)+"i");
        GMT::psxy(outfile,vector<double> {-15,15},vector<double> {0,0},"-JX2i/1.3i -R-15/15/-1/1 -W0.5p,black,. -O -K");
        GMT::psxy(outfile,vector<double> {-10,-5,0,5,10},vector<double> {0,0,0,0,0},"-J -R -Sy0.03i -W0.5p,red -O -K");

        GMT::psxy(outfile,item.dataStack,"-JX -R -W0.5p,black -O -K");
        GMT::psxy(outfile,item.dataStack+item.dataStackStd,"-JX -R -W0.5p,black,- -O -K");
        GMT::psxy(outfile,item.dataStack-item.dataStackStd,"-JX -R -W0.5p,black,- -O -K");

        GMT::psxy(outfile,item.modelStack,"-JX -R -W0.5p,blue -O -K");
        GMT::psxy(outfile,item.modelStack+item.modelStackStd,"-JX -R -W0.5p,black,- -O -K");
        GMT::psxy(outfile,item.modelStack-item.modelStackStd,"-JX -R -W0.5p,black,- -O -K");

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
//         GMT::psxy(outfile,vector<double> {0,15},vector<double> {0,0},"-JX1i/1.3i -R0/15/-0.25/0.25 -W0.5p,black,. -O -K");
//         GMT::psxy(outfile,vector<double> {5,10},vector<double> {0,0},"-J -R -Sy0.03i -W0.5p,red -O -K");
//         GMT::psxy(outfile,item.dataStack,"-JX -R -W0.5p,black -O -K");
//         GMT::psxy(outfile,item.modelStack,"-JX -R -W0.5p,blue -O -K");





        // signal to plot.
        vector<double> dataAmp, modelAmp;

        item.dataStack.CheckAndCutToWindow(0, 15);
        item.modelStack.CheckAndCutToWindow(0, 15);

        if (item.dataStack.GetAmp().size() != item.modelStack.GetAmp().size()) {
            throw runtime_error("Hmmm....");
        }   
        size_t npts = item.dataStack.GetAmp().size();

        for (size_t i = 0; i < npts; ++i) {
            dataAmp.push_back(item.dataStack.GetAmp()[npts - 1 - i] * -1);
            modelAmp.push_back(item.modelStack.GetAmp()[npts - 1 - i] * -1);
        }   

        EvenSampledSignal dataForPlot(dataAmp, item.dataStack.GetDelta(), -15);
        EvenSampledSignal modelForPlot(modelAmp, item.modelStack.GetDelta(), -15);

        GMT::psxy(outfile,vector<double> {-15, 0},vector<double> {0,0},"-JX0.543i/1.4118i -R-15/0/-0.5/0.5 -W0.5p,black -O -K");
//         GMT::psxy(outfile,vector<double> {-10, -5},vector<double> {0,0},"-JX0.543i/1.4118i -R-15/0/-0.5/0.5 -Sy0.03i -W0.5p,red -O -K");
        GMT::psxy(outfile,dataForPlot,"-JX0.543i/1.4118i -R-15/0/-0.5/0.5 -W0.5p,black -O -K");
        GMT::psxy(outfile,modelForPlot,"-JX -R -W0.5p,blue -O -K");























        vector<GMT::Text> texts;
        if (item.away > 0.5) {
            texts.push_back(GMT::Text(14.9,0.16,"CQ "+Float2String(item.cq,2)+", dCQ "+Float2String(item.dcq,2),6,"RB"));
            texts.push_back(GMT::Text(14.9,0.13,Float2String(item.h,0)+"km,"+Float2String(item.dvs,0)+"%dVs,"+Float2String(item.drho,1)+"%d@~\162@~",6,"RB"));
            texts.push_back(GMT::Text(14.9,0.1, "Dist to CMB: " + Float2String(item.away,1)+" km.",6,"RB"));
        }
        else {
            texts.push_back(GMT::Text(14.9,0.16,"CQ "+Float2String(item.cq,2)+", dCQ "+Float2String(item.dcq,2),6,"RB"));
            texts.push_back(GMT::Text(14.9,0.13,Float2String(item.h,0)+"km,"+Float2String(item.dvs,0)+"%dVs,"+Float2String(item.drho,1)+"%d@~\162@~",6,"RB"));
        }
        GMT::pstext(outfile,texts,"-J -R -N -O -K");
    }


    // Plot the map (bin locations).

    GMT::MoveReferencePoint(outfile,"-Xf11i -Yf"+to_string(YSIZE-scale*35)+"i");
    GMT::psbasemap(outfile,"-Jx"+to_string(scale)+"id/"+to_string(scale)+"id -R-99/-64/-5/27 -Bxa10f5 -Bya10f5 -Bx+lLongitude -By+lLatitude -BWSne -O -K");
    GMT::pscoast(outfile,"-J -R -O -A10000 -W0.3p,black -K");

    // Plot the bins.
    sort(Data.begin(),Data.end(),order_cq);



    for (size_t i=0;i<Data.size();++i) {


        double symbolSize = scale * plotBinSize * (Data[i].dcq >= dCQThreshold ? 1 : 0.3);

        if (Data[i].dcq >= dCQThreshold) {

            string penColor = (Data[i].drho >= 0 ? "black" : "50/255/50");

            if (Data[i].dvs < 0) {

                symbolSize *= 0.7;

                GMT::psxy(outfile,vector<double> {Data[i].lon}, vector<double> {Data[i].lat}, "-J -R -Sc" + to_string(symbolSize) + "i -Gred -W6p," + penColor + " -O -K");
            }
            else {
                if (Data[i].h < 45) {
                    symbolSize *= 0.75;
                }
                if (Data[i].h < 40) {
                    symbolSize *= 0.75;
                }
                if (Data[i].h < 35) {
                    symbolSize *= 0.75;
                }
                GMT::psxy(outfile,vector<double> {Data[i].lon}, vector<double> {Data[i].lat}, "-J -R -Sc" + to_string(symbolSize) + "i -G70/70/255 -W6p," + penColor + " -O -K");
            }
        }
        else {

            string penColor = (Data[i].dvs<0.1 ? "orange" : "70/70/255");

            GMT::psxy(outfile,vector<double> {Data[i].lon}, vector<double> {Data[i].lat}, "-J -R -Sx" + to_string(symbolSize) + "i -W3p," + penColor + " -O -K");
        }


//         // Add bin number to significant bins.
//         if (Data[i].dcq > dCQThreshold){
//             vector<GMT::Text> texts;
//             texts.push_back(GMT::Text(Data[i].lon,Data[i].lat,"@;white;"+to_string(Data[i].binN)+"@;;",14,"CM","Helvetica-Bold"));
//             GMT::pstext(outfile,texts,"-J -R -N -O -K");
//         }
    }


    // Plot best fit model properties for the significant bins.
    vector<double> px, py;
    vector<GMT::Text> texts;
    for (size_t i=0;i<Data.size();++i) {
        if (Data[i].dcq>dCQThreshold) {
            texts.push_back(GMT::Text{Data[i].dvs,Data[i].h,to_string(Data[i].binN),12,"CM"});

            px.push_back(Data[i].dvs);
            py.push_back(Data[i].h);
        }
    }


    GMT::MoveReferencePoint(outfile,"-Y-9i");
    GMT::psbasemap(outfile,"-Jx0.14i -R-35/25/0/55 -Bxa5g5f1 -Bya5g5f1 -Bx+l\"dVs (%)\" -By+l\"Thickness (km)\" -BWS -O -K");
    GMT::pstext(outfile,texts,"-J -R -O -K");
    GMT::psxy(outfile, px, py, "-J -R -Sc0.05i -Gblack -O -K");


    // Make PDF.
    GMT::SealPlot(outfile);
    string pdffile=__FILE__;
    ShellExec("ps2pdf "+outfile+" "+pdffile.substr(0,pdffile.find_last_of("."))+".pdf");
    remove(outfile.c_str());

    return 0;
}
