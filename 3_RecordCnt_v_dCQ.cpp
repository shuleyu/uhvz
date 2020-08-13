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

const double dCQThreshold = -1;

const double maxThickness = 50;

const double dRhoMin = -10, dRhoMax = 20;

const bool plotDifferentRanks = false; // false will ignore dCQThreshold:

const string modelingTable = "gen2CA_D.ModelingResult_Subtract";
const string binTable = "gen2CA_D.Bins";
const string propertyTable = "gen2CA_D.Properties";

// ------------------------------------

struct BinResult{

    int binN, traceCnt;

    double lon, lat;

    // Best fit model.
    double ulvzH = -1, ulvzdVs = 0, ulvzdRho = 0;

    double uhvzH = -1, uhvzdVs = 0, uhvzdRho = 0;

    double lamellaH = -1, lamelladVs = 0, lamelladRho = 0, lamellaAway = -1;

    double h = -1, dvs = 0, drho = 0, away = -1;

    // compare quality.
    double ulvzCQ = 0, premCQ = 0, uhvzCQ = 0, lamellaCQ = 0, cq = 0, dcq = 0;

    BinResult (int num,double la,double lo){

        binN=num;
        lat=la;
        lon=lo;
    }

};

int main(){


    // For each bin, get the location and bin number.
    auto binInfo = MariaDB::Select("bin,lon,lat from "+binTable);

    vector<BinResult> Data;
    for (size_t i=0;i<binInfo.NRow();++i) {
        Data.push_back(BinResult(binInfo.GetInt("bin")[i],binInfo.GetDouble("lat")[i],binInfo.GetDouble("lon")[i]));
    }


    // Get model properties.
    auto modelInfo=MariaDB::Select("modelName, thickness, dvs, drho, awayFromCMB from "+propertyTable+" order by modelName");


    // For each bin, get the best fit model.
    for (size_t i=0;i<Data.size();++i) {

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


            // Already done for this bin.
            if (Data[i].ulvzH >= 0 && Data[i].uhvzH >= 0 && Data[i].lamellaH >= 0) {
                break;
            }

            string modelName = res.GetString("modelName")[j];

            // Find the same kind before (since we alread sort by cq desc, no need to carry on.)
            if ( (modelName.find("UHVZ") != string::npos && Data[i].uhvzH >= 0)       ||
                 (modelName.find("ULVZ") != string::npos && Data[i].ulvzH >= 0)       ||
                 (modelName.find("Lamella") != string::npos && Data[i].lamellaH >=0 ) ){
                continue;
            }


            // get model property.
            auto it = lower_bound(modelInfo.GetString("modelName").begin(), modelInfo.GetString("modelName").end(), modelName);

            if (it == modelInfo.GetString("modelName").end() || *it != res.GetString("modelName")[j]) {

                cerr << "Can't find model property for " << res.GetString("modelName")[j];
            }
            else {

                size_t index = distance(modelInfo.GetString("modelName").begin(), it);

                if ( dRhoMin <= modelInfo.GetDouble("drho")[index] && 
                     modelInfo.GetDouble("drho")[index] <= dRhoMax &&
                     modelInfo.GetDouble("thickness")[index] <= maxThickness ){


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


                    // Update best fit model (doesn't count Lamella and PREM).
                    if (modelName.find("Lamella") == string::npos && modelName.find("PREM") == string::npos) {

                        double oldCQ = Data[i].dcq;
                        Data[i].dcq = max(Data[i].ulvzCQ, Data[i].uhvzCQ) - Data[i].premCQ;

                        if (oldCQ < Data[i].dcq) {

                            Data[i].cq = res.GetDouble("cq")[j];
                            Data[i].h = modelInfo.GetDouble("thickness")[index];
                            Data[i].dvs = modelInfo.GetDouble("dvs")[index];
                            Data[i].drho = modelInfo.GetDouble("drho")[index];
                            Data[i].away = modelInfo.GetDouble("awayFromCMB")[index];

                            // Get the waveform data for this bin.
                            Data[i].traceCnt=res.GetInt("stackTraceCnt")[j];
                        }
                    }
                }
            }
        }
    }

    // Some ways to sort the bins.

    // Plot
    double XSIZE = 15, YSIZE = 15;
    string outfile = GMT::BeginEasyPlot(XSIZE,YSIZE);
    GMT::MoveReferencePoint(outfile,"-X1i -Y1i");
    GMT::psbasemap(outfile,"-JX10i/10i -R0/1/0/2000 -Bxa0.25g0.25 -Bya250g250 -Bx+ldCQ -By+l\"number of record\" -BWS -O -K");

    // X,Y Plot.

    for (auto &item: Data) {

        GMT::psxy(outfile,vector<double> {item.dcq}, vector<double> {double(item.traceCnt)},"-J -R -Sc0.10i -Gred -O -K");
    }

    for (auto &item: Data) {

        vector<GMT::Text> texts;
        texts.push_back(GMT::Text(item.dcq, item.traceCnt + 25, to_string(item.binN), 4,"CM"));
        GMT::pstext(outfile,texts,"-J -R -N -O -K");
    }

    // Make PDF.
    GMT::SealPlot(outfile);
    GMT::ps2pdf(outfile, __FILE__);

    return 0;
}
