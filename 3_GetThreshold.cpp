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

#include<EvenSampledSignal.hpp>
#include<MariaDB.hpp>

using namespace std;

// Find the PREM FRS amplitude and use that number as a threshold.

// Inputs. -----------------------------

const string modelingTable="gen2CA_D.ModelingResult_test";
const string binTable="gen2CA_D.Bins";
const string propertyTable="gen2CA_D.Properties";
const double dRhoMin=0,dRhoMax=10;

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


    // Result.
    double maxAmp=0,minAmp=1;
    size_t maxBin=1,minBin=1;


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
                if (dRhoMin<=modelInfo.GetDouble("drho")[index] && modelInfo.GetDouble("drho")[index]<=dRhoMax) {
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

        if (maxAmp<Data[i].dataFR.MaxAmp()) {
            maxBin=Data[i].binN;
            maxAmp=Data[i].dataFR.MaxAmp();
        }

        if (minAmp>Data[i].dataFR.MaxAmp()) {
            minBin=Data[i].binN;
            minAmp=Data[i].dataFR.MaxAmp();
        }

        if (0.025>Data[i].dataFR.MaxAmp()) {
            cout << "Bin " << Data[i].binN << " has lower amplitude than threshold .." << endl;
        }

    }

    cout << "The maximum amplitude from FRS (as a threshold) is: " << maxAmp << " from bin " << maxBin << endl;
    cout << "The minimum amplitude from FRS (as a check) is: " << minAmp << " from bin " << minBin << endl;

    return 0;
}
