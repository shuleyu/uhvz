#include<iostream>
#include<map>
#include<queue>
#include<thread>
#include<atomic>
#include<mutex>
#include<condition_variable>

#include<MariaDB.hpp>
#include<EvenSampledSignal.hpp>
#include<ShellExec.hpp>
#include<PNormErr.hpp>
#include<RampFunction.hpp>
#include<GetHomeDir.hpp>

#include "CalculateCQ.hpp"

using namespace std;

mutex mtx;
condition_variable cv;
queue<size_t> emptySlot;

/*

This code stack the after-S-subtraction data/model traces for each bin-model pair.

For now, we assume the subtraction of synthetics S are always flat.

For each bin,

    select the correct distance
    stack the after-S-subtraction data (using same weighted).

    for each model,

        stack model (using same weight as data).

        make comparison between two stack results.
        calculate comparison quality for the all the bin-model pairs.
*/

// Inputs. ------------------------


const string targetModelType = "Lamella";         // "PREM" or "ULVZ" or "UHVZ" or "Lamella" or "All"(ignore beginIndex/endIndex)
const size_t beginIndex = 1, endIndex = 1584;
const bool reCreateTable = false;

const size_t nThread = 5;

const double distanceCutOff = 70; // To eiliminating Scd possiblility, do a hard distance cut-off.
const size_t cntThreshold = 20;
const double binEdgeWeight = 0.3, snrQuantile = 0.1, compareLen = 10;
const double weightSigma = sqrt(-1.0 / 2 / log(binEdgeWeight));

const string homeDir = GetHomeDir();
const string infoTable = "gen2CA_D.Master_a14";
const string dataTable = "gen2CA_D.Subtract";
const string binTable = "gen2CA_D.Bins";
const string binCenterDistTable = "gen2CA_D.CenterDists";
const string propertyTable = "gen2CA_D.Properties";
const string premTable = "REFL_PREM.Subtract";
const string ulvzTable = "REFL_ULVZ.Subtract";
const string uhvzTable = "REFL_UHVZ.Subtract";
const string lamellaTable = "REFL_Lamella.Subtract";

// Outputs. ------------------------

const string outputDB = "gen2CA_D";
const string outputTable = "ModelingResult_Subtract";
const string dirPrefix = homeDir + "/PROJ/t013.ScS_NextGen/Subtract";


// --------------------------------


void modelThese(size_t num, size_t mySlot,

                const string &modelName, const map<string, double> &criticalDistance,
                const vector<EvenSampledSignal> &dataWaveform, const map<string,size_t> &dataPairNameToIndex,
                const map<double,string> &gcarcSTNM,
                const vector<vector<string>> &binPairnames,
                const vector<vector<double>> &dataBinCenterDists, const vector<vector<double>> &dataBinGcarc,
                const vector<vector<double>> &dataBinSNR,
                const vector<double> &binRadius);

int main(){

    // Update table.
    if (reCreateTable) {

        MariaDB::Query("drop table if exists " + outputDB + "." + outputTable);
        MariaDB::Query("create table " + outputDB + "." + outputTable + " (pairname varchar(40) not null unique primary key, bin integer, modelName varchar(30), CQ double, CQ2 double comment \"Should be this?\", dataScSStack varchar(200), modelScSStack varchar(200), stackTraceCnt integer, weightSum double, dataScSStackStd varchar(200), modelScSStackStd varchar(200), dirPrefix varchar(200), index (bin), index(cq))");
    }


    // Make a model space in the format: "ULVZ_2015XXXXXXXX".
    // Get critical distances.
    vector<string> modelNames;
    if (targetModelType != "All") {
        for (size_t i = beginIndex; i <= endIndex; ++i) {

            modelNames.push_back(targetModelType + "_" + to_string(201500000000 + i));
        }
    }

    //
    auto critInfo = MariaDB::Select("modelName, criticalDist from " + propertyTable);
    map<string, double> criticalDistance;
    for (size_t i = 0; i < critInfo.NRow(); ++i) {
        criticalDistance[critInfo.GetString("modelName")[i]]=critInfo.GetDouble("criticalDist")[i];
    }
    if (modelNames.empty()) {
        modelNames = critInfo.GetString("modelName");
    }


    // Get data info, make a map between pairname and the index.
    // Read in data waveform, cut to -30 ~ 30 sec.
    auto dataInfo = MariaDB::Select("A.pairname as pn, concat(A.dirPrefix,'/',A.SStripped) as SFile, concat(A.dirPrefix,'/',A.ScSStripped) as ScSFile, B.eq as eq, B.stnm as stnm, B.shift_gcarc as shift_gcarc, B.SNR2_ScS as snr from " + dataTable + " as A join " + infoTable + " as B on A.pairname=B.pairname");

    vector<EvenSampledSignal> dataWaveform;
    map<string,size_t> dataPairNameToIndex;

    for (size_t i = 0; i < dataInfo.NRow(); ++i) {
        dataWaveform.push_back(EvenSampledSignal(dataInfo.GetString("ScSFile")[i]) - EvenSampledSignal(dataInfo.GetString("SFile")[i]));

        dataWaveform.back().CheckAndCutToWindow(-30, 30);
        dataWaveform.back().Mask(0, 30);
        dataWaveform.back().FlipReverseSum(0);
        dataPairNameToIndex[dataInfo.GetString("pn")[i]] = i;
    }


    // Make a map between gcarc and stnm (for synthetics selection)
    auto premInfo = MariaDB::Select("pairname, gcarc from " + premTable);
    map<double,string> gcarcSTNM;

    for (size_t i = 0; i < premInfo.NRow(); ++i) {
        gcarcSTNM[premInfo.GetDouble("gcarc")[i]] = premInfo.GetString("pairname")[i].substr(13);
    }


    // Get the list of data inside each bin
    // Get data distance to bin center for each bin.
    // Get bin radius for each bin.
    // Get the gcarc distances in each bin.
    auto binInfo = MariaDB::Select("bin, radius from " + binTable);

    vector<vector<string>> binPairnames;
    vector<vector<double>> dataBinCenterDists, dataBinGcarc, dataBinSNR;
    vector<double> binRadius;

    for (size_t i = 0; i < binInfo.NRow(); ++i) {
        const string binN = to_string(binInfo.GetInt("bin")[i]);
        auto binDataInfo = MariaDB::Select("pairname, dist_" + binN + " as centerDist from " + binCenterDistTable + " where dist_" + binN + " >= 0");

        binPairnames.push_back(binDataInfo.GetString("pairname"));
        dataBinCenterDists.push_back(binDataInfo.GetDouble("centerDist"));
        binRadius.push_back(binInfo.GetDouble("radius")[i]);

        dataBinGcarc.push_back(vector<double> ());
        dataBinSNR.push_back(vector<double> ());
        for (auto pn: binPairnames.back()) {
            dataBinGcarc.back().push_back(dataInfo.GetDouble("shift_gcarc")[dataPairNameToIndex[pn]]);
            dataBinSNR.back().push_back(dataInfo.GetDouble("snr")[dataPairNameToIndex[pn]]);
        }
    }


    // Start modeling.
    // For each model, for each bin, do modeling.
    vector<thread> allThreads(nThread);
    for (size_t i = 0; i < nThread; ++i) {
        emptySlot.push(i);
    }

    for (size_t runThisModel = 0; runThisModel < modelNames.size(); ++runThisModel){

        unique_lock<mutex> lck(mtx);
        while (emptySlot.empty()) {
            cv.wait(lck);
        }
        if (allThreads[emptySlot.front()].joinable()) {
            allThreads[emptySlot.front()].join();
        }
        allThreads[emptySlot.front()] = thread(modelThese, runThisModel, emptySlot.front(),

                                               cref(modelNames[runThisModel]), cref(criticalDistance),
                                               cref(dataWaveform), cref(dataPairNameToIndex),
                                               cref(gcarcSTNM),
                                               cref(binPairnames),
                                               cref(dataBinCenterDists), cref(dataBinGcarc), cref(dataBinSNR),
                                               cref(binRadius));


        emptySlot.pop();
    }

    for (auto &item: allThreads) {
        if (item.joinable()){
            item.join();
        }
    }

    return 0;
}

void modelThese(size_t num, size_t mySlot,

                const string &modelName, const map<string, double> &criticalDistance,
                const vector<EvenSampledSignal> &dataWaveform, const map<string,size_t> &dataPairNameToIndex,
                const map<double,string> &gcarcSTNM,
                const vector<vector<string>> &binPairnames,
                const vector<vector<double>> &dataBinCenterDists, const vector<vector<double>> &dataBinGcarc,
                const vector<vector<double>> &dataBinSNR,
                const vector<double> &binRadius){


    // total reflection distance for this model.

    const string modelEQ=modelName.substr(modelName.find("_")+1);
    const string modelType=modelName.substr(0,modelName.find("_"));
    const string modelTable=( modelType == "PREM" ? premTable : ( modelType == "ULVZ" ? ulvzTable : ( modelType == "UHVZ" ? uhvzTable : lamellaTable)));


    // distance selection.
    double critDist=criticalDistance.at(modelName);
    critDist = min(critDist, distanceCutOff);

    unique_lock<mutex> lck(mtx);
    cout << "Modeling against " << modelName << ", Num: " << num << " ... " << endl;
    auto modelInfo=MariaDB::Select("pairname, concat(dirPrefix,'/',ScSStripped) as fn from "+modelTable+" where eq="+modelEQ);

    // Read in model waveforms, cut to -30 ~ 30 sec.
    vector<EvenSampledSignal> modelWaveform;
    map<string,size_t> modelPairNameToIndex;

    for (size_t i = 0; i < modelInfo.NRow(); ++i) {
        modelWaveform.push_back(EvenSampledSignal(modelInfo.GetString("fn")[i]));
        if (!modelWaveform.back().CheckAndCutToWindow(-30,30)) {
cout << "Data corrupted: " << modelType << " " << modelEQ <<  " " << modelInfo.GetString("fn")[i] << endl ;
        }
        modelWaveform.back().Mask(0,30);
        modelWaveform.back().FlipReverseSum(0);
        modelPairNameToIndex[modelInfo.GetString("pairname")[i]]=i;
    }

    lck.unlock();

    // Model each bin.

    vector<double> weightSum(binRadius.size(),0), stackTraceCnt = weightSum, cqResult(binRadius.size(), 0.0/0.0), cqResult2 = cqResult;
    vector<string> dataScSStackFilename(binRadius.size()), modelScSStackFilename(binRadius.size()), dataScSStackStdFilename(binRadius.size()), modelScSStackStdFilename(binRadius.size());

    for (size_t i=0; i<binRadius.size(); ++i) {

        const string binN=to_string(i+1);

        // Get the SNR threshold using the quantile.
        vector<double> tmpArray;
        for (size_t j=0; j<binPairnames[i].size(); ++j) {
            if (dataBinGcarc[i][j]>=critDist) continue;
            tmpArray.push_back(dataBinSNR[i][j]);
        }
        sort(tmpArray.begin(),tmpArray.end());
        double critSNR=(tmpArray.empty()? -1 : tmpArray[(size_t)(tmpArray.size()*snrQuantile)]);

        // select the data waveform.
        // get the stack weight.
        vector<EvenSampledSignal> binDataWaveform, binModelWaveform;
        vector<double> binStackWeight;

        for (size_t j=0; j<binPairnames[i].size(); ++j) {

            binDataWaveform.push_back(dataWaveform[dataPairNameToIndex.at(binPairnames[i][j])]);

            if (dataBinGcarc[i][j]>=critDist) {
                binDataWaveform.pop_back();
                continue;
            }

            // find the correct distance synthetics data.

            auto it=gcarcSTNM.lower_bound(dataBinGcarc[i][j]);
            if (it==gcarcSTNM.end()) {
                it=prev(it);
            }
            string stnm=it->second;
            binModelWaveform.push_back(modelWaveform[modelPairNameToIndex[modelEQ+"_"+stnm]]);

            // Get weights.

            // 1. gaussian cap.
            binStackWeight.push_back(GaussianFunction(dataBinCenterDists[i][j]/binRadius[i], weightSigma, 0)*sqrt(2*M_PI)*weightSigma);

            // 2. SNR? or just a cut-off at the threshold.
            //if (dataBinSNR[i]<=critSNR) binStackWeight.back()=0;
            binStackWeight.back()*=RampFunction(dataBinSNR[i][j], 0, critSNR);
        }


        weightSum[i]=accumulate(binStackWeight.begin(),binStackWeight.end(),0.0);

        if (weightSum[i] <= 1 || binDataWaveform.size() < cntThreshold) {
            continue;
        }
        stackTraceCnt[i]=binDataWaveform.size();


        // Stack data, normalize stack and its std.
        auto binDataStack=StackSignals(binDataWaveform,binStackWeight);
        binDataStack.first.CheckAndCutToWindow(-29,29);
        binDataStack.second.CheckAndCutToWindow(-29,29);


        // Stack model, normalize stack and its std.
        auto binModelStack=StackSignals(binModelWaveform,binStackWeight);
        binModelStack.first.CheckAndCutToWindow(-29,29);
        binModelStack.second.CheckAndCutToWindow(-29,29);


        // Output to files.
        ShellExec("mkdir -p "+dirPrefix+"/dataScSStack/"+modelName+" "
                             +dirPrefix+"/modelScSStack/"+modelName);


        dataScSStackFilename[i]="dataScSStack/"+modelName+"/"+binN+".signal";
        dataScSStackStdFilename[i]="dataScSStack/"+modelName+"/"+binN+".std";
        modelScSStackFilename[i]="modelScSStack/"+modelName+"/"+binN+".signal";
        modelScSStackStdFilename[i]="modelScSStack/"+modelName+"/"+binN+".std";

        binDataStack.first.OutputToFile(dirPrefix+"/"+dataScSStackFilename[i]);
        binDataStack.second.OutputToFile(dirPrefix+"/"+dataScSStackStdFilename[i]);
        binModelStack.first.OutputToFile(dirPrefix+"/"+modelScSStackFilename[i]);
        binModelStack.second.OutputToFile(dirPrefix+"/"+modelScSStackStdFilename[i]);

        // Compare.
        auto compareResult = CalculateCQ(binDataStack.first, binModelStack.first, compareLen);
        cqResult[i] = compareResult[0] * compareResult[1];
        cqResult2[i] = compareResult[0] * compareResult[2];
    }

    // update database.

    lck.lock();

    vector<string> columnNames{"pairname", "bin", "modelName", "CQ", "CQ2", "dataScSStack", "modelScSStack","stackTraceCnt", "weightSum", "dataScSStackStd", "modelScSStackStd", "dirPrefix"};
    vector<vector<string>> sqlData(columnNames.size());

    for (size_t i=0; i<binRadius.size(); ++i) {
        sqlData[0].push_back(to_string(i+1)+"_"+modelName);
        sqlData[1].push_back(to_string(i+1));
        sqlData[2].push_back(modelName);
        sqlData[3].push_back(isnan(cqResult[i])?"NULL":to_string(cqResult[i]));
        sqlData[4].push_back(isnan(cqResult2[i])?"NULL":to_string(cqResult2[i]));
        sqlData[7].push_back(to_string(stackTraceCnt[i]));
        sqlData[8].push_back(to_string(weightSum[i]));
        sqlData[11].push_back(dirPrefix);
    }
    swap(sqlData[5], dataScSStackFilename);
    swap(sqlData[6], modelScSStackFilename);
    swap(sqlData[9], dataScSStackStdFilename);
    swap(sqlData[10], modelScSStackStdFilename);

    MariaDB::LoadData(outputDB, outputTable, columnNames, sqlData);

    emptySlot.push(mySlot);
    cv.notify_one();

    return;
}
