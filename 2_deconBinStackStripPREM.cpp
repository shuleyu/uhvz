#include <iostream>
#include <map>
#include <queue>
#include <thread>
#include <atomic>

#include <MariaDB.hpp>
#include <EvenSampledSignal.hpp>
#include <ShellExec.hpp>
#include <PNormErr.hpp>
#include <RampFunction.hpp>
#include <GetHomeDir.hpp>

#include "CalculateCQ.hpp"

using namespace std;

mutex mtx;
condition_variable cv;
queue<size_t> emptySlot;

/*

This code stack the after-decon data/PREM/model traces for each bin.

The PREM stack is subtracted from the data and model stack for precursor comparison.

For each bin,

    select the correct distance
    stack the data and PREM (using same weighted).
    modify PREM or data to achive best fit (subjective step).
    subtract PREM from data.

    for each model,

        stack model (using same weight as data).
        modify PREM or model to achive best fit (subjective step).
        subtract PREM from model.

        make comparison between two subtraction results.
        calculate comparison quality for the all the bin-model pairs.

*/

// Inputs. ------------------------

const string targetModelType="Lamella";         // "PREM" or "ULVZ" or "UHVZ" or "Lamella" or "All"(ignore beginIndex/endIndex)
const size_t beginIndex=1, endIndex=1;
const bool reCreateTable=false;

const size_t nThread=5;

const size_t cntThreshold=20;
const double binEdgeWeight=0.3, snrQuantile=0.1, compareLen=15;
double weightSigma=sqrt(-1.0/2/log(binEdgeWeight));

const string homeDir=GetHomeDir();
const string dataTable="gen2CA_D.Master_a14";
const string binTable="gen2CA_D.Bins";
const string binCenterDistTable="gen2CA_D.CenterDists";
const string propertyTable="gen2CA_D.Properties";
const string premTable="REFL_PREM.Decon";
const string ulvzTable="REFL_ULVZ.Decon";
const string uhvzTable="REFL_UHVZ.Decon";
const string lamellaTable="REFL_Lamella.Decon";

// Outputs. ------------------------

const string outputDB="gen2CA_D";
const string outputTable="ModelingResult_Decon";
const string dirPrefix=homeDir+"/PROJ/t013.ScS_NextGen/Decon";


// --------------------------------

pair<EvenSampledSignal,double> matchHalfHeightWidth(const EvenSampledSignal &target, const EvenSampledSignal &varying);

void modelThese(size_t num, size_t mySlot,

                const string &modelName, const map<string, double> &criticalDistance,
                const vector<EvenSampledSignal> &dataWaveform, const map<string,size_t> &dataPairNameToIndex,
                const map<double,string> &gcarcSTNM,
                const vector<vector<string>> &binPairnames,
                const vector<vector<double>> &dataBinCenterDists, const vector<vector<double>> &dataBinGcarc,
                const vector<vector<double>> &dataBinSNR,
                const vector<double> &binRadius,
                const vector<EvenSampledSignal> &premWaveform, const map<string,size_t> premPairNameToIndex);

int main(){

    // Update table.
    if (reCreateTable) {
        MariaDB::Query("drop table if exists "+outputDB+"."+outputTable);
        MariaDB::Query("create table "+outputDB+"."+outputTable+" (pairname varchar(40) not null unique primary key, bin integer, modelName varchar(30), CQ double, CQ2 double comment \"Should be this?\", stackTraceCnt integer, weightSum double, dataAlterFactor double, modelAlterFactor double, dataStack varchar(200), modelStack varchar(200), premStack varchar(200), dataStackStd varchar(200), modelStackStd varchar(200), premStackStd varchar(200), dataAlteredPremStack varchar(200), modelAlteredPremStack varchar(200), premStrippedDataStack varchar(200), premStrippedModelStack varchar(200), dataFR varchar(200), modelFR varchar(200), dirPrefix varchar(200), index (bin), index(cq))");
    }


    // Make a model space in the format: "ULVZ_2015XXXXXXXX".
    // Get critical distances.
    vector<string> modelNames;
    if (targetModelType != "All") {
        for (size_t i=beginIndex; i<=endIndex; ++i) {
            modelNames.push_back(targetModelType + "_" + to_string(201500000000 + i));
        }
    }

    //
    auto critInfo=MariaDB::Select("modelName, criticalDist from "+propertyTable);
    map<string, double> criticalDistance;
    for (size_t i=0; i<critInfo.NRow(); ++i) {
        criticalDistance[critInfo.GetString("modelName")[i]]=critInfo.GetDouble("criticalDist")[i];
    }
    if (modelNames.empty()) {
        modelNames=critInfo.GetString("modelName");
    }


    // Get data info, make a map between pairname and the index.
    // Read in data waveform, cut to -30 ~ 30 sec.
    auto dataInfo=MariaDB::Select("pairname, eq, stnm, shift_gcarc, SNR2_ScS as snr, concat(dirPrefix,'/',DeconResult) as fn from "+dataTable);
    vector<EvenSampledSignal> dataWaveform;
    map<string,size_t> dataPairNameToIndex;

    size_t discardCnt=0;
    for (size_t i=0; i<dataInfo.NRow(); ++i) {
        dataWaveform.push_back(EvenSampledSignal(dataInfo.GetString("fn")[i]));
        dataWaveform.back().FindPeakAround(0);

        // Mark the waveform if its peak is changing (usually this is because on synthetics ScS is too close to SS).
        if (fabs(dataWaveform.back().PeakTime()) > dataWaveform.back().GetDelta()*1.5) {
            dataWaveform.back().SetTag(1);
            cout << "check decon of : " << dataInfo.GetString("eq")[i] << " - " << dataInfo.GetString("stnm")[i] << " peak now at: " << dataWaveform.back().PeakTime() << endl;
            ++discardCnt;
        }

        dataWaveform.back().ShiftTimeReferenceToPeak();
        dataWaveform.back().NormalizeToPeak();
        dataWaveform.back().CheckAndCutToWindow(-29.5,29.5);

        dataPairNameToIndex[dataInfo.GetString("pairname")[i]]=i;
    }

    cout << "---------- Will discard " << discardCnt << " / " << dataInfo.NRow() << " traces due to decon peak finding error..." << endl;


    // Make a map between gcarc and stnm (for synthetics selection)
    auto premInfo=MariaDB::Select("pairname, gcarc, concat(dirPrefix, '/', ScSDeconed) as fn from "+premTable);
    map<double,string> gcarcSTNM;

    for (size_t i=0; i<premInfo.NRow(); ++i) {
        gcarcSTNM[premInfo.GetDouble("gcarc")[i]]=premInfo.GetString("pairname")[i].substr(13);
    }


    // Get the list of data inside each bin
    // Get data distance to bin center for each bin.
    // Get bin radius for each bin.
    // Get the gcarc distances in each bin.
    auto binInfo=MariaDB::Select("bin, radius from "+binTable);

    vector<vector<string>> binPairnames;
    vector<vector<double>> dataBinCenterDists, dataBinGcarc, dataBinSNR;
    vector<double> binRadius;

    for (size_t i=0; i<binInfo.NRow(); ++i) {
        const string binN=to_string(binInfo.GetInt("bin")[i]);
        auto binDataInfo=MariaDB::Select("pairname, dist_"+binN+" as centerDist from "+binCenterDistTable+" where dist_"+binN+">=0");

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


    // Read in prem waveform, cut to -30 ~ 30 sec.
    vector<EvenSampledSignal> premWaveform;
    map<string,size_t> premPairNameToIndex;

    for (size_t i=0; i<premInfo.NRow(); ++i) {
        premWaveform.push_back(EvenSampledSignal(premInfo.GetString("fn")[i]));
        premWaveform.back().FindPeakAround(0);
        premWaveform.back().ShiftTimeReferenceToPeak();
        premWaveform.back().NormalizeToPeak();
        premWaveform.back().CheckAndCutToWindow(-29.5,29.5);

        premPairNameToIndex[premInfo.GetString("pairname")[i]]=i;
    }


    // Start modeling.
    // For each model, for each bin, do modeling.
    vector<thread> allThreads(nThread);
    for (size_t i=0; i<nThread; ++i) {
        emptySlot.push(i);
    }

    for (size_t runThisModel=0; runThisModel<modelNames.size(); ++runThisModel){

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
                                               cref(binRadius),
                                               cref(premWaveform), cref(premPairNameToIndex));

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
                const vector<double> &binRadius,
                const vector<EvenSampledSignal> &premWaveform, const map<string,size_t> premPairNameToIndex){


    // total reflection distance for this model.
    const double critDist=criticalDistance.at(modelName);

    const string modelEQ=modelName.substr(modelName.find("_")+1);
    const string modelType=modelName.substr(0,modelName.find("_"));
    const string modelTable=( modelType == "PREM" ? premTable : ( modelType == "ULVZ" ? ulvzTable : ( modelType == "UHVZ" ? uhvzTable : lamellaTable)));

    unique_lock<mutex> lck(mtx);
    cout << "Modeling against " << modelName << ", Num: " << num << " ... " << endl;
    auto modelInfo=MariaDB::Select("pairname, concat(dirPrefix,'/',ScSDeconed) as fn from "+modelTable+" where eq="+modelEQ);

    // Read in model waveforms, cut to -30 ~ 30 sec.
    vector<EvenSampledSignal> modelWaveform;
    map<string,size_t> modelPairNameToIndex;

    for (size_t i=0; i<modelInfo.NRow(); ++i) {
        modelWaveform.push_back(EvenSampledSignal(modelInfo.GetString("fn")[i]));
        if (!modelWaveform.back().CheckAndCutToWindow(-30,30)) {
cout << "Data corrupted: " << modelType << " " << modelEQ <<  " " << modelInfo.GetString("fn")[i] << endl ;
        }
        modelWaveform.back().FindPeakAround(0);
        modelWaveform.back().ShiftTimeReferenceToPeak();
        modelWaveform.back().NormalizeToPeak();
        modelWaveform.back().CheckAndCutToWindow(-29.5,29.5);

        modelPairNameToIndex[modelInfo.GetString("pairname")[i]]=i;
    }

    lck.unlock();


    // Model each bin.


    vector<double> weightSum(binRadius.size(),0), stackTraceCnt = weightSum, dataAlterFactor = weightSum, modelAlterFactor = weightSum, cqResult(binRadius.size(), 0.0/0.0), cqResult2 = cqResult;
    vector<string> dataStackFilename(binRadius.size()), modelStackFilename(binRadius.size()), premStackFilename(binRadius.size());
    vector<string> dataStackStdFilename(binRadius.size()), modelStackStdFilename(binRadius.size()), premStackStdFilename(binRadius.size());
    vector<string> dataAlteredPremStackFilename(binRadius.size()), modelAlteredPremStackFilename(binRadius.size());
    vector<string> premStrippedDataStackFileName(binRadius.size()), premStrippedModelStackFileName(binRadius.size()), dataFRFileName(binRadius.size()), modelFRFileName(binRadius.size());

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
        vector<EvenSampledSignal> binDataWaveform, binPremWaveform, binModelWaveform;
        vector<double> binStackWeight;

        for (size_t j=0; j<binPairnames[i].size(); ++j) {

            binDataWaveform.push_back(dataWaveform[dataPairNameToIndex.at(binPairnames[i][j])]);

            if (dataBinGcarc[i][j]>=critDist || binDataWaveform.back().GetTag()!=0) {
                binDataWaveform.pop_back();
                continue;
            }

            // find the correct distance synthetics data.

            auto it=gcarcSTNM.lower_bound(dataBinGcarc[i][j]);
            if (it==gcarcSTNM.end()) {
                it=prev(it);
            }
            string stnm=it->second;
            binPremWaveform.push_back(premWaveform[premPairNameToIndex.at("201500000000_"+stnm)]);
            binModelWaveform.push_back(modelWaveform[modelPairNameToIndex[modelEQ+"_"+stnm]]);


            // Get weights.

            // 1. gaussian cap.
            binStackWeight.push_back(GaussianFunction(dataBinCenterDists[i][j]/binRadius[i], weightSigma, 0)*sqrt(2*M_PI)*weightSigma);

            // 2. SNR? or just a cut-off at the threshold.
            //if (dataBinSNR[i][j]<=critSNR) binStackWeight.back()=0;
            binStackWeight.back()*=RampFunction(dataBinSNR[i][j],0,critSNR);
        }

        weightSum[i]=accumulate(binStackWeight.begin(),binStackWeight.end(),0.0);

        if (weightSum[i] <= 1 || binDataWaveform.size() < cntThreshold) {
            continue;
        }
        stackTraceCnt[i]=binDataWaveform.size();


        // Stack data, normalize stack and its std.
        auto binDataStack=StackSignals(binDataWaveform,binStackWeight);
        binDataStack.first.FindPeakAround(0,1);
        binDataStack.first.ShiftTimeReferenceToPeak();

        double amp=fabs(binDataStack.first.PeakAmp());
        binDataStack.second/=amp;
        binDataStack.first.NormalizeToPeak();

        // Stack prem, normalize stack and its std.
        auto binPremStack=StackSignals(binPremWaveform,binStackWeight);
        binPremStack.first.FindPeakAround(0);
        binPremStack.first.ShiftTimeReferenceToPeak();

        amp=fabs(binPremStack.first.PeakAmp());
        binPremStack.second/=amp;
        binPremStack.first.NormalizeToPeak();

        // Stack model, normalize stack and its std.
        auto binModelStack=StackSignals(binModelWaveform,binStackWeight);
        binModelStack.first.FindPeakAround(0,1);
        binModelStack.first.ShiftTimeReferenceToPeak();

        amp=fabs(binModelStack.first.PeakAmp());
        binModelStack.second/=amp;
        binModelStack.first.NormalizeToPeak();


        // Modify.

        auto dataRes=matchHalfHeightWidth(binDataStack.first, binPremStack.first);
        auto modelRes=matchHalfHeightWidth(binModelStack.first, binPremStack.first);

        auto &alteredDataPREM=dataRes.first;
        auto &alteredModelPREM=modelRes.first;

        dataAlterFactor[i]=dataRes.second;
        modelAlterFactor[i]=modelRes.second;

        binPremStack.first.CheckAndCutToWindow(-29,29);
        binPremStack.second.CheckAndCutToWindow(-29,29);
        binDataStack.first.CheckAndCutToWindow(-29,29);
        binDataStack.second.CheckAndCutToWindow(-29,29);
        binModelStack.first.CheckAndCutToWindow(-29,29);
        binModelStack.second.CheckAndCutToWindow(-29,29);
        alteredModelPREM.CheckAndCutToWindow(-29,29);
        alteredDataPREM.CheckAndCutToWindow(-29,29);

        // Output to files.
        ShellExec("mkdir -p "+dirPrefix+"/dataStack/"+modelName+" "
                             +dirPrefix+"/modelStack/"+modelName+" "
                             +dirPrefix+"/premStack/"+modelName+" "
                             +dirPrefix+"/modelAlteredPremStack/"+modelName+" "
                             +dirPrefix+"/dataAlteredPremStack/"+modelName+" "
                             +dirPrefix+"/premStrippedDataStack/"+modelName+" "
                             +dirPrefix+"/premStrippedModelStack/"+modelName+" "
                             +dirPrefix+"/dataFR/"+modelName+" "
                             +dirPrefix+"/modelFR/"+modelName+" "
                 );


        dataStackFilename[i]="dataStack/"+modelName+"/"+binN+".signal";
        dataStackStdFilename[i]="dataStack/"+modelName+"/"+binN+".std";
        modelStackFilename[i]="modelStack/"+modelName+"/"+binN+".signal";
        modelStackStdFilename[i]="modelStack/"+modelName+"/"+binN+".std";
        premStackFilename[i]="premStack/"+modelName+"/"+binN+".signal";
        premStackStdFilename[i]="premStack/"+modelName+"/"+binN+".std";
        dataAlteredPremStackFilename[i]="dataAlteredPremStack/"+modelName+"/"+binN+".signal";
        modelAlteredPremStackFilename[i]="modelAlteredPremStack/"+modelName+"/"+binN+".signal";

        binDataStack.first.OutputToFile(dirPrefix+"/"+dataStackFilename[i]);
        binDataStack.second.OutputToFile(dirPrefix+"/"+dataStackStdFilename[i]);
        binModelStack.first.OutputToFile(dirPrefix+"/"+modelStackFilename[i]);
        binModelStack.second.OutputToFile(dirPrefix+"/"+modelStackStdFilename[i]);
        binPremStack.first.OutputToFile(dirPrefix+"/"+premStackFilename[i]);
        binPremStack.second.OutputToFile(dirPrefix+"/"+premStackStdFilename[i]);
        alteredDataPREM.OutputToFile(dirPrefix+"/"+dataAlteredPremStackFilename[i]);
        alteredModelPREM.OutputToFile(dirPrefix+"/"+modelAlteredPremStackFilename[i]);


        // Strip prem from data.
        // Strip prem from model.
        auto dataFR=binDataStack.first-alteredDataPREM;
        auto modelFR=binModelStack.first-alteredModelPREM;


        // Output prem-subtraced bin data stack and bin model stack.

        premStrippedDataStackFileName[i]="premStrippedDataStack/"+modelName+"/"+binN+".signal";
        premStrippedModelStackFileName[i]="premStrippedModelStack/"+modelName+"/"+binN+".signal";
        dataFR.OutputToFile(dirPrefix+"/"+premStrippedDataStackFileName[i]);
        modelFR.OutputToFile(dirPrefix+"/"+premStrippedModelStackFileName[i]);


        // Flip and Reverse.

        dataFR.Mask(0,30);
        dataFR.FlipReverseSum(0);

        modelFR.Mask(0,30);
        modelFR.FlipReverseSum(0);


        dataFRFileName[i]="dataFR/"+modelName+"/"+binN+".signal";
        modelFRFileName[i]="modelFR/"+modelName+"/"+binN+".signal";
        dataFR.OutputToFile(dirPrefix+"/"+dataFRFileName[i]);
        modelFR.OutputToFile(dirPrefix+"/"+modelFRFileName[i]);


        // Compare.
        auto compareResult = CalculateCQ(dataFR, modelFR, compareLen);
        cqResult[i] = compareResult[0] * compareResult[1];
        cqResult2[i] = compareResult[0] * compareResult[2];

    }

    // Update table.

    lck.lock();

    vector<string> columnNames{"pairname", "bin", "modelName", "CQ", "CQ2","stackTraceCnt", "weightSum", "dataAlterFactor", "modelAlterFactor", "dataStack", "modelStack", "premStack", "dataStackStd", "modelStackStd", "premStackStd", "dataAlteredPremStack", "modelAlteredPremStack", "premStrippedDataStack", "premStrippedModelStack", "dataFR", "modelFR", "dirPrefix"};
    vector<vector<string>> sqlData(columnNames.size());

    for (size_t i=0; i<binRadius.size(); ++i) {
        sqlData[0].push_back(to_string(i+1)+"_"+modelName);
        sqlData[1].push_back(to_string(i+1));
        sqlData[2].push_back(modelName);
        sqlData[3].push_back(isnan(cqResult[i])?"NULL":to_string(cqResult[i]));
        sqlData[4].push_back(isnan(cqResult2[i])?"NULL":to_string(cqResult2[i]));
        sqlData[5].push_back(to_string(stackTraceCnt[i]));
        sqlData[6].push_back(to_string(weightSum[i]));
        sqlData[7].push_back(isnan(dataAlterFactor[i])?"NULL":to_string(dataAlterFactor[i]));
        sqlData[8].push_back(isnan(modelAlterFactor[i])?"NULL":to_string(modelAlterFactor[i]));
        sqlData[21].push_back(dirPrefix);
    }
    swap(sqlData[9],  dataStackFilename);
    swap(sqlData[10], modelStackFilename);
    swap(sqlData[11], premStackFilename);
    swap(sqlData[12], dataStackStdFilename);
    swap(sqlData[13], modelStackStdFilename);
    swap(sqlData[14], premStackStdFilename);
    swap(sqlData[15], dataAlteredPremStackFilename);
    swap(sqlData[16], modelAlteredPremStackFilename);
    swap(sqlData[17], premStrippedDataStackFileName);
    swap(sqlData[18], premStrippedModelStackFileName);
    swap(sqlData[19], dataFRFileName);
    swap(sqlData[20], modelFRFileName);

    MariaDB::LoadData(outputDB, outputTable, columnNames, sqlData);

    emptySlot.push(mySlot);
    cv.notify_one();

    return;
}

pair<EvenSampledSignal,double> matchHalfHeightWidth(const EvenSampledSignal &target, const EvenSampledSignal &varying){

    // Measure half-height width on target
    auto targetNptsLen=target.FindAmplevel(0.5);
    double targetHalfHeightWidth=target.GetDelta()*(targetNptsLen.second-targetNptsLen.first);

    // Measure half-height width on varying
    auto varyingNptsLen=varying.FindAmplevel(0.5);
    double varyingHalfHeightWidth=varying.GetDelta()*(varyingNptsLen.second-varyingNptsLen.first);


    // Shrink PREM when PREM is broader.
    if (varyingHalfHeightWidth>=targetHalfHeightWidth) {
        auto ans=varying*0;
        double sFactor=targetHalfHeightWidth/varyingHalfHeightWidth;
        ans.AddSignal(varying.Stretch(sFactor));
        return {ans, sFactor-1};
    }
    else {

        // When PREM is skinner,
        // either stretch, or gaussian PREM.

        // Choice 1. Stretch PREM (same code as above).
        auto ans=varying*0;
        double sFactor=targetHalfHeightWidth/varyingHalfHeightWidth;
        ans.AddSignal(varying.Stretch(sFactor));
        return {ans,sFactor-1};

        // Choice 2. Convolve varying with gaussian to find the best fit to target.
//         auto ans=varying;
//         double l=0.01, r=2, gFactor;
//
//         while (r-l>1e-3){
//
//             gFactor=l+(r-l)/2;
//
//             ans=varying;
//             ans.GaussianBlur(gFactor);
//             ans.FindPeakAround(0);
//             ans.ShiftTimeReferenceToPeak();
//             ans.NormalizeToPeak();
//
//             auto ansLen=ans.FindAmplevel(0.5);
//             double ansHalfHeightWidth=ans.GetDelta()*(ansLen.second-ansLen.first);
//
//             if (ansHalfHeightWidth==targetHalfHeightWidth) return {ans,gFactor};
//             else if (ansHalfHeightWidth<targetHalfHeightWidth) l=gFactor;
//             else r=gFactor;
//         }
//         return {ans,gFactor};
    }
}
