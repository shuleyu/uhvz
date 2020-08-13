#include<iostream>
#include<vector>
#include<queue>
#include<thread>
#include<atomic>
#include<algorithm>

#include<ShellExecVec.hpp>
#include<ShellExec.hpp>
#include<EvenSampledSignal.hpp>
#include<SACSignals.hpp>
#include<GMTPlotSignal.hpp>
#include<MariaDB.hpp>
#include<Float2String.hpp>
#include<GetHomeDir.hpp>

/**********************************************************************************************************
 *
 * Run this code on t039.PREM or t039.ULVZ or t039.UHVZ -- Use subtraction instead of deconvolution.
 *
 * The source is S stack then modified according to ScS Stack (for each events).
 *
 * This will generate subtracted ScS (*.ScSStripped) and FRSed (*.frs).
 *
 **********************************************************************************************************/

using namespace std;

mutex mtx;
condition_variable cv;
queue<size_t> emptySlot;

// Inputs. ------------------------------------

const string homeDir = GetHomeDir();
const string synDataDir = homeDir + "/PROJ/t039.UHVZ";
const string premDataDir = homeDir + "/PROJ/t039.PREM/201500000000";

const size_t beginIndex = 600, endIndex = 600, TraceCnt = 451; // [beginIndex, endIndex] inclusive.
const size_t nThread = 5;

const bool reCreateTable = false, makePlots = true;
const double plotIndex = 600;

const double filterCornerLow = 0.033, filterCornerHigh = 0.3, dt = 0.025;
const double cutSourceT1 = -100, cutSourceT2 = 100;
const double cutBeforeStripT1 = -100, cutBeforeStripT2 = 100;
const double cutResultT1 = -50, cutResultT2 = 50;


// Outputs. ------------------------------------

const string dirPrefix=homeDir+"/PROJ/t041.REFL_UHVZ/Subtract";
const string outputDB="REFL_UHVZ", outputTable="Subtract";


// --------------------------------------------

void processThis(const size_t Index, int mySlot, const EvenSampledSignal &sESW);

int main(){

    // Update table.
    if (reCreateTable) {
        MariaDB::Query("drop table if exists "+outputDB+"."+outputTable);
        MariaDB::Query("create table "+outputDB+"."+outputTable+" (PairName varchar(30) not null unique primary key, EQ varchar(20), Gcarc double, ScSStripped varchar(200), dirPrefix varchar(200))");
    }


    // Make ESW one time.

    /***********************************
    * 1. Read in PREM S waveform.      *
    ***********************************/

    SACSignals sESWData;

    sESWData=SACSignals (ShellExecVec("ls "+premDataDir+"/*.THT.sac"));

    if(sESWData.Size()!=TraceCnt) throw runtime_error("Reading error: PREM.");


    sESWData.SortByGcarc();
    sESWData.Interpolate(dt);
    sESWData.RemoveTrend();
    sESWData.HannTaper(20);
    sESWData.Butterworth(filterCornerLow,filterCornerHigh);


    // find S peak and shift time reference to the peak.
    sESWData.FindPeakAround(sESWData.GetTravelTimes("S"),10);
    sESWData.ShiftTimeReferenceToPeak();
    sESWData.FlipPeakUp();
    sESWData.NormalizeToPeak();

    /********************************************
    *
    * 2. Make S ESW stack.
    *
    *******************************************/

    sESWData.CheckAndCutToWindow(cutSourceT1-10,cutSourceT2+10);

    auto sESW=sESWData.XCorrStack(0,-15,15,2).second.first;
    sESW.FindPeakAround(0,10);
    sESW.ShiftTimeReferenceToPeak();
    sESW.CheckAndCutToWindow(cutSourceT1, cutSourceT2);
    sESW.FlipPeakUp();
    sESW.NormalizeToPeak();
    sESW.HannTaper(20);


    // Optional: make S ESW again, this time, stretch/shrink each S to match S ESW, then stack.

    sESWData.StretchToFit(sESW,-13,13,-0.3,0.3,0.25,true); // stretch, then compare waveform Amp_WinDiff.

    sESW=sESWData.XCorrStack(0,-15,15,2).second.first;
    sESW.FindPeakAround(0,10);
    sESW.ShiftTimeReferenceToPeak();
    sESW.CheckAndCutToWindow(cutSourceT1, cutSourceT2);
    sESW.FlipPeakUp();
    sESW.NormalizeToPeak();
    sESW.HannTaper(20);

    // Run the threads.

    vector<thread> allThreads(nThread);
    for (size_t i=0; i<nThread; ++i) {
        emptySlot.push(i);
    }

    for (size_t Index=0; Index<endIndex-beginIndex+1; ++Index) {
        unique_lock<mutex> lck(mtx);
        while (emptySlot.empty()) {
            cv.wait(lck);
        }
        if (allThreads[emptySlot.front()].joinable()) {
            allThreads[emptySlot.front()].join();
        }
        allThreads[emptySlot.front()] = thread(processThis, Index, emptySlot.front(), cref(sESW));
        emptySlot.pop();
    }

    for (auto &item: allThreads) {
        if (item.joinable()){
            item.join();
        }
    }

    return 0;
}

void processThis(const size_t Index, int mySlot, const EvenSampledSignal &sESW){


    /*************************************************
    * For each synthetic model...                    *
    * 1. Read in synthetic S and ScS waveform.       *
    * 2. Subtract S ESW from S waveform.             *
    * 3. Modifiy S ESW to match ScS waveforms.       *
    * 4. Strip modified S ESW from each ScS waveform.*
    * 5. Output waveforms.                           *
    * 6. UpdateTables. (Optional)                    *
    *************************************************/

    unique_lock<mutex> lck(mtx);
    const string modelName=to_string(201500000000+beginIndex+Index);

    // Properly modify the S ESW to look like ScS ESW: make the decon source.
    cout << "Processing synthetics: " << modelName << " ... " << endl;


    /*************************************
    * 1. Read in synthetic S waveform.   *
    *************************************/

    const string modelFolder=synDataDir+"/"+modelName;

    SACSignals Data;

    Data=SACSignals (ShellExecVec("ls "+modelFolder+"/*.THT.sac"));
    lck.unlock();

    if(Data.Size()!=TraceCnt) throw runtime_error("Reading error: "+modelName);


    Data.SortByGcarc();
    Data.Interpolate(dt);
    Data.RemoveTrend();
    Data.HannTaper(20);
    Data.Butterworth(filterCornerLow,filterCornerHigh);


    // find S peak and shift time reference to the peak.
    Data.FindPeakAround(Data.GetTravelTimes("S"),10);
    Data.ShiftTimeReferenceToPeak();
    Data.FlipPeakUp();
    Data.NormalizeToPeak();

    /************************************************
     *
     * 2. Subtract S waveforms.
     *
    ************************************************/

    SACSignals beforeSStrip=Data;

    // Properly modify the S ESW to look like ScS.
    vector<EvenSampledSignal> modifiedToFitS;
    for (size_t i=0; i<beforeSStrip.Size(); ++i) {
        //modifiedToFitS.push_back(sESW.StretchToFit(beforeSStrip.GetData()[i],-13,13,-0.3,0.3,0.25,true)); // stretch, then compare waveform Amp_WinDiff.
        //modifiedToFitScS.push_back(sESW.StretchToFit(beforeSStrip.GetData()[i],-13,13,-0.3,0.3,0.25,true,1)); // stretch, then compare waveform Amp_WinDiff.
        modifiedToFitS.push_back(sESW.StretchToFitHalfWidth(beforeSStrip.GetData()[i])); // stretch to fit the half-width.
    }

    auto SXCTimeShift=beforeSStrip.CrossCorrelation(-10,10,modifiedToFitS,-10,10).first;
    auto afterSStrip=beforeSStrip;
    afterSStrip.StripSignal(sESW,SXCTimeShift);

    beforeSStrip.CheckAndCutToWindow(cutBeforeStripT1, cutBeforeStripT2);


    /************************************************
     *
     * 3. Modifiy S ESW to match ScS waveforms.
     *
    ************************************************/

    SACSignals beforeScSStrip=afterSStrip;


    beforeScSStrip.FindPeakAround(beforeScSStrip.GetTravelTimes("ScS"),10);
    beforeScSStrip.ShiftTimeReferenceToPeak();
    beforeScSStrip.FlipPeakUp();
    beforeScSStrip.NormalizeToPeak();


    /************************************************
     *
     * 4. Strip modified S ESW from each ScS waveform.
     *
    ************************************************/


    // Properly modify the S ESW to look like ScS.
    vector<EvenSampledSignal> modifiedToFitScS;
    for (size_t i=0; i<beforeScSStrip.Size(); ++i) {
        //modifiedToFitScS.push_back(sESW.StretchToFit(beforeScSStrip.GetData()[i],-13,13,-0.3,0.3,0.25,true)); // stretch, then compare waveform Amp_WinDiff.
        //modifiedToFitScS.push_back(sESW.StretchToFit(beforeScSStrip.GetData()[i],-13,13,-0.3,0.3,0.25,true,1)); // stretch, then compare waveform Amp_WinDiff.
        modifiedToFitScS.push_back(sESW.StretchToFitHalfWidth(beforeScSStrip.GetData()[i])); // stretch to fit the half-width.
    }


    // Align at the best cross-correlation fit, before subtraction.
    // auto ScSXCTimeShift=beforeScSStrip.CrossCorrelation(-10,10,modifiedToFitScS,-10,10).first;

    // Align at peak before subtraction.
    vector<double> ScSXCTimeShift = vector<double> (beforeScSStrip.Size(), 0);

    auto afterScSStrip=beforeScSStrip;
    afterScSStrip.StripSignal(modifiedToFitScS,ScSXCTimeShift);

    for (size_t i=0; i<modifiedToFitScS.size(); ++i) {
        modifiedToFitScS[i].ShiftTime(ScSXCTimeShift[i]);
    }

    beforeScSStrip.CheckAndCutToWindow(cutBeforeStripT1, cutBeforeStripT2);
    afterScSStrip.CheckAndCutToWindow(cutResultT1,cutResultT2);

    /******************
     *
     * 5. Output.
     *
    ******************/

    // Output Decon waveforms.
    ShellExec("mkdir -p "+dirPrefix+"/"+modelName);
    afterScSStrip.DumpWaveforms(dirPrefix+"/"+modelName,"StationName","","","ScSStripped");



    // Plot
    if (makePlots && beginIndex+Index==plotIndex) {

        vector<string> outfiles;
        string outfile;

        for (size_t i=0; i<afterScSStrip.Size(); ++i) {

            if (i%17==0) { // A New page.
                outfile=GMT::BeginEasyPlot(28,40);
                outfiles.push_back(outfile);
                GMT::MoveReferencePoint(outfile,"-Xf1i -Yf37.2i");
            }
            else {
                GMT::MoveReferencePoint(outfile,"-Y-2.3i");
            }

            // Plot to verify sESW.
            GMT::psbasemap(outfile,"-JX13i/2i -R-100/100/-1/1 -Bxa10 -Bya0.5 -BWSne -O -K -Xf1i");
            GMT::psxy(outfile,vector<double> {-100,100},vector<double> {0,0},"-J -R -W1p,gray,- -O -K");
            GMT::psxy(outfile,sESW,"-J -R -W1p,yellow -O -K");
            GMT::psxy(outfile,Data, i, "-J -R -W1p,black -O -K");


            // Plot to verify the stripping.
            GMT::psbasemap(outfile,"-JX13i/2i -R-100/100/-1/1 -Bxa10 -Bya0.5 -BWSne -O -K -Xf14.5i");
            GMT::psxy(outfile,vector<double> {-100,100},vector<double> {0,0},"-J -R -W1p,gray,- -O -K");
            GMT::psxy(outfile,beforeScSStrip,i,"-J -R -W1p,black -O -K");
            GMT::psxy(outfile,modifiedToFitScS[i],"-J -R -W1p,yellow -O -K");
            GMT::psxy(outfile,afterScSStrip,i,"-J -R -W1p,green -O -K");
            vector<GMT::Text> texts{GMT::Text(-30,0.9,Float2String(afterScSStrip.GetMData()[i].gcarc,2),12,"LT")};
            GMT::pstext(outfile,texts,"-J -R -N -O -K");
        }

        for (auto file: outfiles) {
            GMT::SealPlot(file);
            ShellExec("cat "+file+" >> tmp.ps");
            remove(file.c_str());
        }
        string pdffile=__FILE__;
        ShellExec("ps2pdf tmp.ps "+pdffile.substr(0,pdffile.find_last_of("."))+".pdf");
        remove("tmp.ps");

    }

    // Finished. Update database and notify the main thread to continue;
    lck.lock();

    auto stationNames=Data.GetStationNames();
    auto gcarcs=Data.GetDistances();
    vector<vector<string>> sqlData(5,vector<string> ());
    for (size_t i=0; i<Data.Size(); ++i) {
        sqlData[0].push_back(modelName);
        sqlData[1].push_back(modelName+"_"+stationNames[i]);
        sqlData[2].push_back(Float2String(gcarcs[i],2));
        sqlData[3].push_back(modelName+"/"+stationNames[i]+".ScSStripped");
        sqlData[4].push_back(dirPrefix);
    }
    MariaDB::LoadData(outputDB,outputTable,vector<string> {"eq", "pairname", "gcarc", "ScSStripped", "dirPrefix"},sqlData);

    emptySlot.push(mySlot);
    cv.notify_one();

    return;
}
