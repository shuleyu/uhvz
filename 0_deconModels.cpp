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
 * Run this code on t039.PREM or t039.ULVZ or t039.UHVZ
 *
 * The source is S stack then modified according to ScS Stack (for each events).
 *
 * This will generate after-decon ScS (*.traces) and FRSed (*.frs).
 *
 **********************************************************************************************************/

using namespace std;

mutex mtx;
condition_variable cv;
queue<size_t> emptySlot;

// Inputs. ------------------------------------

const string homeDir=GetHomeDir();

const string synDataDir=homeDir+"/PROJ/t039.Lamella";
const size_t beginIndex=1279, endIndex=1279, TraceCnt=451; // [beginIndex, endIndex] inclusive.

const size_t nThread=5;
const bool reCreateTable=false;

const bool makePlots=true;
const double plotIndex=1279;

const double filterCornerLow=0.033, filterCornerHigh=0.3;
const double cutDeconSourceT1=-60, cutDeconSourceT2=60;
const double cutBeforeDeconT1=-120, cutBeforeDeconT2=120;
const double cutDeconResultT1=-50, cutDeconResultT2=50;
const double waterLevel=0.1, sigma=1.27398, deconFilterCornerLow=0.03, dt=0.025;
const string premDataDir=homeDir+"/PROJ/t039.PREM/201500000000";


// Outputs. ------------------------------------

const string dirPrefix=homeDir+"/PROJ/t041.REFL_Lamella";
const string outputDB="REFL_Lamella", outputTable="Decon";


// --------------------------------------------

void processThis(const size_t Index, int mySlot, const EvenSampledSignal &sESW);

int main(){

    // Update table.
    if (reCreateTable) {
        MariaDB::Query("drop table if exists "+outputDB+"."+outputTable);
        MariaDB::Query("create table "+outputDB+"."+outputTable+" (PairName varchar(30) not null unique primary key, EQ varchar(20), Gcarc double, FRSAmp double, FRSTime double, ScSFRSed varchar(200), ScSDeconed varchar(200), dirPrefix varchar(200))");
    }


    // Make ESW one time.

    /***********************************
    * 1. Read in PREM S waveform.      *
    ***********************************/

    SACSignals sESWData;

    sESWData=SACSignals (ShellExecVec("ls "+premDataDir+"/*.THT.sac"));

    if (sESWData.Size() != TraceCnt) {
        throw runtime_error("Reading error: PREM.");
    }


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
    * 2. Make S ESW stack from PREM.
    *
    *******************************************/

    sESWData.CheckAndCutToWindow(cutDeconSourceT1-10,cutDeconSourceT2+10);

    auto sESW=sESWData.XCorrStack(0,-15,15,2).second.first;
    sESW.FindPeakAround(0,10);
    sESW.ShiftTimeReferenceToPeak();
    sESW.CheckAndCutToWindow(cutDeconSourceT1, cutDeconSourceT2);
    sESW.FlipPeakUp();
    sESW.NormalizeToPeak();
    sESW.HannTaper(20);


    // Optional: make S ESW again, this time, stretch/shrink each S to match S ESW, then stack.

    sESWData.StretchToFit(sESW,-13,13,-0.3,0.3,0.25,true); // stretch, then compare waveform Amp_WinDiff.

    sESW=sESWData.XCorrStack(0,-15,15,2).second.first;
    sESW.FindPeakAround(0,10);
    sESW.ShiftTimeReferenceToPeak();
    sESW.CheckAndCutToWindow(cutDeconSourceT1, cutDeconSourceT2);
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
    * 2. Make ScS ESW, modified input S ESW.         *
    * 3. Preprocess it (Band-pass, strip S and etc.) *
    * 4. Decon.                                      *
    * 5. FRS.                                        *
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

    if(Data.Size() != TraceCnt) {
        throw runtime_error("Reading error: "+modelName);
    }


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


    /********************************************
    * 2. Make decon source:                     *
    * stack ScS                                 *
    * stretch/shrink S esw to look like ScS esw *
    ********************************************/

    SACSignals tmpData;

    tmpData=Data;
    tmpData.FindPeakAround(tmpData.GetTravelTimes("ScS"),10);
    tmpData.ShiftTimeReferenceToPeak();
    tmpData.FlipPeakUp();
    tmpData.NormalizeToPeak();

    tmpData.CheckAndCutToWindow(cutDeconSourceT1-10,cutDeconSourceT2+10);
    auto scsESW=tmpData.XCorrStack(0,-15,15,2).second.first;

    scsESW.FindPeakAround(0,10);
    scsESW.ShiftTimeReferenceToPeak();
    scsESW.CheckAndCutToWindow(cutDeconSourceT1,cutDeconSourceT2);
    scsESW.FlipPeakUp();
    scsESW.NormalizeToPeak();
    scsESW.HannTaper(10);


    EvenSampledSignal deconSource=sESW.StretchToFit(scsESW,-13,13,-0.1,0.2,0.25);

    deconSource.FindPeakAround(0,4);
    deconSource.ShiftTimeReferenceToPeak();
    deconSource.CheckAndCutToWindow(cutDeconSourceT1,cutDeconSourceT2);
    deconSource.NormalizeToPeak();
    deconSource.FlipPeakUp();
    deconSource.HannTaper(10);


    /****************************************************
    *                                                 
    * 3. Strip S waveform using the modified S ESW.    
    *
    *    Whether (and how to) fit the modeified S ESW to each
    *    S waveform is a choice.
    *                                                
    ****************************************************/

    SACSignals beforeStrip, afterStrip;

    if (makePlots && beginIndex+Index==plotIndex) {
        beforeStrip=Data;
    }

    vector<EvenSampledSignal> stripSources;
    for (size_t i=0; i<Data.Size(); ++i) {
        // stripSources.push_back(sESW);
        stripSources.push_back(sESW.StretchToFitHalfWidth(Data.GetData()[i])); // stretch to fit the half-width.
    }

    auto XCTimeShift=Data.CrossCorrelation(-15,15,stripSources,-15,15).first;
    Data.StripSignal(stripSources,XCTimeShift);

    // For plotting, need to explicitely shift it.
    for (size_t i=0; i<stripSources.size(); ++i) {
        stripSources[i].ShiftTime(XCTimeShift[i]);
    }

    if (makePlots && beginIndex+Index == plotIndex) {
        afterStrip=Data;
    }


    /***********
    * 4. Decon *
    ***********/

    // Find ScS from the S-stripped waveform.
    Data.FindPeakAround(Data.GetTravelTimes("ScS"),10);
    Data.ShiftTimeReferenceToPeak();
    Data.FlipPeakUp();
    Data.NormalizeToPeak();
    Data.CheckAndCutToWindow(cutBeforeDeconT1,cutBeforeDeconT2);

    SACSignals afterDecon, beforeDecon;

    if (makePlots && beginIndex+Index==plotIndex) {
        beforeDecon=Data;
    }

    // Decon.
    fftw_make_planner_thread_safe();
    Data.WaterLevelDecon(deconSource,waterLevel);
    Data.GaussianBlur(sigma);
    Data.Butterworth(deconFilterCornerLow,10000);

    Data.FindPeakAround(0,10);
    Data.ShiftTimeReferenceToPeak();
    Data.FlipPeakUp();
    Data.NormalizeToPeak();
    Data.CheckAndCutToWindow(cutDeconResultT1,cutDeconResultT2);


    // Output Decon waveforms.
    ShellExec("mkdir -p "+dirPrefix+"/Decon/"+modelName);
    Data.DumpWaveforms(dirPrefix+"/Decon/"+modelName,"StationName","","","trace");

    if (makePlots && beginIndex+Index==plotIndex) {
        afterDecon=Data;
    }


    /*********
    * 5. FRS *
    *********/

    // FRS.
    Data.FlipReverseSum(0);
    Data.CheckAndCutToWindow(0,15-dt*0.8);

    // Measure the peak time and amplitude.
    Data.FindPeakAround(7.5,7.5);

    // Output FRS waveforms.
    ShellExec("mkdir -p "+dirPrefix+"/DeconFRS/"+modelName);
    Data.DumpWaveforms(dirPrefix+"/DeconFRS/"+modelName,"StationName","","","frs");

    // Make plots.
    if (makePlots && beginIndex+Index==plotIndex) {

        string outfile;
        vector<GMT::Text> texts;

        for (size_t i=0; i<afterDecon.Size(); ++i) {

            if (i%17==0) { // A New page.
                if (outfile.empty()) outfile=GMT::BeginEasyPlot(69.5,40);
                else GMT::NewPage(outfile);
                GMT::MoveReferencePoint(outfile,"-Xf1i -Yf37.2i");
            }
            else {
                GMT::MoveReferencePoint(outfile,"-Y-2.3i");
            }

        

            // Plot to verify sESW.
            GMT::psbasemap(outfile, "-JX13i/2i -R-100/100/-1/1 -Bxa10 -Bya0.5 -BWSne -O -K -Xf1i");
            GMT::psxy(outfile, vector<double> {-100,100},vector<double> {0,0},"-J -R -W1p,gray,- -O -K");
            GMT::psxy(outfile, sESW,"-J -R -W1p,yellow -O -K");
            texts.clear();
            texts.push_back(GMT::Text(-50,0.9,"Un-modified S ESW in yellow",12,"LT"));
            GMT::pstext(outfile, texts,"-J -R -N -O -K");

            // Plot to verify scsESW.
            GMT::psbasemap(outfile,"-JX13i/2i -R-100/100/-1/1 -Bxa10 -Bya0.5 -BWSne -O -K -Xf14.5i");
            GMT::psxy(outfile,vector<double> {-100,100},vector<double> {0,0},"-J -R -W1p,gray,- -O -K");
            GMT::psxy(outfile,scsESW,"-J -R -W1p,black -O -K");
            texts.clear();
            texts.push_back(GMT::Text(-50,0.9,"Un-modified ScS ESW in black",12,"LT"));
            GMT::pstext(outfile, texts,"-J -R -N -O -K");

            // Plot to verify the decon source.
            GMT::psbasemap(outfile,"-JX13i/2i -R-100/100/-1/1 -Bxa10 -Bya0.5 -BWSne -O -K -Xf28i");
            GMT::psxy(outfile,vector<double> {-100,100},vector<double> {0,0},"-J -R -W1p,gray,- -O -K");
            GMT::psxy(outfile,scsESW,"-J -R -W1p,black -O -K");
            GMT::psxy(outfile,sESW,"-J -R -W1p,yellow -O -K");
            GMT::psxy(outfile,deconSource,"-J -R -W1p,red -O -K");
            texts.clear();
            texts.push_back(GMT::Text(-50,0.9,"Un-modified S ESW in yellow",12,"LT"));
            texts.push_back(GMT::Text(-50,0.75,"Un-modified ScS ESW in black",12,"LT"));
            texts.push_back(GMT::Text(-50,0.6,"Modified S ESW (Decon source) in red",12,"LT"));
            GMT::pstext(outfile, texts,"-J -R -N -O -K");

            // Plot to verify the stripping.
            GMT::psbasemap(outfile,"-JX13i/2i -R-100/100/-1/1 -Bxa10 -Bya0.5 -BWSne -O -K -Xf41.5i");
            GMT::psxy(outfile,beforeStrip,i,"-J -R -W1p,black -O -K");
            GMT::psxy(outfile,stripSources[i],"-J -R -W1p,yellow -O -K");
            GMT::psxy(outfile,afterStrip,i,"-J -R -W1p,green -O -K");
            texts.clear();
            texts.push_back(GMT::Text(-30,0.9,Float2String(afterStrip.GetMData()[i].gcarc,2),12,"LT"));
            texts.push_back(GMT::Text(-70,0.9,"After S-Subtraction in green",12,"LT"));
            texts.push_back(GMT::Text(-70,0.75,"Modified S ESW (Subtract source) in yellow",12,"LT"));
            texts.push_back(GMT::Text(-70,0.6,"Before S-Subtraction in black",12,"LT"));
            GMT::pstext(outfile,texts,"-J -R -N -O -K");

            // Plot to verify deconed trace.
            GMT::psbasemap(outfile,"-JX13i/2i -R-100/100/-1/1 -Bxa10 -Bya0.5 -BWSne -O -K -Xf55i");
            GMT::psxy(outfile,beforeDecon,i,"-J -R -W1p,black -O -K");
            GMT::psxy(outfile,deconSource,"-J -R -W1p,red -O -K");
            GMT::psxy(outfile,afterDecon,i,"-J -R -W1p,purple -O -K");
            texts.clear();
            texts.push_back(GMT::Text(-30,0.9,Float2String(afterDecon.GetMData()[i].gcarc,2),12,"LT"));
            texts.push_back(GMT::Text(-70,0.9,"After Decon in purple",12,"LT"));
            texts.push_back(GMT::Text(-70,0.75,"Modified S ESW (Decon source) in red",12,"LT"));
            texts.push_back(GMT::Text(-70,0.6,"Before Decon in black",12,"LT"));
            GMT::pstext(outfile,texts,"-J -R -N -O -K");

        }

        GMT::SealPlot(outfile);
        GMT::ps2pdf(outfile,__FILE__);
    }

    // Finished. Update database and notify the main thread to continue;
    lck.lock();

    auto stationNames=Data.GetStationNames();
    auto gcarcs=Data.GetDistances();
    auto peakAmps=Data.PeakAmp();
    auto peakTimes=Data.PeakTime();
    vector<vector<string>> sqlData(8,vector<string> ());
    for (size_t i=0; i<Data.Size(); ++i) {
        sqlData[0].push_back(modelName);
        sqlData[1].push_back(modelName+"_"+stationNames[i]);
        sqlData[2].push_back(Float2String(gcarcs[i],2));
        sqlData[3].push_back(Float2String(peakAmps[i],2));
        sqlData[4].push_back(Float2String(peakTimes[i],2));
        sqlData[5].push_back("DeconFRS/"+modelName+"/"+stationNames[i]+".frs");
        sqlData[6].push_back("Decon/"+modelName+"/"+stationNames[i]+".trace");
        sqlData[7].push_back(dirPrefix);
    }
    MariaDB::LoadData(outputDB,outputTable,vector<string> {"eq", "pairname", "gcarc", "FRSAmp", "FRSTime", "ScSFRSed", "ScSDeconed", "dirPrefix"},sqlData);

    emptySlot.push(mySlot);
    cv.notify_one();

    return;
}
