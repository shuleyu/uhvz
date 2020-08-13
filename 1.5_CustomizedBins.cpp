#include<iostream>
#include<vector>
#include<map>
#include<set>

#include<PointInPolygon.hpp>
#include<MariaDB.hpp>
#include<GcpDistance.hpp>
#include<GMT.hpp>
#include<ShellExec.hpp>

using namespace std;

int main(){

    // Inputs. ------------------------

    // This file will create bins using the binning result or modeling result.

    const string DB="gen2CA_D",inputTable1="gen2CA_D.Master_a14",inputTable2="CenterDists_46";
    const string outputTable1="Bins",outputTable2="CenterDists";

    // bool: use gaussian cap weighting, or not.
    const vector<pair<bool,vector<size_t>>> newBinsFromOldBins {
        {true,{5,10,11,12,13,16,17,18,19,20,23,24,25,26,30}},
        {false,{5,10,11,12,13,16,17,18,19,20,23,24,25,26,30}},
        {true,{44,41,46}},
        {false,{44,41,46}},
        {true,{27,28}},
        {false,{27,28}},
        {false,{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,30,31,34,35,36,39,40,41,43,44,45,46}},
        {false,{27,28,29,32,33,37,38,42}},
    };
    const vector<pair<bool,string>> newBinsFromPolygonFiles {
        {true,"1.5_Bin1.txt"},
        {false,"1.5_Bin1.txt"},
        {true,"1.5_Bin2.txt"},
        {false,"1.5_Bin2.txt"},
        {true,"1.5_Bin3.txt"},
        {false,"1.5_Bin3.txt"}
    };

    // --------------------------------

    // get the scs hit locations.
    auto dataInfo=MariaDB::Select("pairname,hitlo,hitla from "+inputTable1+" order by pairname");
    const auto &pairNames=dataInfo.GetString("pairname");


    // get new bins.
    vector<pair<bool,set<string>>> newBins;

    for (auto item: newBinsFromOldBins) {
        set<string> selectedPairNames;
        for (size_t binN: item.second) {
            auto binData=MariaDB::Select("pairname from "+DB+"."+inputTable2+" where dist_"+to_string(binN)+">-0.5");
            for (auto pn: binData.GetString("pairname")) selectedPairNames.insert(pn);
        }
        newBins.push_back({item.first,selectedPairNames});
    }

    for (auto item: newBinsFromPolygonFiles) {
        set<string> selectedPairNames;
        vector<pair<double,double>> polygon;
        double x,y;
        ifstream fpin(item.second);
        while (fpin >> x >> y) polygon.push_back({x,y});
        fpin.close();
        for (size_t i=0; i<dataInfo.NRow(); ++i)
            if (PointInPolygon(polygon,make_pair(dataInfo.GetDouble("hitlo")[i], dataInfo.GetDouble("hitla")[i]),1))
                selectedPairNames.insert(pairNames[i]);
        newBins.push_back({item.first,selectedPairNames});
    }


    // Make outputs.
    vector<vector<string>> output1Data(7),output2Data{pairNames};
    vector<string> newColumns{"pairname"},outfiles;
    size_t binN=1;

    for (auto item: newBins) {
        auto useWeight=item.first;

        map<string,pair<double,double>> selectedData;
        for (auto pn: item.second) {
            auto it=lower_bound(pairNames.begin(),pairNames.end(),pn);
            auto index=distance(pairNames.begin(),it);
            selectedData[pn]=make_pair(dataInfo.GetDouble("hitlo")[index], dataInfo.GetDouble("hitla")[index]);
        }

        // get averaged location.
        double binLon=0,binLat=0;
        for (auto item2: selectedData) {
            binLon+=item2.second.first;
            binLat+=item2.second.second;
        }
        binLon/=selectedData.size();
        binLat/=selectedData.size();

        // get distances to center.
        double maxDist=0;
        output2Data.push_back({});
        for (auto pn: output2Data[0]) {
            if (selectedData.find(pn)==selectedData.end()) output2Data.back().push_back("-1");
            else {
                auto item2=selectedData[pn];
                double dist=GcpDistance(item2.first, item2.second, binLon, binLat);
                maxDist=max(maxDist,dist);

                if (!useWeight) output2Data.back().push_back("0");
                else output2Data.back().push_back(to_string(dist));
            }
        }

        newColumns.push_back("dist_"+to_string(binN));
        output1Data[0].push_back(to_string(binN++));
        output1Data[1].push_back(to_string(item.second.size()));
        output1Data[2].push_back(to_string(maxDist));
        output1Data[3].push_back(to_string(binLon));
        output1Data[4].push_back(to_string(binLat));
        output1Data[5].push_back(to_string(binLon));
        output1Data[6].push_back(to_string(binLat));


        // plot.
        string outfile=GMT::BeginEasyPlot(15,15,ShellExec("pwd",true)+"/"+string(__FILE__));
        outfiles.push_back(outfile);

        double scale=13.0/37;
        GMT::MoveReferencePoint(outfile,"-Xf1i -Yf1i");
        GMT::psbasemap(outfile,"-Jx"+to_string(scale)+"id/"+to_string(scale)+"id -R-101/-64/-5/27 -Bxa10f5 -Bya10f5 -Bx+lLongitude -By+lLatitude -BWSne -O -K");
        GMT::pscoast(outfile,"-J -R -W0.3p,black -O -K");

        // plot hit locations.
        vector<pair<double,double>> hitLocations;
        for (auto item: selectedData) 
            hitLocations.push_back({item.second.first, item.second.second});
        GMT::psxy(outfile,hitLocations,"-J -R -Sc0.05i -Ggray -O -K");

        // plot texts.
        vector<GMT::Text> texts;
        texts.push_back(GMT::Text(binLon,binLat,to_string(selectedData.size()),12,"CM"));
        GMT::pstext(outfile,texts,"-J -R -N -O -K");

        GMT::SealPlot(outfile);
    }

    // Output to database.
    MariaDB::Query("create database if not exists "+DB);
    MariaDB::Query("drop table if exists "+DB+"."+outputTable1);
    MariaDB::Query("drop table if exists "+DB+"."+outputTable2);


    MariaDB::Query("create table "+DB+"."+outputTable1+" (bin integer, nRecord integer, radius double, lon_before double, lat_before double, lon double, lat double)");
    string ss="";
    for (size_t i=1; i<newColumns.size(); ++i) ss+=", "+newColumns[i]+" double";
    MariaDB::Query("create table "+DB+"."+outputTable2+" (pairname varchar(30) not null unique primary key"+ss+")");


    MariaDB::LoadData(DB,outputTable1,vector<string> {"bin","nRecord","radius","lon_before","lat_before","lon","lat"},output1Data);
    MariaDB::LoadData(DB,outputTable2,newColumns,output2Data);


    // Make PDF.
    
    remove("tmp.ps");
    for (auto file: outfiles) {
        ShellExec("cat "+file+" >> tmp.ps");
        remove(file.c_str());
    }
    string pdffile=__FILE__;
    ShellExec("ps2pdf tmp.ps "+pdffile.substr(0,pdffile.find_last_of("."))+".pdf");
    remove("tmp.ps");

    return 0;
}
