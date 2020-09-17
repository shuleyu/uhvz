#include<iostream>
#include<vector>

#include<MariaDB.hpp>
#include<GcpDistance.hpp>
#include<MeshGrid.hpp>
#include<GMT.hpp>
#include<WayPoint.hpp>

using namespace std;

/*

This code bin the sampling points for a given bin radius and a given bin center increment.

*/


// Inputs. ------------------------

const string inputTable = "gen2CA_D.Master_a14";
const double lonMin = -99, lonMax = -64, latMin = -5, latMax = 27;
const double binRadius = 4, binInc = binRadius;
const size_t recordCntThreshold = 100;
const bool makePlot = true, recreateTable = false;

// Outputs. --------------------------------

const string outputDB="gen2CA_D";
const string outputTable1="Bins",outputTable2="CenterDists";

// -----------------------------------------

int main(){
    
    auto dataInfo=MariaDB::Select("pairname, hitlo, hitla, shift_gcarc from " + inputTable);

    // Make bins.
    vector<vector<double>> p{{latMax,latMin,binInc},{lonMin,lonMax+1e-5,binInc}};
    vector<vector<string>> allData{dataInfo.GetString("pairname")};
    auto grid=MeshGrid(p,1);
    vector<string> nRecord,lon_before,lat_before,lon,lat,bin,radius,averagedDist, columnNames{"pairname"};
    size_t binN=1;
    string newColumns="";

    for (auto item:grid) {

        size_t recordCnt = 0;
        double newLon = 0, newLat = 0, dist = 0;
        vector<double> curData(dataInfo.NRow(),0);
        for (size_t i=0; i<dataInfo.NRow(); ++i) {
            curData[i]=GcpDistance(dataInfo.GetDouble("hitlo")[i], dataInfo.GetDouble("hitla")[i], item[1], item[0]);
            if (curData[i]<=binRadius) {
                ++recordCnt;
                dist += dataInfo.GetDouble("shift_gcarc")[i];
                newLon += dataInfo.GetDouble("hitlo")[i];
                newLat += dataInfo.GetDouble("hitla")[i];
            }
            else curData[i]=-1;
        }

        if (recordCnt>=recordCntThreshold) {
            newLon/=recordCnt;
            newLat/=recordCnt;
            dist /= recordCnt;

            bin.push_back(to_string(binN++));
            radius.push_back(to_string(binRadius));
            lon_before.push_back(to_string(item[1]));
            lat_before.push_back(to_string(item[0]));
            lon.push_back(to_string(newLon));
            lat.push_back(to_string(newLat));
            nRecord.push_back(to_string(recordCnt));
            averagedDist.push_back(to_string(dist));

            columnNames.push_back("dist_"+bin.back());
            newColumns+=", "+columnNames.back()+" double ";
            allData.push_back({});
            for (size_t i=0; i<dataInfo.NRow(); ++i) {
                if (curData[i]<0) allData.back().push_back("-1");
                else allData.back().push_back(to_string(GcpDistance(dataInfo.GetDouble("hitlo")[i], dataInfo.GetDouble("hitla")[i], newLon, newLat)));
            }
        }
    }

    // Output.

    if (recreateTable) {

        MariaDB::Query("create database if not exists "+outputDB);
        MariaDB::Query("drop table if exists "+outputDB+"."+outputTable1);
        MariaDB::Query("drop table if exists "+outputDB+"."+outputTable2);

        MariaDB::Query("create table "+outputDB+"."+outputTable1+" (bin integer, nRecord integer, radius double, lon_before double, lat_before double, lon double, lat double, averagedDist double)");
        MariaDB::Query("create table "+outputDB+"."+outputTable2+" (pairname varchar(30) not null unique primary key"+newColumns+")");

        MariaDB::LoadData(outputDB,outputTable1,vector<string> {"bin","nRecord","radius","lon_before","lat_before","lon","lat", "averagedDist"},vector<vector<string>> {bin,nRecord,radius,lon_before,lat_before,lon,lat, averagedDist});
        MariaDB::LoadData(outputDB,outputTable2,columnNames,allData);
    }


    // plot.
    if (makePlot) {

        string outfile=GMT::BeginEasyPlot(15,15,ShellExec("pwd",true)+"/"+string(__FILE__));
        double scale=13.0/37;
        GMT::MoveReferencePoint(outfile,"-Xf1i -Yf1i");
        GMT::psbasemap(outfile,"-Jx"+to_string(scale)+"id/"+to_string(scale)+"id -R-101/-64/-5/27 -Bxa10f5 -Bya10f5 -Bx+lLongitude -By+lLatitude -BWSne -O -K");
        GMT::pscoast(outfile,"-J -R -W0.3p,black -A10000 -O -K");

        // plot bins.
        auto binInfo=MariaDB::Select("bin, lon, lat, lon_before, lat_before, nRecord, radius, averagedDist from "+outputDB+"."+outputTable1);

        //GMT::psxy(outfile, dataInfo.GetDouble("hitlo"), dataInfo.GetDouble("hitla"), "-J -R -Sc0.02i -Gblue -O -K");

        for (size_t i=0; i<binInfo.NRow(); ++i){
            vector<pair<double,double>> circle;
            for (double az=0; az<360; az+=0.05) {
                circle.push_back(WayPoint(binInfo.GetDouble("lon")[i], binInfo.GetDouble("lat")[i], az, binInfo.GetDouble("radius")[i]));
            }
            int grayScale = 65 + (80.0 - binInfo.GetDouble("averagedDist")[i]) * (255 - 65) / 80;
            // GMT::psxy(outfile,circle,"-J -R -W0.75p,black -O -K");
            vector<GMT::Text> texts;
            texts.push_back(GMT::Text(binInfo.GetDouble("lon")[i],binInfo.GetDouble("lat")[i],to_string(binInfo.GetInt("nRecord")[i]),12,"CM"));
            // GMT::pstext(outfile,texts,"-J -R -N -O -K");

            GMT::psxy(outfile, binInfo.GetDouble("lon")[i], binInfo.GetDouble("lat")[i], "-J -R -Sc0.5i -G" + to_string(grayScale) + " -O -K");
            GMT::psxy(outfile, binInfo.GetDouble("lon")[i], binInfo.GetDouble("lat")[i], "-J -R -S+0.15i -W2p,green -O -K");
        }

        // Make PDF.
        GMT::SealPlot(outfile);
        GMT::ps2pdf(outfile,__FILE__);

    }


    return 0;
}
