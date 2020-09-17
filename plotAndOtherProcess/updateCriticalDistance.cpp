#include<iostream>
#include<fstream>
#include<map>

#include<MariaDB.hpp>
#include<Float2String.hpp>

using namespace std;

const string premTable="REFL_PREM.Properties";
const string ulvzTable="REFL_ULVZ.Properties";
const string uhvzTable="REFL_UHVZ.Properties";

int main (){


    auto res1=MariaDB::Select("concat('PREM_',modelname) as modelname, thickness, dvs, drho from " + premTable);
    auto res2=MariaDB::Select("concat('ULVZ_',modelname) as modelname, thickness, dvs, drho from " + ulvzTable);
    auto res3=MariaDB::Select("concat('UHVZ_',modelname) as modelname, thickness, dvs, drho from " + uhvzTable);
    auto res=res1+res2+res3;

    cout << res.GetString("modelname").size() << endl;


    ifstream fpin("criticalDistance.txt");
    map<pair<double,double>,double> M;
    double t,v,d;
    string tmpstr;
    while (fpin >> tmpstr >> t >> v >> d){
        M[make_pair(t,v)]=d;
    }
    fpin.close();

    vector<vector<string>> sqlData(5,vector<string> ());
    for (size_t i=0; i<res.NRow(); ++i) {
        sqlData[0].push_back(res.GetString("modelname")[i]);
        sqlData[1].push_back(Float2String(res.GetDouble("thickness")[i],1));
        sqlData[2].push_back(Float2String(res.GetDouble("dvs")[i],1));
        sqlData[3].push_back(Float2String(res.GetDouble("drho")[i],1));

        string s=sqlData[0].back();
        if (s.find("UHVZ") != string::npos) {
            sqlData[4].push_back(Float2String(M[make_pair(res.GetDouble("thickness")[i], res.GetDouble("dvs")[i])],3));
        }
        else {
            sqlData[4].push_back("100");
        }
    }

    MariaDB::LoadData("gen2CA_D","Properties",vector<string> {"modelname","thickness","dvs","drho","criticalDist"},sqlData);

    return 0;
}
