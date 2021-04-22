#include "ROOT/RDataFrame.hxx"
#include "basefunctions.hxx"
/*
class Foo {
    Channel mChannel;

public:
    Foo(const Channel channel){
        mChannel = foo; 
    }

    float operator()(float x){
        if (mChannel == Channel::MT) {
            return x > 42;
        } else {
            ...;
        }
    }
};

ROOT::RDataFrame df(...);
Foo foo(...);
df.Define("y", foo, {"x"});

float MyFoo(float x, Channel channel) {
    if (channel == Channel::MT) {
        return 42;
    }
}
*/


namespace physicsobject {
    auto CutPt(auto df, const std::string quantity, const std::string maskname, const float ptThreshold){
        auto df1 = df.Define(maskname, basefunctions::FilterMin(ptThreshold), {quantity});
        return df1;
        }

    auto CutEta(auto df, const std::string quantity, const std::string maskname, const float EtaThreshold){
        auto df1 = df.Define(maskname, basefunctions::FilterAbsMax(EtaThreshold), {quantity});
        return df1;
        }

    auto CutDz(auto df, const std::string quantity, const std::string maskname, const float Threshold){
        auto df1 = df.Define(maskname, basefunctions::FilterMax(Threshold), {quantity});
        return df1;
        }

    auto CutDxy(auto df, const std::string quantity, const std::string maskname, const float Threshold){
        auto df1 = df.Define(maskname, basefunctions::FilterMax(Threshold), {quantity});
        return df1;
        }

    auto CombineMasks(auto df, const std::string maskname, std::vector<std::string> MaskList){
        // if(MaskList.size() == 0)
        // {
        //     std::cout << "No masks to combine for filtering\n";
        //     return df;
        // }
        // std::string former_mask = MaskList.pop();
        // std::string latter_mask = "";
        // while(MaskList.size() > 0){
        // }
        std::string filterstring;
        for (auto mask: MaskList){
            filterstring.append(mask + "*");
        }
        filterstring.pop_back(); // removing the last * from the string
        return df.Define(maskname, filterstring); // TODO make this a compiled version
    }

    auto FilterMasks(auto df, const std::string maskname){
        auto df1 = df.Filter([](const ROOT::RVec<Int_t>& mask){
            auto result = Any(mask);
            return result;
        },{maskname});
        return df1;
    }

    auto FilterObjects(auto df, const std::string objectcounter, const int minThreshold){
            return df.Filter([minThreshold](const UInt_t& nobject){return nobject>=minThreshold;}, {objectcounter});
        }

    namespace muon {

        auto FilterID(auto df, const std::string maskname, const std::string nameID){
            auto df1 = df.Define(maskname, [](const ROOT::RVec<Bool_t>& id){return id;}, {nameID});
            return df1;
        }

        auto FilterIsolation(auto df, const std::string maskname, const std::string isolationName, const float Threshold){
            auto df1 = df.Define(maskname, basefunctions::FilterMax(Threshold), {isolationName});
            return df1;
        }

    } // end namespace muon

    namespace tau {
        auto FilterDecayModes(auto df, const std::string maskname, const std::vector<int> SelectedDecayModes){
            auto df1 = df.Define(maskname,
                [SelectedDecayModes](const ROOT::RVec<Int_t>& decaymodes){
                    ROOT::RVec<int> mask;
                    for (auto n: decaymodes){
                        mask.push_back(int(std::find(SelectedDecayModes.begin(), SelectedDecayModes.end(), n) != SelectedDecayModes.end()));
                    }
                    return mask;},
                {"Tau_decayMode"});
            return df1;
        }

        auto FilterTauID(auto df, const std::string maskname, const std::string nameID, const int idxID){
            auto df1 = df.Define(maskname, basefunctions::FilterID(idxID), {nameID});
            return df1;
        }


    } // end namespace tau

    namespace electron {

    } // end namespace electron

    namespace jet {

    } // end namespace jet


} // end namespace objects