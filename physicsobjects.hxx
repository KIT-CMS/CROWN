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
        auto df1 = df.Define(maskname, basefunctions::FilterPt(ptThreshold), {quantity});
        return df1;
        }

    auto CutEta(auto df, const std::string quantity, const std::string maskname, const float EtaThreshold){
        auto df1 = df.Define(maskname, basefunctions::FilterEta(EtaThreshold), {quantity});
        return df1;
        }

    auto CutMax(auto df, const std::string quantity, const std::string maskname, const float Threshold){
        auto df1 = df.Define(maskname, basefunctions::FilterMax(Threshold), {quantity});
        return df1;
        }

    namespace muon {

    } // ent namespace muon

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


    } // ent namespace tau

    namespace electron {

    } // ent namespace electron

    namespace jet {

    } // ent namespace jet


} // end namespace objects