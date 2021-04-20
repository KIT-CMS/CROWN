#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"

enum Channel {
    MT = 0,
    ET = 1,
    TT = 2,
    EM = 3
};


namespace basefunctions {

    auto FilterMax(float cut){
        return [cut](const ROOT::RVec<float>& values){
            ROOT::RVec<int> mask = values < cut;
            return mask;
        };
    }

    auto FilterPt(float ptCut){
        return [ptCut](const ROOT::RVec<float>& pts){
            ROOT::RVec<int> mask = pts > ptCut;
            return mask;
        };
    }

    auto FilterEta(float etaCut){
        return [etaCut](const ROOT::RVec<float>& etas){
            ROOT::RVec<int> mask = abs(etas) < etaCut;
            return mask;
        };
    }

    auto FilterID(int index){
        return [index](const ROOT::RVec<UChar_t>& IDs){
            ROOT::RVec<int> mask;
            for(auto const ID : IDs){
                mask.push_back(std::min(1, int(ID & 1 << index - 1)));
            }
            return mask;
        };
    }
} // namespace basefunctions