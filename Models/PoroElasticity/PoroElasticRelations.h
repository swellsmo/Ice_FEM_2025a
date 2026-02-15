#ifndef POROELASTICRELATIONS_H
#define POROELASTICRELATIONS_H

#include <vector>
#include <iostream>
#include <unordered_map>

#include "../BaseModel/BaseModel.h"

class PoroElasticRelations{
    public:
        PoroElasticRelations(inputData& inputs);
        ~PoroElasticRelations();

        double GetSaturation(double p, double& dSw, double& ddSw);
        double GetRelativePermeability(double p, double por, double& dk_dp, double& dk_dpor);
        double GetBiotCoefficient(double poros, double& dAlpha);

        enum BiotModel{PorosityBasedBiot};
        std::unordered_map<std::string, BiotModel> BiotModels = {
            {"Porosity", PorosityBasedBiot}
        };       

        enum SaturationModel{ExponentSaturation, ExponentSmoothSaturation, TanhSaturation, FullSaturation, LinearSaturation};
        std::unordered_map<std::string, SaturationModel> SaturationModels = {
            {"Exponent", ExponentSaturation},
            {"ExponentSmooth", ExponentSmoothSaturation},
            {"Tanh", TanhSaturation},
            {"Saturated", FullSaturation},
            {"Linear", LinearSaturation}
        };

        SaturationModel satModel;
        double sat_exp_pref, sat_exp_diffexp;
        double relperm_min, S_0;

        BiotModel BiotMethod;
    protected:

    private:

};
#endif

