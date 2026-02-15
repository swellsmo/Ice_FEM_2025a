#include "PoroElasticRelations.h"

PoroElasticRelations::PoroElasticRelations(inputData& inputs){

    std::string SatModelName;
    inputs.GetRequired(SatModelName, {"properties", "Porosity", "Saturation","Model"});
    if (SaturationModels.count(SatModelName)){
        satModel = SaturationModels[SatModelName];
    } else {
        throw std::invalid_argument("Porosity model type "+SatModelName+" not defined,");
    }

    if (satModel == ExponentSaturation || satModel == ExponentSmoothSaturation || satModel == TanhSaturation){
        inputs.GetRequired(sat_exp_pref, {"properties", "Porosity", "Saturation","pref"});
        inputs.GetRequired(sat_exp_diffexp, {"properties", "Porosity", "Saturation","diffExp"});
    }

    if (satModel == LinearSaturation){
        inputs.GetRequired(sat_exp_pref, {"properties", "Porosity", "Saturation","pref"});
        inputs.GetRequired(sat_exp_diffexp, {"properties", "Porosity", "Saturation","diffExp"});
    }

    std::string BiotMethodName;
    inputs.GetRequired(BiotMethodName, {"properties", "Porosity", "bulk_t_model"});
    if (BiotModels.count(BiotMethodName)){
        BiotMethod = BiotModels[BiotMethodName];
    } else {
        throw std::invalid_argument("Biot porosity model type "+BiotMethodName+" not defined,");
    }

    inputs.GetRequired(relperm_min, {"properties", "Porosity","relperm_min"});
    inputs.GetRequired(S_0, {"properties", "Porosity","S_min"});
}

PoroElasticRelations::~PoroElasticRelations(){

}

double PoroElasticRelations::GetSaturation(double p, double& dSw, double& ddSw){
    double Sw;

    switch (satModel){
        case FullSaturation:{ 
            Sw = 1.0;
            dSw = 0.0;
            ddSw = 0.0;
        } break;
        case ExponentSmoothSaturation:{
            double Offset = 0.0;
            double x = p / sat_exp_pref + Offset;
            Sw   = S_0 + (1.0-S_0)*sigmoid(x);
            dSw  = (1.0-S_0)*(1.0 / sat_exp_pref) * sigmoid_derivative(x);
            ddSw =-(1.0-S_0)*(1.0 / (sat_exp_pref * sat_exp_pref)) * sigmoid_second_derivative(x);

        } break;
        case TanhSaturation:{ //note: functionally the same as ExponentSmooth, but this one is used in other literature
            Sw   = S_0 + (1.0-S_0)*(0.5+0.5*tanh(p/sat_exp_pref));
            dSw  = (1.0-S_0)*0.5/sat_exp_pref*std::pow(cosh(p/sat_exp_pref), -2);
            ddSw =-(1.0-S_0)/sat_exp_pref/sat_exp_pref*std::pow(cosh(p/sat_exp_pref), -3)*sinh(p/sat_exp_pref);
        } break;
        case ExponentSaturation:{
            if (p<0){
                Sw = S_0 + (1.0-S_0)*std::exp(p/sat_exp_pref);
                dSw = (1.0-S_0)*1.0/sat_exp_pref*std::exp(p/sat_exp_pref);
                ddSw = (1.0-S_0)*1.0/sat_exp_pref/sat_exp_pref*std::exp(p/sat_exp_pref);
            } else {
                Sw = 1.0;
                dSw = 0.0;
                ddSw = 0.0;
            }
        } break;
        case LinearSaturation:{
            if (p<-sat_exp_pref){
                Sw = S_0;
                dSw = 0.0;
                ddSw = 0.0;
            } else if (p<0){
                Sw = S_0 + (1.0-S_0)*(p+sat_exp_pref)/sat_exp_pref;
                dSw = (1.0-S_0)/sat_exp_pref;
                ddSw = 0.0;
            } else {
                Sw = 1.0;
                dSw = 0.0;
                ddSw = 0.0;
            }
        } break;
        default:{
            throw std::invalid_argument("Saturation function not defined in PoroElasticRelations.cpp,");
        }
    }

    if (std::isnan(Sw) || std::isnan(dSw) || std::isnan(ddSw)){
        std::stringstream  ErrString;
        ErrString << "Saturnation function is NaN in PoroElasticRelations.cpp,\n" << "Sw=" << Sw << ", dSw=" << dSw << ", ddSw=" << ddSw << ", pw=" << p << "\n";
        std::cout << ErrString.str();
        throw std::invalid_argument(ErrString.str());
    }

    return Sw;
}

double PoroElasticRelations::GetRelativePermeability(double p, double por, double& dk_dp, double& dk_dpor){
    if (por<0) por = 0.0;
    
    double krel, krel0;
    double Sw, dSw, ddSw;
    Sw = GetSaturation(p ,dSw, ddSw);
    switch (satModel){
        case FullSaturation:{ 
            krel0 = 1.0;
            dk_dp = 0.0;
        } break;
        case ExponentSaturation:{
            krel0 = std::pow(Sw-S_0, sat_exp_diffexp);
            dk_dp = sat_exp_diffexp*std::pow(Sw-S_0, sat_exp_diffexp-1.0)*dSw;
        } break;
        case ExponentSmoothSaturation:{
            krel0 = std::pow(Sw-S_0, sat_exp_diffexp);
            dk_dp = sat_exp_diffexp*std::pow(Sw-S_0, sat_exp_diffexp-1.0)*dSw;
        } break;
        case TanhSaturation:{
            krel0 = std::pow(Sw-S_0, sat_exp_diffexp);
            dk_dp = sat_exp_diffexp*std::pow(Sw-S_0, sat_exp_diffexp-1.0)*dSw;
        } break;
        case LinearSaturation:{
            krel0 = std::pow(Sw-S_0, sat_exp_diffexp);
            dk_dp = sat_exp_diffexp*std::pow(Sw-S_0, sat_exp_diffexp-1.0)*dSw;
        } break;
        default:{
            throw std::invalid_argument("Saturation function not defined in PoroElasticRelations.cpp,");
        }
    }
    if (std::isnan(krel0) || std::isnan(dk_dp)){
        std::stringstream  ErrString;
        ErrString << "Relative permeability function is NaN in PoroElasticRelations.cpp," << "Sw=" << Sw << "\n";
        throw std::invalid_argument(ErrString.str());
    }

    krel = krel0*por*(1.0-relperm_min)+relperm_min*(Sw-S_0);
    dk_dp = dk_dp*por*(1.0-relperm_min)+relperm_min*dSw;
    dk_dpor = krel0*(1.0-relperm_min);

    return krel;
}

double PoroElasticRelations::GetBiotCoefficient(double poros, double& dAlpha){
    double alpha;

    switch (BiotMethod){ //normal relation alpha=1-K_S/K_T   (skeleton over total)
        case PorosityBasedBiot:{ 
            alpha = 1.0-(1.0-std::max(0.0,poros));
            dAlpha = 1.0;
        } break;
        default:{
            throw std::invalid_argument("Saturation function not defined in PoroElasticRelations.cpp,");
        }
    }
    return alpha;
}

