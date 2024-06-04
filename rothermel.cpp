#include <iostream>
#include <cmath>

struct FireSpreadResult {
    double rate_of_spread;
    double reaction_intensity;
    double fireline_intensity;
};

FireSpreadResult GetSimpleFireSpread(double fuelload, double fueldepth, double windspeed, double slope, double fuelmoisture, double fuelsav) {
    FireSpreadResult result = {0.0, 0.0, 0.0};

    if (fuelload > 0) {
        double wo = fuelload; // Ovendry fuel loading in (lb/ft^2)
        double fd = fueldepth; // Fuel depth (ft)
        double wv = windspeed * 88; // Wind velocity at midflame height (ft/minute)
        double fpsa = fuelsav; // Fuel Particle surface area to volume ratio (1/ft)
        double mf = fuelmoisture; // Fuel particle moisture content
        double h = 8000; // Fuel particle low heat content
        double pp = 32.0; // Ovendry particle density
        double st = 0.0555; // Fuel particle mineral content
        double se = 0.010; // Fuel Particle effective mineral content
        double mois_ext = 0.12; // Moisture content of extinction or 0.3 if dead

        double slope_rad = atan(slope);
        double tan_slope = tan(slope_rad);

        // Betas Packing ratio
        double Beta_op = 3.348 * pow(fpsa, -0.8189); // Optimum packing ratio
        double ODBD = wo / fd; // Ovendry bulk density
        double Beta = ODBD / pp; // Packing ratio
        double Beta_rel = Beta / Beta_op;

        // Reaction Intensity
        double WN = wo / (1 + st); // Net fuel loading
        double A = 133.0 / pow(fpsa, 0.7913); // Updated A
        double T_max = pow(fpsa, 1.5) * pow(495.0 + 0.0594 * pow(fpsa, 1.5), -1.0); // Maximum reaction velocity
        double T = T_max * pow((Beta / Beta_op), A) * exp(A * (1 - Beta / Beta_op)); // Optimum reaction velocity

        // Moisture damping coefficient
        double NM = 1.0 - 2.59 * (mf / mois_ext) + 5.11 * pow(mf / mois_ext, 2.0) - 3.52 * pow(mf / mois_ext, 3.0);

        // Mineral damping coefficient
        double NS = 0.174 * pow(se, -0.19);

        double RI = T * WN * h * NM * NS;

        // Propagating flux ratio
        double PFR = pow(192.0 + 0.2595 * fpsa, -1) * exp((0.792 + 0.681 * sqrt(fpsa)) * (Beta + 0.1));

        // Wind Coefficient
        double B = 0.02526 * pow(fpsa, 0.54);
        double C = 7.47 * exp(-0.1333 * pow(fpsa, 0.55));
        double E = 0.715 * exp(-3.59 * 10e-4 * fpsa);

        if (wv > (0.9 * RI)) {
            wv = 0.9 * RI;
        }

        double WC = (C * pow(wv, B)) * pow((Beta / Beta_op), (-E));

        // Slope coefficient
        double SC = 5.275 * pow(Beta, -0.3) * pow(tan_slope, 2);

        // Effective Heating Number
        double EHN = exp(-138.0 / fpsa);

        // Heat of preignition
        double QIG = 250.0 + 1116.0 * mf;

        // Rate of spread (ft per minute)
        double numerator = (RI * PFR * (1 + WC + SC));
        double denominator = (ODBD * EHN * QIG);
        double R = numerator / denominator;

        double RT = 384.0 / fpsa;
        double HA = RI * RT;

        // Fireline intensity
        double FI = (384.0 / fpsa) * RI * R;

        if (RI <= 0) {
            return result;
        }

        result.rate_of_spread = R;
        result.reaction_intensity = RI;
        result.fireline_intensity = FI;
    }

    return result;
}

int main() {
    double fuelload = 0.033976; // lbs/ft^2
    double fueldepth = 1; // feet
    double windspeed = 5; // mph
    double slope = 10; // degrees
    double fuelmoisture = 0.05; // as a proportion
    double fuelsav = 3500; // 1/ft

    FireSpreadResult result = GetSimpleFireSpread(fuelload, fueldepth, windspeed, slope, fuelmoisture, fuelsav);

    std::cout << "Rate of spread (ft/min): " << result.rate_of_spread << std::endl;
    std::cout << "Reaction Intensity (Btu/ft/min): " << result.reaction_intensity << std::endl;
    std::cout << "Fireline Intensity (Btu/ft): " << result.fireline_intensity << std::endl;

    return 0;
}
