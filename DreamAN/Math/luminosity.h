// luminosity.hpp
#pragma once
#include <stdexcept>

// Physical constants (SI unless noted)
inline constexpr double kElementaryCharge_C = 1.602176634e-19;     // e [C]
inline constexpr double kAvogadro_per_mol   = 6.02214076e23;       // N_A [1/mol]

// Bundle optional corrections
struct Corrections {
    double live_time   = 1.0;  // f_live   (0..1)
    double prescale    = 1.0;  // f_PS     (>=1; if you prefer 1/PS, put it here)
    double boiling     = 1.0;  // f_boil   (density correction factor)
    double good_frac   = 1.0;  // f_good   (fraction of good running)
};

// Areal number density [1/cm^2]
inline double areal_number_density_cm2(double density_g_cm3,
                                       double length_cm,
                                       double A_eff_g_per_mol = 1.0079)
{
    if (density_g_cm3 <= 0.0 || length_cm <= 0.0 || A_eff_g_per_mol <= 0.0)
        throw std::invalid_argument("Inputs must be positive.");
    return (density_g_cm3 * length_cm * kAvogadro_per_mol) / A_eff_g_per_mol;
}

// Integrated luminosity from delivered charge Q.
// Returns L_int in cm^-2.
inline double integrated_luminosity_cm2(double charge_C,
                                        double density_g_cm3,
                                        double length_cm,
                                        double A_eff_g_per_mol = 1.0079,
                                        Corrections corr = {})
{
    if (charge_C < 0.0) throw std::invalid_argument("Charge cannot be negative.");

    const double Nt = areal_number_density_cm2(density_g_cm3, length_cm, A_eff_g_per_mol);
    const double electrons = charge_C / kElementaryCharge_C;   // dimensionless
    const double f = corr.live_time * corr.prescale * corr.boiling * corr.good_frac;
    return electrons * Nt * f;  // cm^-2
}

// Instantaneous luminosity from beam current I.
// Returns L in cm^-2 s^-1.
inline double instantaneous_luminosity_cm2s(double current_A,
                                            double density_g_cm3,
                                            double length_cm,
                                            double A_eff_g_per_mol = 1.0079)
{
    if (current_A < 0.0) throw std::invalid_argument("Current cannot be negative.");

    const double Nt = areal_number_density_cm2(density_g_cm3, length_cm, A_eff_g_per_mol);
    const double electrons_per_s = current_A / kElementaryCharge_C; // s^-1
    return electrons_per_s * Nt; // cm^-2 s^-1
}

// Convenience converters
inline double cm2_to_inv_nb(double L_cm2) { return L_cm2 * 1e-33; } // nb^-1
inline double cm2_to_inv_pb(double L_cm2) { return L_cm2 * 1e-36; } // pb^-1


