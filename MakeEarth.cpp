#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace self_oxidation {

enum class GrowthPhase {
  kPrecomputeAccretion,
  kGiantImpact,
  kLateVeneer,
};

enum class ImpactorClass {
  kECLikeUndifferentiated,
  kECLikeDifferentiated,
  kCILikeDifferentiated,
};

enum class CarbonGasSpecies {
  kCO2,
  kCO,
  kCH4,
};

enum class HydrogenGasSpecies {
  kH2O,
  kH2,
};

enum class NitrogenGasSpecies {
  kN2,
};

enum class InventorySource {
  kNone,
  kEcBulk,
  kCiBulk,
  kPrecomputeOutputAt0p10EarthMass,
  kPrecomputeCoreAt0p10EarthMass,
};

struct ElementPpm {
  double carbon_ppm = 0.0;
  double hydrogen_ppm = 0.0;
  double nitrogen_ppm = 0.0;
};

struct ElementMassKg {
  double carbon_kg = 0.0;
  double hydrogen_kg = 0.0;
  double nitrogen_kg = 0.0;
};

struct AlloyComposition {
  double x_s = 0.02;
  double x_c = 0.01;
  double x_ni = 0.05;
  double x_si = 0.06;

  [[nodiscard]] double x_fe() const {
    return 1.0 - x_s - x_c - x_ni - x_si;
  }

  void Validate() const {
    auto validate_fraction = [](double value, const char* label) {
      if (value < 0.0 || value > 1.0) {
        throw std::runtime_error(std::string(label) + " must be in [0, 1].");
      }
    };

    validate_fraction(x_s, "x_s");
    validate_fraction(x_c, "x_c");
    validate_fraction(x_ni, "x_ni");
    validate_fraction(x_si, "x_si");
    if (x_fe() < 0.0 || x_fe() > 1.0) {
      throw std::runtime_error(
          "Alloy composition is invalid: x_Fe became negative.");
    }
  }
};

struct RedoxState {
  double fe3_over_fe_total = 0.0;
  double x_feo_silicate = 0.056;
  double delta_iw_eq = -2.0;
  double delta_iw_surface = std::numeric_limits<double>::quiet_NaN();
  double oxidized_fe0_mass_kg = 0.0;
};

struct ReservoirState {
  double mass_kg = 0.0;
  ElementPpm volatiles_ppm{};
  double delta15n_permil = 0.0;
};

struct AtmosphereState {
  ElementMassKg elements_kg{};
  double delta15n_permil = 0.0;
  CarbonGasSpecies carbon_species = CarbonGasSpecies::kCO2;
  HydrogenGasSpecies hydrogen_species = HydrogenGasSpecies::kH2O;
  NitrogenGasSpecies nitrogen_species = NitrogenGasSpecies::kN2;
};

struct MeltPoolState {
  double target_radius_m = 0.0;
  double core_radius_m = 0.0;
  double mantle_thickness_m = 0.0;
  double combined_body_mass_kg = 0.0;
  double combined_body_radius_m = 0.0;
  double impactor_density_kg_m3 = 0.0;
  double impactor_diameter_m = 0.0;
  double impact_velocity_m_s = 0.0;
  double gravity_m_s2 = 0.0;
  double melt_energy_j_kg = 0.0;
  double target_melt_volume_m3 = 0.0;
  double target_melt_mass_kg = 0.0;
  double impactor_silicate_mass_kg = 0.0;
  double impactor_metal_mass_kg = 0.0;
  double equilibrating_silicate_mass_kg = 0.0;
  double raw_depth_m = 0.0;
  double depth_m = 0.0;
  double pressure_gpa = 0.0;
  double temperature_k = 0.0;
  int iteration_count = 0;
};

struct PlanetState {
  ReservoirState mantle{};
  ReservoirState core{};
  AtmosphereState atmosphere{};
  RedoxState redox{};

  [[nodiscard]] double total_mass_kg() const {
    return mantle.mass_kg + core.mass_kg;
  }

  void Validate() const {
    auto validate_non_negative = [](double value, const char* label) {
      if (value < 0.0) {
        throw std::runtime_error(std::string(label) + " must be non-negative.");
      }
    };

    validate_non_negative(mantle.mass_kg, "mantle.mass_kg");
    validate_non_negative(core.mass_kg, "core.mass_kg");
    validate_non_negative(atmosphere.elements_kg.carbon_kg,
                          "atmosphere.elements_kg.carbon_kg");
    validate_non_negative(atmosphere.elements_kg.hydrogen_kg,
                          "atmosphere.elements_kg.hydrogen_kg");
    validate_non_negative(atmosphere.elements_kg.nitrogen_kg,
                          "atmosphere.elements_kg.nitrogen_kg");
  }
};

struct Impactor {
  std::string label;
  GrowthPhase phase = GrowthPhase::kPrecomputeAccretion;
  ImpactorClass impactor_class = ImpactorClass::kECLikeUndifferentiated;
  double mass_kg = 0.0;
  double metallic_fe_fraction = 0.325;
  bool is_differentiated = false;
  InventorySource silicate_inventory_source = InventorySource::kEcBulk;
  InventorySource metal_inventory_source = InventorySource::kNone;
  ElementPpm silicate_volatiles_ppm{};
  ElementPpm metal_volatiles_ppm{};
  double delta15n_silicate_permil = 0.0;
  double delta15n_metal_permil = 0.0;
  AlloyComposition alloy{};

  void Validate() const {
    if (label.empty()) {
      throw std::runtime_error("Impactor label must not be empty.");
    }
    if (mass_kg <= 0.0) {
      throw std::runtime_error("Impactor mass must be positive.");
    }
    if (metallic_fe_fraction < 0.0 || metallic_fe_fraction > 1.0) {
      throw std::runtime_error("Impactor metallic_fe_fraction must be in [0, 1].");
    }
    alloy.Validate();
  }
};

struct StepInput {
  std::size_t index = 0;
  Impactor impactor{};
};

struct StepResult {
  std::size_t index = 0;
  Impactor resolved_impactor{};
  PlanetState planet_after{};
  MeltPoolState melt_pool{};
  double dn_metal_silicate = std::numeric_limits<double>::quiet_NaN();
  double dh_metal_silicate = std::numeric_limits<double>::quiet_NaN();
  double dc_metal_silicate = std::numeric_limits<double>::quiet_NaN();
  double xo_metal = std::numeric_limits<double>::quiet_NaN();
  double delta15n_silicate_permil = std::numeric_limits<double>::quiet_NaN();
  double delta15n_metal_permil = std::numeric_limits<double>::quiet_NaN();
  double erosion_fraction = 0.0;
  ElementMassKg eroded_volatiles_kg{};
  ElementMassKg conservation_error_kg{};
};

struct DifferentiatedImpactorSnapshot {
  ElementPpm silicate_volatiles_ppm{};
  ElementPpm metal_volatiles_ppm{};
  double delta15n_silicate_permil = 0.0;
  double delta15n_metal_permil = 0.0;
  AlloyComposition alloy{};
};

struct SimplifiedRunResult {
  PlanetState final_planet{};
  std::size_t completed_steps = 0;
  bool has_precompute_snapshot = false;
  DifferentiatedImpactorSnapshot precompute_snapshot{};
  ElementMassKg cumulative_eroded_kg{};
  ElementMassKg max_abs_conservation_error_kg{};
  StepResult first_step{};
  StepResult last_step{};
};

struct ElementPartitionKg {
  double atmosphere_kg = 0.0;
  double silicate_kg = 0.0;
  double metal_kg = 0.0;
};

struct Constants {
  double earth_mass_kg = 5.9722e24;
  double earth_radius_m = 6.371e6;
  double gravitational_constant_si = 6.67430e-11;
  double mantle_density_kg_m3 = 3300.0;
  double core_density_kg_m3 = 8000.0;
  double melt_scaling_k = 0.42;
  double melt_scaling_mu = 0.56;
  double melt_scaling_p_sin = 0.66;
  double impact_angle_deg = 45.0;
  double melt_energy_reference_j_kg = 9.0e6;
  double heat_capacity_j_kg_k = 1300.0;
  double latent_heat_j_kg = 718.0e3;
  double preimpact_surface_temperature_k = 1750.0;
  double liquidus_surface_temperature_k = 1950.0;
  double liquidus_slope_low_pressure_k_gpa = 28.3;
  double liquidus_slope_high_pressure_k_gpa = 14.0;
  double liquidus_slope_break_gpa = 60.0;
  double geotherm_gradient_k_m = 1.0e-4;
  double melt_pool_convergence_tolerance_m = 1.0;
  int melt_pool_max_iterations = 100;
  double placeholder_delta_iw_surface = -2.0;

  double fixed_nbo_over_t = 2.5;

  double dn_intercept = -0.38;
  double dn_temperature_coeff = 3370.0;
  double dn_pressure_coeff = 75.0;
  double dn_nbo_coeff = 0.56;
  double dn_delta_iw_coeff = 0.47;
  double dn_xs_coeff = 1.97;
  double dn_xc_coeff = 5.73;
  double dn_xni_coeff = 3.88;
  double dn_xsi_coeff = 5.52;

  double delta15n_a = -1.6e7;
  double delta15n_b = -1.1e7;

  double dh_intercept = 2.31;
  double dh_temperature_coeff = -1542.0;
  double dh_pressure_coeff = -38.9;
  double dh_xc_coeff = 6.69;

  double oxygen_lnK_temperature_coeff = -6910.0;
  double oxygen_lnK_intercept = 3.52;
  double oxygen_lnK_pressure_coeff = 71.0;

  double dc_intercept = 1.49;
  double dc_temperature_coeff = 3000.0;
  double dc_pressure_coeff = -235.0;
  double dc_xs_coeff = 9.6;
  double dc_xo_coeff = -19.5;
  double dc_nbo_coeff = -0.118;
  double dc_delta_iw_coeff = -0.238;
  double placeholder_atmosphere_fraction_carbon = 0.10;
  double placeholder_atmosphere_fraction_hydrogen = 0.25;
  double placeholder_atmosphere_fraction_nitrogen = 0.10;
  double placeholder_erosion_fraction_precompute = 1.0e-4;
  double placeholder_erosion_fraction_giant_impact = 0.05;
  double placeholder_erosion_fraction_late_veneer = 0.01;

  double initial_planet_mass_earth = 0.01;
  double precompute_end_mass_earth = 0.10;
  double precompute_step_mass_earth = 1.0e-5;
  double standard_gi_mass_earth = 0.10;
  double theia_mass_earth = 0.095;
  double late_veneer_mass_earth = 0.005;
  double final_planet_mass_earth = 1.0;
  double representative_metallic_fe_fraction = 0.325;

  ElementPpm ec_bulk_volatiles_ppm{
      4000.0,
      400.0,
      250.0,
  };
  double ec_delta15n_permil = -30.0;

  ElementPpm ci_bulk_volatiles_ppm{
      35000.0,
      6900.0,
      1500.0,
  };
  double ci_delta15n_permil = 40.0;

  AlloyComposition representative_ec_alloy{};
};

[[nodiscard]] double EarthMassToKg(double mass_earth, const Constants& constants) {
  return mass_earth * constants.earth_mass_kg;
}

[[nodiscard]] double KgToEarthMass(double mass_kg, const Constants& constants) {
  return mass_kg / constants.earth_mass_kg;
}

[[nodiscard]] double PpmToMassFraction(double concentration_ppm) {
  return concentration_ppm * 1.0e-6;
}

[[nodiscard]] double MassFractionToPpm(double mass_fraction) {
  return mass_fraction * 1.0e6;
}

[[nodiscard]] ElementMassKg PpmToElementMassKg(const ElementPpm& concentrations_ppm,
                                               double host_mass_kg) {
  if (host_mass_kg < 0.0) {
    throw std::runtime_error("PpmToElementMassKg requires non-negative host mass.");
  }

  return {
      PpmToMassFraction(concentrations_ppm.carbon_ppm) * host_mass_kg,
      PpmToMassFraction(concentrations_ppm.hydrogen_ppm) * host_mass_kg,
      PpmToMassFraction(concentrations_ppm.nitrogen_ppm) * host_mass_kg,
  };
}

[[nodiscard]] ElementPpm ElementMassKgToPpm(const ElementMassKg& element_masses_kg,
                                            double host_mass_kg) {
  if (host_mass_kg < 0.0) {
    throw std::runtime_error("ElementMassKgToPpm requires non-negative host mass.");
  }
  if (host_mass_kg == 0.0) {
    return {};
  }

  return {
      MassFractionToPpm(element_masses_kg.carbon_kg / host_mass_kg),
      MassFractionToPpm(element_masses_kg.hydrogen_kg / host_mass_kg),
      MassFractionToPpm(element_masses_kg.nitrogen_kg / host_mass_kg),
  };
}

[[nodiscard]] ElementMassKg AddElementMassKg(const ElementMassKg& left,
                                             const ElementMassKg& right) {
  return {
      left.carbon_kg + right.carbon_kg,
      left.hydrogen_kg + right.hydrogen_kg,
      left.nitrogen_kg + right.nitrogen_kg,
  };
}

[[nodiscard]] ElementMassKg SubtractElementMassKg(const ElementMassKg& left,
                                                  const ElementMassKg& right) {
  return {
      left.carbon_kg - right.carbon_kg,
      left.hydrogen_kg - right.hydrogen_kg,
      left.nitrogen_kg - right.nitrogen_kg,
  };
}

[[nodiscard]] ElementMassKg MaxAbsElementMassKg(const ElementMassKg& left,
                                                const ElementMassKg& right) {
  return {
      std::max(std::fabs(left.carbon_kg), std::fabs(right.carbon_kg)),
      std::max(std::fabs(left.hydrogen_kg), std::fabs(right.hydrogen_kg)),
      std::max(std::fabs(left.nitrogen_kg), std::fabs(right.nitrogen_kg)),
  };
}

[[nodiscard]] double MaxComponentAbs(const ElementMassKg& element_masses_kg) {
  return std::max({std::fabs(element_masses_kg.carbon_kg),
                   std::fabs(element_masses_kg.hydrogen_kg),
                   std::fabs(element_masses_kg.nitrogen_kg)});
}

[[nodiscard]] double Pi() {
  return std::acos(-1.0);
}

[[nodiscard]] double DegreesToRadians(double degrees) {
  return degrees * Pi() / 180.0;
}

[[nodiscard]] double SafeLog10(double value, const char* label) {
  if (!(value > 0.0)) {
    throw std::runtime_error(std::string(label) + " must be positive before log10.");
  }
  return std::log10(value);
}

[[nodiscard]] double SafeOneMinus(double value, const char* label) {
  const double remainder = 1.0 - value;
  if (!(remainder > 0.0)) {
    throw std::runtime_error(std::string(label) + " must be < 1.");
  }
  return remainder;
}

[[nodiscard]] double ComputeNitrogenPartitionCoefficient(double temperature_k,
                                                         double pressure_gpa,
                                                         double delta_iw_eq,
                                                         const AlloyComposition& alloy,
                                                         const Constants& constants) {
  alloy.Validate();
  if (temperature_k <= 0.0) {
    throw std::runtime_error("ComputeNitrogenPartitionCoefficient requires T > 0.");
  }

  const double log10_dn = constants.dn_intercept +
                          constants.dn_temperature_coeff / temperature_k +
                          constants.dn_pressure_coeff * pressure_gpa / temperature_k +
                          constants.dn_nbo_coeff * constants.fixed_nbo_over_t +
                          constants.dn_delta_iw_coeff * delta_iw_eq -
                          constants.dn_xs_coeff * alloy.x_s -
                          constants.dn_xc_coeff * alloy.x_c -
                          constants.dn_xni_coeff * alloy.x_ni -
                          constants.dn_xsi_coeff * alloy.x_si;
  return std::pow(10.0, log10_dn);
}

[[nodiscard]] double ComputeNitrogenIsotopeFractionationPermil(double temperature_k,
                                                               double delta_iw_eq,
                                                               const Constants& constants) {
  if (temperature_k <= 0.0) {
    throw std::runtime_error("ComputeNitrogenIsotopeFractionationPermil requires T > 0.");
  }
  return constants.delta15n_a / (temperature_k * temperature_k) +
         constants.delta15n_b * delta_iw_eq / (temperature_k * temperature_k);
}

[[nodiscard]] double ComputeHydrogenPartitionCoefficient(double temperature_k,
                                                         double pressure_gpa,
                                                         double delta_iw_eq,
                                                         const AlloyComposition& alloy,
                                                         const Constants& constants) {
  alloy.Validate();
  if (temperature_k <= 0.0) {
    throw std::runtime_error("ComputeHydrogenPartitionCoefficient requires T > 0.");
  }

  const double log10_dh = constants.dh_intercept +
                          constants.dh_temperature_coeff / temperature_k +
                          constants.dh_pressure_coeff * pressure_gpa / temperature_k -
                          0.5 * delta_iw_eq +
                          constants.dh_xc_coeff * alloy.x_c;
  return std::pow(10.0, log10_dh);
}

[[nodiscard]] double ComputeOxygenMoleFractionInMetal(double temperature_k,
                                                      double pressure_gpa,
                                                      double x_feo_silicate,
                                                      const AlloyComposition& alloy,
                                                      const Constants& constants) {
  alloy.Validate();
  if (temperature_k <= 0.0) {
    throw std::runtime_error("ComputeOxygenMoleFractionInMetal requires T > 0.");
  }
  if (x_feo_silicate < 0.0) {
    throw std::runtime_error("ComputeOxygenMoleFractionInMetal requires X_FeO >= 0.");
  }

  const double x_fe = alloy.x_fe();
  const double ln_k = constants.oxygen_lnK_temperature_coeff / temperature_k +
                      constants.oxygen_lnK_intercept +
                      constants.oxygen_lnK_pressure_coeff * pressure_gpa / temperature_k;
  const double equilibrium_constant = std::exp(ln_k);
  return equilibrium_constant * x_feo_silicate / std::max(x_fe, 1.0e-12);
}

[[nodiscard]] double ComputeCarbonPartitionCoefficient(double temperature_k,
                                                       double pressure_gpa,
                                                       double delta_iw_eq,
                                                       double x_s_metal,
                                                       double x_o_metal,
                                                       const Constants& constants) {
  if (temperature_k <= 0.0) {
    throw std::runtime_error("ComputeCarbonPartitionCoefficient requires T > 0.");
  }
  if (x_s_metal < 0.0 || x_o_metal < 0.0) {
    throw std::runtime_error("ComputeCarbonPartitionCoefficient requires non-negative X_S and X_O.");
  }

  const double log10_dc = constants.dc_intercept +
                          constants.dc_temperature_coeff / temperature_k +
                          constants.dc_pressure_coeff * pressure_gpa / temperature_k +
                          constants.dc_xs_coeff * x_s_metal +
                          constants.dc_xo_coeff * x_o_metal +
                          constants.dc_nbo_coeff * constants.fixed_nbo_over_t +
                          constants.dc_delta_iw_coeff * delta_iw_eq;
  return std::pow(10.0, log10_dc);
}

[[nodiscard]] AlloyComposition UpdatePrecomputeAlloyFromCoreCarbonPpm(
    double core_carbon_ppm, const AlloyComposition& reference_alloy) {
  reference_alloy.Validate();
  if (core_carbon_ppm < 0.0) {
    throw std::runtime_error("UpdatePrecomputeAlloyFromCoreCarbonPpm requires non-negative C ppm.");
  }

  constexpr double kAtomicMassCarbon = 12.011;
  constexpr double kAtomicMassSulfur = 32.06;
  constexpr double kAtomicMassNickel = 58.6934;
  constexpr double kAtomicMassSilicon = 28.085;
  constexpr double kAtomicMassIron = 55.845;

  const double fixed_minor_sum =
      reference_alloy.x_s + reference_alloy.x_ni + reference_alloy.x_si;
  const double max_x_c = std::max(0.0, 1.0 - fixed_minor_sum - 1.0e-12);
  const double carbon_mass_fraction = std::clamp(PpmToMassFraction(core_carbon_ppm), 0.0, 0.20);

  const double non_carbon_molar_mass =
      reference_alloy.x_s * kAtomicMassSulfur +
      reference_alloy.x_ni * kAtomicMassNickel +
      reference_alloy.x_si * kAtomicMassSilicon +
      (1.0 - fixed_minor_sum) * kAtomicMassIron;
  const double denominator =
      kAtomicMassCarbon - carbon_mass_fraction * (kAtomicMassCarbon - kAtomicMassIron);
  if (!(denominator > 0.0)) {
    throw std::runtime_error("Carbon mole-fraction denominator became non-positive.");
  }

  AlloyComposition updated = reference_alloy;
  updated.x_c = std::clamp(carbon_mass_fraction * non_carbon_molar_mass / denominator,
                           0.0, max_x_c);
  updated.Validate();
  return updated;
}

[[nodiscard]] double RadiusFromMassKg(double mass_kg, const Constants& constants) {
  if (mass_kg <= 0.0) {
    throw std::runtime_error("RadiusFromMassKg requires positive mass.");
  }

  // Sakuraba.C uses the Seager et al. rocky-planet fit in cgs units.
  const double mass_ratio_to_earth = KgToEarthMass(mass_kg, constants);
  const double scaled_mass = mass_ratio_to_earth / 6.41;
  const double log10_radius_in_earth_units =
      -0.20945 + (1.0 / 3.0) * std::log10(scaled_mass) -
      0.0804 * std::pow(scaled_mass, 0.394) + std::log10(3.29);
  return std::pow(10.0, log10_radius_in_earth_units) * constants.earth_radius_m;
}

[[nodiscard]] double Gravity(double mass_kg, double radius_m,
                             const Constants& constants) {
  if (mass_kg <= 0.0) {
    throw std::runtime_error("Gravity requires positive mass.");
  }
  if (radius_m <= 0.0) {
    throw std::runtime_error("Gravity requires positive radius.");
  }

  return constants.gravitational_constant_si * mass_kg / (radius_m * radius_m);
}

[[nodiscard]] double EscapeVelocity(double mass_kg, double radius_m,
                                    const Constants& constants) {
  if (mass_kg <= 0.0) {
    throw std::runtime_error("EscapeVelocity requires positive mass.");
  }
  if (radius_m <= 0.0) {
    throw std::runtime_error("EscapeVelocity requires positive radius.");
  }

  return std::sqrt(2.0 * constants.gravitational_constant_si * mass_kg / radius_m);
}

[[nodiscard]] double SphereRadiusFromMassAndDensity(double mass_kg,
                                                    double density_kg_m3) {
  if (mass_kg < 0.0) {
    throw std::runtime_error("SphereRadiusFromMassAndDensity requires non-negative mass.");
  }
  if (density_kg_m3 <= 0.0) {
    throw std::runtime_error(
        "SphereRadiusFromMassAndDensity requires positive density.");
  }
  if (mass_kg == 0.0) {
    return 0.0;
  }

  const double volume_m3 = mass_kg / density_kg_m3;
  return std::cbrt((3.0 * volume_m3) / (4.0 * std::acos(-1.0)));
}

[[nodiscard]] double MantleThickness(const PlanetState& planet,
                                     const Constants& constants) {
  const double planet_radius_m = RadiusFromMassKg(planet.total_mass_kg(), constants);
  const double core_radius_m =
      SphereRadiusFromMassAndDensity(planet.core.mass_kg, constants.core_density_kg_m3);
  return std::max(0.0, planet_radius_m - core_radius_m);
}

[[nodiscard]] double BulkDensityFromMetalFraction(double metallic_fe_fraction,
                                                  const Constants& constants) {
  if (metallic_fe_fraction < 0.0 || metallic_fe_fraction > 1.0) {
    throw std::runtime_error("BulkDensityFromMetalFraction requires f in [0, 1].");
  }

  const double metal_volume_per_kg = metallic_fe_fraction / constants.core_density_kg_m3;
  const double silicate_volume_per_kg =
      (1.0 - metallic_fe_fraction) / constants.mantle_density_kg_m3;
  return 1.0 / (metal_volume_per_kg + silicate_volume_per_kg);
}

[[nodiscard]] double SphereDiameterFromMassAndDensity(double mass_kg,
                                                      double density_kg_m3) {
  return 2.0 * SphereRadiusFromMassAndDensity(mass_kg, density_kg_m3);
}

[[nodiscard]] double PressureGPaFromDepth(double depth_m, double gravity_m_s2,
                                          double density_kg_m3) {
  if (depth_m < 0.0) {
    throw std::runtime_error("PressureGPaFromDepth requires non-negative depth.");
  }
  if (gravity_m_s2 < 0.0) {
    throw std::runtime_error("PressureGPaFromDepth requires non-negative gravity.");
  }
  if (density_kg_m3 <= 0.0) {
    throw std::runtime_error("PressureGPaFromDepth requires positive density.");
  }
  return density_kg_m3 * gravity_m_s2 * depth_m * 1.0e-9;
}

[[nodiscard]] double LiquidusSlopeKPerGPa(double pressure_gpa,
                                          const Constants& constants) {
  return pressure_gpa < constants.liquidus_slope_break_gpa
             ? constants.liquidus_slope_low_pressure_k_gpa
             : constants.liquidus_slope_high_pressure_k_gpa;
}

[[nodiscard]] double LiquidusTemperatureK(double pressure_gpa,
                                          const Constants& constants) {
  return constants.liquidus_surface_temperature_k +
         LiquidusSlopeKPerGPa(pressure_gpa, constants) * pressure_gpa;
}

[[nodiscard]] double EffectiveMeltEnergyJPerKg(double depth_m, double gravity_m_s2,
                                               const Constants& constants) {
  const double pressure_gpa =
      PressureGPaFromDepth(depth_m, gravity_m_s2, constants.mantle_density_kg_m3);
  const double geotherm_temperature_k =
      constants.preimpact_surface_temperature_k +
      constants.geotherm_gradient_k_m * depth_m;
  const double liquidus_temperature_k = LiquidusTemperatureK(pressure_gpa, constants);
  const double denominator =
      constants.heat_capacity_j_kg_k * liquidus_temperature_k +
      constants.latent_heat_j_kg;
  if (denominator <= 0.0) {
    throw std::runtime_error("EffectiveMeltEnergyJPerKg encountered non-positive denominator.");
  }

  const double melt_energy =
      constants.melt_energy_reference_j_kg *
      (1.0 - constants.heat_capacity_j_kg_k * geotherm_temperature_k / denominator);
  if (melt_energy <= 0.0) {
    throw std::runtime_error("EffectiveMeltEnergyJPerKg became non-positive.");
  }
  return melt_energy;
}

[[nodiscard]] double MeltVolumeM3(double melt_energy_j_kg,
                                  double impactor_density_kg_m3,
                                  double target_density_kg_m3,
                                  double impactor_diameter_m,
                                  double impact_velocity_m_s,
                                  const Constants& constants) {
  if (melt_energy_j_kg <= 0.0) {
    throw std::runtime_error("MeltVolumeM3 requires positive melt energy.");
  }
  if (impactor_density_kg_m3 <= 0.0 || target_density_kg_m3 <= 0.0) {
    throw std::runtime_error("MeltVolumeM3 requires positive densities.");
  }
  if (impactor_diameter_m <= 0.0 || impact_velocity_m_s <= 0.0) {
    throw std::runtime_error("MeltVolumeM3 requires positive diameter and velocity.");
  }

  return (Pi() / 6.0) * constants.melt_scaling_k *
         std::pow(melt_energy_j_kg, -1.5 * constants.melt_scaling_mu) *
         (impactor_density_kg_m3 / target_density_kg_m3) *
         std::pow(impactor_diameter_m, 3.0) *
         std::pow(impact_velocity_m_s, 3.0 * constants.melt_scaling_mu) *
         std::pow(std::sin(DegreesToRadians(constants.impact_angle_deg)),
                  2.0 * constants.melt_scaling_p_sin);
}

[[nodiscard]] double SphericalCapVolumeM3(double depth_m, double radius_m) {
  if (radius_m <= 0.0) {
    throw std::runtime_error("SphericalCapVolumeM3 requires positive radius.");
  }
  if (depth_m < 0.0) {
    throw std::runtime_error("SphericalCapVolumeM3 requires non-negative depth.");
  }

  return (Pi() / 3.0) * depth_m * depth_m * (3.0 * radius_m - depth_m);
}

[[nodiscard]] double SolveSphericalCapDepthM(double volume_m3, double radius_m,
                                             const Constants& constants) {
  if (volume_m3 < 0.0) {
    throw std::runtime_error("SolveSphericalCapDepthM requires non-negative volume.");
  }
  if (radius_m <= 0.0) {
    throw std::runtime_error("SolveSphericalCapDepthM requires positive radius.");
  }
  if (volume_m3 == 0.0) {
    return 0.0;
  }

  const double full_sphere_volume_m3 = (4.0 / 3.0) * Pi() * std::pow(radius_m, 3.0);
  if (volume_m3 >= full_sphere_volume_m3) {
    return 2.0 * radius_m;
  }

  double lower_depth_m = 0.0;
  double upper_depth_m = 2.0 * radius_m;
  double depth_m = std::min(
      upper_depth_m,
      std::max(constants.melt_pool_convergence_tolerance_m,
               std::sqrt(volume_m3 / (Pi() * radius_m))));

  for (int iteration = 0; iteration < constants.melt_pool_max_iterations; ++iteration) {
    const double volume_residual_m3 = SphericalCapVolumeM3(depth_m, radius_m) - volume_m3;
    if (volume_residual_m3 > 0.0) {
      upper_depth_m = depth_m;
    } else {
      lower_depth_m = depth_m;
    }

    const double derivative_m2 = Pi() * depth_m * (2.0 * radius_m - depth_m);
    double next_depth_m = 0.5 * (lower_depth_m + upper_depth_m);
    if (derivative_m2 > 0.0) {
      const double newton_depth_m = depth_m - volume_residual_m3 / derivative_m2;
      if (std::isfinite(newton_depth_m) && newton_depth_m > lower_depth_m &&
          newton_depth_m < upper_depth_m) {
        next_depth_m = newton_depth_m;
      }
    }

    if (std::fabs(next_depth_m - depth_m) <
        constants.melt_pool_convergence_tolerance_m) {
      return next_depth_m;
    }
    depth_m = next_depth_m;
  }

  throw std::runtime_error("SolveSphericalCapDepthM did not converge.");
}

[[nodiscard]] MeltPoolState ComputeMeltPoolState(const PlanetState& planet,
                                                 const Impactor& impactor,
                                                 const Constants& constants) {
  planet.Validate();
  impactor.Validate();

  MeltPoolState melt_pool{};
  melt_pool.target_radius_m = RadiusFromMassKg(planet.total_mass_kg(), constants);
  melt_pool.core_radius_m = SphereRadiusFromMassAndDensity(
      planet.core.mass_kg, constants.core_density_kg_m3);
  melt_pool.mantle_thickness_m = MantleThickness(planet, constants);
  melt_pool.combined_body_mass_kg = planet.total_mass_kg() + impactor.mass_kg;
  melt_pool.combined_body_radius_m =
      RadiusFromMassKg(melt_pool.combined_body_mass_kg, constants);
  melt_pool.gravity_m_s2 = Gravity(melt_pool.combined_body_mass_kg,
                                   melt_pool.combined_body_radius_m, constants);
  melt_pool.impact_velocity_m_s =
      EscapeVelocity(melt_pool.combined_body_mass_kg,
                     melt_pool.combined_body_radius_m, constants);
  melt_pool.impactor_density_kg_m3 =
      BulkDensityFromMetalFraction(impactor.metallic_fe_fraction, constants);
  melt_pool.impactor_diameter_m = SphereDiameterFromMassAndDensity(
      impactor.mass_kg, melt_pool.impactor_density_kg_m3);
  melt_pool.impactor_silicate_mass_kg =
      (1.0 - impactor.metallic_fe_fraction) * impactor.mass_kg;
  melt_pool.impactor_metal_mass_kg =
      impactor.metallic_fe_fraction * impactor.mass_kg;

  double depth_guess_m = 0.0;
  for (int iteration = 0; iteration < constants.melt_pool_max_iterations; ++iteration) {
    const double melt_energy_j_kg =
        EffectiveMeltEnergyJPerKg(depth_guess_m, melt_pool.gravity_m_s2, constants);
    const double target_melt_volume_m3 =
        MeltVolumeM3(melt_energy_j_kg, melt_pool.impactor_density_kg_m3,
                     constants.mantle_density_kg_m3, melt_pool.impactor_diameter_m,
                     melt_pool.impact_velocity_m_s, constants);
    const double next_depth_m = SolveSphericalCapDepthM(
        target_melt_volume_m3, melt_pool.target_radius_m, constants);

    melt_pool.iteration_count = iteration + 1;
    melt_pool.melt_energy_j_kg = melt_energy_j_kg;
    melt_pool.target_melt_volume_m3 = target_melt_volume_m3;
    melt_pool.raw_depth_m = next_depth_m;

    if (std::fabs(next_depth_m - depth_guess_m) <
        constants.melt_pool_convergence_tolerance_m) {
      break;
    }

    depth_guess_m = next_depth_m;
    if (iteration == constants.melt_pool_max_iterations - 1) {
      throw std::runtime_error("ComputeMeltPoolState fixed-point iteration did not converge.");
    }
  }

  melt_pool.depth_m =
      std::min(melt_pool.raw_depth_m, melt_pool.mantle_thickness_m);
  melt_pool.pressure_gpa = PressureGPaFromDepth(
      melt_pool.depth_m, melt_pool.gravity_m_s2, constants.mantle_density_kg_m3);
  melt_pool.temperature_k = LiquidusTemperatureK(melt_pool.pressure_gpa, constants);
  melt_pool.target_melt_mass_kg =
      std::min(melt_pool.target_melt_volume_m3 * constants.mantle_density_kg_m3,
               planet.mantle.mass_kg);
  melt_pool.equilibrating_silicate_mass_kg =
      melt_pool.target_melt_mass_kg + melt_pool.impactor_silicate_mass_kg;
  return melt_pool;
}

void ApplyPlaceholderAtmosphereSpeciation(AtmosphereState* atmosphere,
                                          const Constants& constants) {
  if (atmosphere == nullptr) {
    throw std::runtime_error("ApplyPlaceholderAtmosphereSpeciation received null atmosphere.");
  }

  if (constants.placeholder_delta_iw_surface > -1.0) {
    atmosphere->carbon_species = CarbonGasSpecies::kCO2;
    atmosphere->hydrogen_species = HydrogenGasSpecies::kH2O;
    atmosphere->nitrogen_species = NitrogenGasSpecies::kN2;
  } else if (constants.placeholder_delta_iw_surface > -3.0) {
    atmosphere->carbon_species = CarbonGasSpecies::kCO;
    atmosphere->hydrogen_species = HydrogenGasSpecies::kH2;
    atmosphere->nitrogen_species = NitrogenGasSpecies::kN2;
  } else {
    atmosphere->carbon_species = CarbonGasSpecies::kCH4;
    atmosphere->hydrogen_species = HydrogenGasSpecies::kH2;
    atmosphere->nitrogen_species = NitrogenGasSpecies::kN2;
  }
}

[[nodiscard]] double PlaceholderErosionFraction(GrowthPhase phase,
                                                const Constants& constants) {
  switch (phase) {
    case GrowthPhase::kPrecomputeAccretion:
      return constants.placeholder_erosion_fraction_precompute;
    case GrowthPhase::kGiantImpact:
      return constants.placeholder_erosion_fraction_giant_impact;
    case GrowthPhase::kLateVeneer:
      return constants.placeholder_erosion_fraction_late_veneer;
  }
  throw std::runtime_error("Unknown GrowthPhase in PlaceholderErosionFraction.");
}

[[nodiscard]] Impactor ResolveImpactorForStep(
    const Impactor& impactor,
    const std::optional<DifferentiatedImpactorSnapshot>& precompute_snapshot) {
  Impactor resolved = impactor;

  switch (resolved.silicate_inventory_source) {
    case InventorySource::kNone:
      resolved.silicate_volatiles_ppm = {};
      resolved.delta15n_silicate_permil = 0.0;
      break;
    case InventorySource::kEcBulk:
    case InventorySource::kCiBulk:
      break;
    case InventorySource::kPrecomputeOutputAt0p10EarthMass:
      if (!precompute_snapshot.has_value()) {
        throw std::runtime_error("Impactor requires precompute silicate snapshot.");
      }
      resolved.silicate_volatiles_ppm = precompute_snapshot->silicate_volatiles_ppm;
      resolved.delta15n_silicate_permil =
          precompute_snapshot->delta15n_silicate_permil;
      break;
    case InventorySource::kPrecomputeCoreAt0p10EarthMass:
      throw std::runtime_error(
          "Silicate inventory source cannot be precompute core snapshot.");
  }

  switch (resolved.metal_inventory_source) {
    case InventorySource::kNone:
      resolved.metal_volatiles_ppm = {};
      resolved.delta15n_metal_permil = std::numeric_limits<double>::quiet_NaN();
      break;
    case InventorySource::kEcBulk:
    case InventorySource::kCiBulk:
      break;
    case InventorySource::kPrecomputeCoreAt0p10EarthMass:
      if (!precompute_snapshot.has_value()) {
        throw std::runtime_error("Impactor requires precompute core snapshot.");
      }
      resolved.metal_volatiles_ppm = precompute_snapshot->metal_volatiles_ppm;
      resolved.delta15n_metal_permil = precompute_snapshot->delta15n_metal_permil;
      resolved.alloy = precompute_snapshot->alloy;
      break;
    case InventorySource::kPrecomputeOutputAt0p10EarthMass:
      throw std::runtime_error(
          "Metal inventory source cannot be precompute silicate snapshot.");
  }

  resolved.Validate();
  return resolved;
}

[[nodiscard]] DifferentiatedImpactorSnapshot CaptureDifferentiatedImpactorSnapshot(
    const PlanetState& planet, const AlloyComposition& alloy) {
  planet.Validate();
  alloy.Validate();
  return {
      planet.mantle.volatiles_ppm,
      planet.core.volatiles_ppm,
      planet.mantle.delta15n_permil,
      planet.core.delta15n_permil,
      alloy,
  };
}

[[nodiscard]] ElementPartitionKg PartitionElementMass(double total_element_kg,
                                                      double silicate_mass_kg,
                                                      double metal_mass_kg,
                                                      double partition_coefficient,
                                                      double atmosphere_fraction) {
  if (total_element_kg < 0.0) {
    throw std::runtime_error("PartitionElementMass requires non-negative total mass.");
  }
  if (silicate_mass_kg < 0.0 || metal_mass_kg < 0.0) {
    throw std::runtime_error("PartitionElementMass requires non-negative host masses.");
  }
  if (partition_coefficient < 0.0) {
    throw std::runtime_error("PartitionElementMass requires non-negative partition coefficient.");
  }

  const double clamped_atmosphere_fraction =
      std::clamp(atmosphere_fraction, 0.0, 1.0);
  ElementPartitionKg partition{};
  partition.atmosphere_kg = total_element_kg * clamped_atmosphere_fraction;

  const double retained_non_atmosphere_kg =
      total_element_kg - partition.atmosphere_kg;
  if (retained_non_atmosphere_kg <= 0.0) {
    return partition;
  }

  if (silicate_mass_kg <= 0.0 && metal_mass_kg <= 0.0) {
    partition.atmosphere_kg = total_element_kg;
    return partition;
  }
  if (metal_mass_kg <= 0.0) {
    partition.silicate_kg = retained_non_atmosphere_kg;
    return partition;
  }
  if (silicate_mass_kg <= 0.0) {
    partition.metal_kg = retained_non_atmosphere_kg;
    return partition;
  }

  const double metal_to_silicate_mass_ratio = metal_mass_kg / silicate_mass_kg;
  const double silicate_share_kg =
      retained_non_atmosphere_kg /
      (1.0 + partition_coefficient * metal_to_silicate_mass_ratio);
  partition.silicate_kg = silicate_share_kg;
  partition.metal_kg = retained_non_atmosphere_kg - silicate_share_kg;
  return partition;
}

[[nodiscard]] StepResult RunOneImpactStep(std::size_t index, const PlanetState& planet_before,
                                          const Impactor& resolved_impactor,
                                          const Constants& constants) {
  planet_before.Validate();
  resolved_impactor.Validate();

  StepResult result{};
  result.index = index;
  result.resolved_impactor = resolved_impactor;
  result.erosion_fraction =
      PlaceholderErosionFraction(resolved_impactor.phase, constants);
  result.melt_pool =
      ComputeMeltPoolState(planet_before, resolved_impactor, constants);
  result.dn_metal_silicate = ComputeNitrogenPartitionCoefficient(
      result.melt_pool.temperature_k, result.melt_pool.pressure_gpa,
      planet_before.redox.delta_iw_eq, resolved_impactor.alloy, constants);
  result.dh_metal_silicate = ComputeHydrogenPartitionCoefficient(
      result.melt_pool.temperature_k, result.melt_pool.pressure_gpa,
      planet_before.redox.delta_iw_eq, resolved_impactor.alloy, constants);
  result.xo_metal = ComputeOxygenMoleFractionInMetal(
      result.melt_pool.temperature_k, result.melt_pool.pressure_gpa,
      planet_before.redox.x_feo_silicate, resolved_impactor.alloy, constants);
  result.dc_metal_silicate = ComputeCarbonPartitionCoefficient(
      result.melt_pool.temperature_k, result.melt_pool.pressure_gpa,
      planet_before.redox.delta_iw_eq, resolved_impactor.alloy.x_s, result.xo_metal,
      constants);

  const double target_unmelted_mantle_mass_kg =
      planet_before.mantle.mass_kg - result.melt_pool.target_melt_mass_kg;
  if (target_unmelted_mantle_mass_kg < -1.0e-6) {
    throw std::runtime_error("Melt pool consumed more mantle mass than available.");
  }

  const ElementMassKg old_mantle_volatiles_kg =
      PpmToElementMassKg(planet_before.mantle.volatiles_ppm,
                         planet_before.mantle.mass_kg);
  const ElementMassKg old_core_volatiles_kg =
      PpmToElementMassKg(planet_before.core.volatiles_ppm, planet_before.core.mass_kg);
  const ElementMassKg old_atmosphere_volatiles_kg =
      planet_before.atmosphere.elements_kg;

  const ElementMassKg target_melt_volatiles_kg =
      PpmToElementMassKg(planet_before.mantle.volatiles_ppm,
                         result.melt_pool.target_melt_mass_kg);
  const ElementMassKg target_unmelted_volatiles_kg =
      PpmToElementMassKg(planet_before.mantle.volatiles_ppm,
                         std::max(0.0, target_unmelted_mantle_mass_kg));
  const ElementMassKg impactor_silicate_volatiles_kg =
      PpmToElementMassKg(resolved_impactor.silicate_volatiles_ppm,
                         result.melt_pool.impactor_silicate_mass_kg);
  const ElementMassKg impactor_metal_volatiles_kg =
      PpmToElementMassKg(resolved_impactor.metal_volatiles_ppm,
                         result.melt_pool.impactor_metal_mass_kg);

  const ElementMassKg mixed_silicate_volatiles_kg = AddElementMassKg(
      target_melt_volatiles_kg, impactor_silicate_volatiles_kg);

  double mixed_delta15n_permil = planet_before.mantle.delta15n_permil;
  const double mixed_nitrogen_mass_kg = mixed_silicate_volatiles_kg.nitrogen_kg;
  if (mixed_nitrogen_mass_kg > 0.0) {
    double numerator = 0.0;
    if (target_melt_volatiles_kg.nitrogen_kg > 0.0) {
      numerator += target_melt_volatiles_kg.nitrogen_kg *
                   planet_before.mantle.delta15n_permil;
    }
    if (impactor_silicate_volatiles_kg.nitrogen_kg > 0.0) {
      numerator += impactor_silicate_volatiles_kg.nitrogen_kg *
                   resolved_impactor.delta15n_silicate_permil;
    }
    mixed_delta15n_permil = numerator / mixed_nitrogen_mass_kg;
  }

  const ElementPartitionKg carbon_partition = PartitionElementMass(
      old_atmosphere_volatiles_kg.carbon_kg + mixed_silicate_volatiles_kg.carbon_kg +
          impactor_metal_volatiles_kg.carbon_kg,
      result.melt_pool.equilibrating_silicate_mass_kg,
      result.melt_pool.impactor_metal_mass_kg, result.dc_metal_silicate,
      constants.placeholder_atmosphere_fraction_carbon);
  const ElementPartitionKg hydrogen_partition = PartitionElementMass(
      old_atmosphere_volatiles_kg.hydrogen_kg +
          mixed_silicate_volatiles_kg.hydrogen_kg +
          impactor_metal_volatiles_kg.hydrogen_kg,
      result.melt_pool.equilibrating_silicate_mass_kg,
      result.melt_pool.impactor_metal_mass_kg, result.dh_metal_silicate,
      constants.placeholder_atmosphere_fraction_hydrogen);
  const ElementPartitionKg nitrogen_partition = PartitionElementMass(
      old_atmosphere_volatiles_kg.nitrogen_kg +
          mixed_silicate_volatiles_kg.nitrogen_kg +
          impactor_metal_volatiles_kg.nitrogen_kg,
      result.melt_pool.equilibrating_silicate_mass_kg,
      result.melt_pool.impactor_metal_mass_kg, result.dn_metal_silicate,
      constants.placeholder_atmosphere_fraction_nitrogen);

  const ElementMassKg silicate_after_equilibrium_kg{
      carbon_partition.silicate_kg,
      hydrogen_partition.silicate_kg,
      nitrogen_partition.silicate_kg,
  };
  const ElementMassKg metal_after_equilibrium_kg{
      carbon_partition.metal_kg,
      hydrogen_partition.metal_kg,
      nitrogen_partition.metal_kg,
  };
  const ElementMassKg atmosphere_after_equilibrium_kg{
      carbon_partition.atmosphere_kg,
      hydrogen_partition.atmosphere_kg,
      nitrogen_partition.atmosphere_kg,
  };

  double delta15n_bulk_permil = mixed_delta15n_permil;
  const double total_n_before_equilibrium_kg =
      old_atmosphere_volatiles_kg.nitrogen_kg +
      mixed_silicate_volatiles_kg.nitrogen_kg +
      impactor_metal_volatiles_kg.nitrogen_kg;
  if (total_n_before_equilibrium_kg > 0.0) {
    double numerator = 0.0;
    if (mixed_silicate_volatiles_kg.nitrogen_kg > 0.0) {
      numerator += mixed_silicate_volatiles_kg.nitrogen_kg * mixed_delta15n_permil;
    }
    if (old_atmosphere_volatiles_kg.nitrogen_kg > 0.0) {
      numerator += old_atmosphere_volatiles_kg.nitrogen_kg *
                   planet_before.atmosphere.delta15n_permil;
    }
    if (impactor_metal_volatiles_kg.nitrogen_kg > 0.0) {
      if (!std::isfinite(resolved_impactor.delta15n_metal_permil)) {
        throw std::runtime_error(
            "Impactor metal carries nitrogen mass but has undefined delta15N.");
      }
      numerator += impactor_metal_volatiles_kg.nitrogen_kg *
                   resolved_impactor.delta15n_metal_permil;
    }
    delta15n_bulk_permil = numerator / total_n_before_equilibrium_kg;
  }

  result.delta15n_silicate_permil = delta15n_bulk_permil;
  result.delta15n_metal_permil = delta15n_bulk_permil;
  if (total_n_before_equilibrium_kg > 0.0) {
    result.delta15n_silicate_permil =
        delta15n_bulk_permil -
        (metal_after_equilibrium_kg.nitrogen_kg / total_n_before_equilibrium_kg) *
            ComputeNitrogenIsotopeFractionationPermil(
                result.melt_pool.temperature_k, planet_before.redox.delta_iw_eq, constants);
    result.delta15n_metal_permil =
        result.delta15n_silicate_permil +
        ComputeNitrogenIsotopeFractionationPermil(
            result.melt_pool.temperature_k, planet_before.redox.delta_iw_eq, constants);
  }

  PlanetState planet_after = planet_before;
  planet_after.redox.delta_iw_surface = constants.placeholder_delta_iw_surface;
  planet_after.redox.oxidized_fe0_mass_kg = 0.0;

  planet_after.core.mass_kg =
      planet_before.core.mass_kg + result.melt_pool.impactor_metal_mass_kg;
  planet_after.mantle.mass_kg =
      planet_before.mantle.mass_kg + result.melt_pool.impactor_silicate_mass_kg;

  const ElementMassKg new_core_volatiles_kg =
      AddElementMassKg(old_core_volatiles_kg, metal_after_equilibrium_kg);
  const ElementMassKg new_mantle_volatiles_kg =
      AddElementMassKg(target_unmelted_volatiles_kg, silicate_after_equilibrium_kg);

  planet_after.core.volatiles_ppm =
      ElementMassKgToPpm(new_core_volatiles_kg, planet_after.core.mass_kg);
  planet_after.mantle.volatiles_ppm =
      ElementMassKgToPpm(new_mantle_volatiles_kg, planet_after.mantle.mass_kg);

  const double new_core_nitrogen_kg = new_core_volatiles_kg.nitrogen_kg;
  if (new_core_nitrogen_kg > 0.0) {
    const double numerator =
        old_core_volatiles_kg.nitrogen_kg * planet_before.core.delta15n_permil +
        metal_after_equilibrium_kg.nitrogen_kg * result.delta15n_metal_permil;
    planet_after.core.delta15n_permil = numerator / new_core_nitrogen_kg;
  }

  const double new_mantle_nitrogen_kg = new_mantle_volatiles_kg.nitrogen_kg;
  if (new_mantle_nitrogen_kg > 0.0) {
    const double numerator =
        target_unmelted_volatiles_kg.nitrogen_kg * planet_before.mantle.delta15n_permil +
        silicate_after_equilibrium_kg.nitrogen_kg * result.delta15n_silicate_permil;
    planet_after.mantle.delta15n_permil = numerator / new_mantle_nitrogen_kg;
  }

  planet_after.atmosphere.elements_kg = atmosphere_after_equilibrium_kg;
  planet_after.atmosphere.delta15n_permil = result.delta15n_silicate_permil;
  ApplyPlaceholderAtmosphereSpeciation(&planet_after.atmosphere, constants);

  result.eroded_volatiles_kg = {
      planet_after.atmosphere.elements_kg.carbon_kg * result.erosion_fraction,
      planet_after.atmosphere.elements_kg.hydrogen_kg * result.erosion_fraction,
      planet_after.atmosphere.elements_kg.nitrogen_kg * result.erosion_fraction,
  };
  planet_after.atmosphere.elements_kg = SubtractElementMassKg(
      planet_after.atmosphere.elements_kg, result.eroded_volatiles_kg);

  const ElementMassKg total_before_kg = AddElementMassKg(
      AddElementMassKg(old_mantle_volatiles_kg, old_core_volatiles_kg),
      AddElementMassKg(old_atmosphere_volatiles_kg,
                       AddElementMassKg(impactor_silicate_volatiles_kg,
                                        impactor_metal_volatiles_kg)));
  const ElementMassKg total_after_with_loss_kg = AddElementMassKg(
      AddElementMassKg(
          PpmToElementMassKg(planet_after.mantle.volatiles_ppm, planet_after.mantle.mass_kg),
          PpmToElementMassKg(planet_after.core.volatiles_ppm, planet_after.core.mass_kg)),
      AddElementMassKg(planet_after.atmosphere.elements_kg, result.eroded_volatiles_kg));
  result.conservation_error_kg =
      SubtractElementMassKg(total_after_with_loss_kg, total_before_kg);

  planet_after.Validate();
  result.planet_after = planet_after;
  return result;
}

[[nodiscard]] SimplifiedRunResult RunScenarioSimplified(
    const PlanetState& initial_planet, const std::vector<Impactor>& scenario,
    const Constants& constants) {
  SimplifiedRunResult run{};
  PlanetState current_planet = initial_planet;
  std::optional<DifferentiatedImpactorSnapshot> precompute_snapshot;
  AlloyComposition current_precompute_alloy = constants.representative_ec_alloy;

  for (std::size_t index = 0; index < scenario.size(); ++index) {
    Impactor resolved_impactor =
        ResolveImpactorForStep(scenario.at(index), precompute_snapshot);
    if (resolved_impactor.phase == GrowthPhase::kPrecomputeAccretion) {
      resolved_impactor.alloy = current_precompute_alloy;
    }
    const StepResult step_result =
        RunOneImpactStep(index, current_planet, resolved_impactor, constants);

    if (index == 0) {
      run.first_step = step_result;
    }
    run.last_step = step_result;
    run.completed_steps = index + 1;
    run.cumulative_eroded_kg =
        AddElementMassKg(run.cumulative_eroded_kg, step_result.eroded_volatiles_kg);
    run.max_abs_conservation_error_kg = MaxAbsElementMassKg(
        run.max_abs_conservation_error_kg, step_result.conservation_error_kg);
    current_planet = step_result.planet_after;

    if (scenario.at(index).phase == GrowthPhase::kPrecomputeAccretion) {
      current_precompute_alloy = UpdatePrecomputeAlloyFromCoreCarbonPpm(
          current_planet.core.volatiles_ppm.carbon_ppm, current_precompute_alloy);
    }

    const bool is_last_precompute_step =
        scenario.at(index).phase == GrowthPhase::kPrecomputeAccretion &&
        (index + 1 == scenario.size() ||
         scenario.at(index + 1).phase != GrowthPhase::kPrecomputeAccretion);
    if (is_last_precompute_step) {
      precompute_snapshot =
          CaptureDifferentiatedImpactorSnapshot(current_planet, current_precompute_alloy);
      run.has_precompute_snapshot = true;
      run.precompute_snapshot = *precompute_snapshot;
    }
  }

  run.final_planet = current_planet;
  return run;
}

[[nodiscard]] PlanetState MakeInitialPlanetState(const Constants& constants) {
  PlanetState initial{};
  initial.mantle.mass_kg =
      EarthMassToKg(constants.initial_planet_mass_earth, constants);
  initial.mantle.volatiles_ppm = constants.ec_bulk_volatiles_ppm;
  initial.mantle.delta15n_permil = constants.ec_delta15n_permil;
  initial.core.mass_kg = 0.0;
  initial.redox.x_feo_silicate = 0.056;
  initial.redox.delta_iw_eq = -2.0;
  initial.Validate();
  return initial;
}

[[nodiscard]] std::string ToString(GrowthPhase phase) {
  switch (phase) {
    case GrowthPhase::kPrecomputeAccretion:
      return "precompute_accretion";
    case GrowthPhase::kGiantImpact:
      return "giant_impact";
    case GrowthPhase::kLateVeneer:
      return "late_veneer";
  }
  throw std::runtime_error("Unknown GrowthPhase.");
}

[[nodiscard]] std::string ToString(ImpactorClass impactor_class) {
  switch (impactor_class) {
    case ImpactorClass::kECLikeUndifferentiated:
      return "ec_like_undifferentiated";
    case ImpactorClass::kECLikeDifferentiated:
      return "ec_like_differentiated";
    case ImpactorClass::kCILikeDifferentiated:
      return "ci_like_differentiated";
  }
  throw std::runtime_error("Unknown ImpactorClass.");
}

[[nodiscard]] std::string ToString(CarbonGasSpecies species) {
  switch (species) {
    case CarbonGasSpecies::kCO2:
      return "CO2";
    case CarbonGasSpecies::kCO:
      return "CO";
    case CarbonGasSpecies::kCH4:
      return "CH4";
  }
  throw std::runtime_error("Unknown CarbonGasSpecies.");
}

[[nodiscard]] std::string ToString(HydrogenGasSpecies species) {
  switch (species) {
    case HydrogenGasSpecies::kH2O:
      return "H2O";
    case HydrogenGasSpecies::kH2:
      return "H2";
  }
  throw std::runtime_error("Unknown HydrogenGasSpecies.");
}

[[nodiscard]] std::string ToString(NitrogenGasSpecies species) {
  switch (species) {
    case NitrogenGasSpecies::kN2:
      return "N2";
  }
  throw std::runtime_error("Unknown NitrogenGasSpecies.");
}

[[nodiscard]] std::string ToString(InventorySource inventory_source) {
  switch (inventory_source) {
    case InventorySource::kNone:
      return "none";
    case InventorySource::kEcBulk:
      return "ec_bulk";
    case InventorySource::kCiBulk:
      return "ci_bulk";
    case InventorySource::kPrecomputeOutputAt0p10EarthMass:
      return "precompute_output_at_0p10";
    case InventorySource::kPrecomputeCoreAt0p10EarthMass:
      return "precompute_core_at_0p10";
  }
  throw std::runtime_error("Unknown InventorySource.");
}

[[nodiscard]] std::string ZeroPadded(std::size_t value, int width) {
  std::ostringstream stream;
  stream << std::setw(width) << std::setfill('0') << value;
  return stream.str();
}

[[nodiscard]] std::size_t ComputePrecomputeStepCount(const Constants& constants) {
  const double growth_interval_earth =
      constants.precompute_end_mass_earth - constants.initial_planet_mass_earth;
  if (growth_interval_earth <= 0.0) {
    throw std::runtime_error("Precompute growth interval must be positive.");
  }
  if (constants.precompute_step_mass_earth <= 0.0) {
    throw std::runtime_error("precompute_step_mass_earth must be positive.");
  }

  const double raw_step_count =
      growth_interval_earth / constants.precompute_step_mass_earth;
  const auto step_count = static_cast<std::size_t>(std::llround(raw_step_count));
  const double reconstructed_growth =
      static_cast<double>(step_count) * constants.precompute_step_mass_earth;
  if (std::fabs(reconstructed_growth - growth_interval_earth) > 1.0e-12) {
    throw std::runtime_error(
        "Precompute step mass does not evenly tile 0.01 -> 0.10 Earth masses.");
  }
  return step_count;
}

[[nodiscard]] Impactor MakePrecomputeImpactor(std::size_t step_index,
                                              const Constants& constants) {
  Impactor impactor{};
  impactor.label = "PRE-" + ZeroPadded(step_index, 5);
  impactor.phase = GrowthPhase::kPrecomputeAccretion;
  impactor.impactor_class = ImpactorClass::kECLikeUndifferentiated;
  impactor.mass_kg = EarthMassToKg(constants.precompute_step_mass_earth, constants);
  impactor.metallic_fe_fraction = constants.representative_metallic_fe_fraction;
  impactor.is_differentiated = false;
  impactor.silicate_inventory_source = InventorySource::kEcBulk;
  impactor.metal_inventory_source = InventorySource::kNone;
  impactor.silicate_volatiles_ppm = constants.ec_bulk_volatiles_ppm;
  impactor.metal_volatiles_ppm = {};
  impactor.delta15n_silicate_permil = constants.ec_delta15n_permil;
  impactor.delta15n_metal_permil = std::numeric_limits<double>::quiet_NaN();
  impactor.alloy = constants.representative_ec_alloy;
  impactor.Validate();
  return impactor;
}

[[nodiscard]] Impactor MakeEcDifferentiatedGi(std::size_t gi_number,
                                              double impactor_mass_earth,
                                              const Constants& constants) {
  Impactor impactor{};
  impactor.label = "GI" + std::to_string(gi_number);
  impactor.phase = GrowthPhase::kGiantImpact;
  impactor.impactor_class = ImpactorClass::kECLikeDifferentiated;
  impactor.mass_kg = EarthMassToKg(impactor_mass_earth, constants);
  impactor.metallic_fe_fraction = constants.representative_metallic_fe_fraction;
  impactor.is_differentiated = true;
  impactor.silicate_inventory_source =
      InventorySource::kPrecomputeOutputAt0p10EarthMass;
  impactor.metal_inventory_source = InventorySource::kPrecomputeCoreAt0p10EarthMass;
  impactor.silicate_volatiles_ppm = {};
  impactor.metal_volatiles_ppm = {};
  impactor.delta15n_silicate_permil = constants.ec_delta15n_permil;
  impactor.delta15n_metal_permil = std::numeric_limits<double>::quiet_NaN();
  impactor.alloy = constants.representative_ec_alloy;
  impactor.Validate();
  return impactor;
}

[[nodiscard]] Impactor MakeGi8CcEvent(const Constants& constants) {
  Impactor impactor{};
  impactor.label = "GI8";
  impactor.phase = GrowthPhase::kGiantImpact;
  impactor.impactor_class = ImpactorClass::kCILikeDifferentiated;
  impactor.mass_kg = EarthMassToKg(constants.standard_gi_mass_earth, constants);
  impactor.metallic_fe_fraction = constants.representative_metallic_fe_fraction;
  impactor.is_differentiated = true;
  impactor.silicate_inventory_source = InventorySource::kCiBulk;
  impactor.metal_inventory_source = InventorySource::kPrecomputeCoreAt0p10EarthMass;
  impactor.silicate_volatiles_ppm = constants.ci_bulk_volatiles_ppm;
  impactor.metal_volatiles_ppm = {};
  impactor.delta15n_silicate_permil = constants.ci_delta15n_permil;
  impactor.delta15n_metal_permil = std::numeric_limits<double>::quiet_NaN();
  impactor.alloy = constants.representative_ec_alloy;
  impactor.Validate();
  return impactor;
}

[[nodiscard]] Impactor MakeTheiaGi9(const Constants& constants) {
  Impactor impactor = MakeEcDifferentiatedGi(9, constants.theia_mass_earth, constants);
  impactor.label = "GI9-Theia";
  impactor.Validate();
  return impactor;
}

[[nodiscard]] Impactor MakeLateVeneerImpactor(const Constants& constants) {
  Impactor impactor{};
  impactor.label = "LateVeneer";
  impactor.phase = GrowthPhase::kLateVeneer;
  impactor.impactor_class = ImpactorClass::kECLikeUndifferentiated;
  impactor.mass_kg = EarthMassToKg(constants.late_veneer_mass_earth, constants);
  impactor.metallic_fe_fraction = constants.representative_metallic_fe_fraction;
  impactor.is_differentiated = false;
  impactor.silicate_inventory_source = InventorySource::kEcBulk;
  impactor.metal_inventory_source = InventorySource::kNone;
  impactor.silicate_volatiles_ppm = constants.ec_bulk_volatiles_ppm;
  impactor.metal_volatiles_ppm = {};
  impactor.delta15n_silicate_permil = constants.ec_delta15n_permil;
  impactor.delta15n_metal_permil = std::numeric_limits<double>::quiet_NaN();
  impactor.alloy = constants.representative_ec_alloy;
  impactor.Validate();
  return impactor;
}

[[nodiscard]] std::vector<Impactor> BuildScenario(const Constants& constants) {
  const std::size_t precompute_steps = ComputePrecomputeStepCount(constants);

  std::vector<Impactor> scenario;
  scenario.reserve(precompute_steps + 10);

  for (std::size_t step = 1; step <= precompute_steps; ++step) {
    scenario.push_back(MakePrecomputeImpactor(step, constants));
  }

  for (std::size_t gi_number = 1; gi_number <= 7; ++gi_number) {
    scenario.push_back(
        MakeEcDifferentiatedGi(gi_number, constants.standard_gi_mass_earth, constants));
  }

  scenario.push_back(MakeGi8CcEvent(constants));
  scenario.push_back(MakeTheiaGi9(constants));
  scenario.push_back(MakeLateVeneerImpactor(constants));

  double total_added_mass_earth = 0.0;
  for (const Impactor& impactor : scenario) {
    impactor.Validate();
    total_added_mass_earth += KgToEarthMass(impactor.mass_kg, constants);
  }

  const double expected_added_mass_earth =
      constants.final_planet_mass_earth - constants.initial_planet_mass_earth;
  if (std::fabs(total_added_mass_earth - expected_added_mass_earth) > 1.0e-12) {
    throw std::runtime_error("Scenario masses do not sum to the expected growth.");
  }

  return scenario;
}

[[nodiscard]] double ComputeScenarioAddedMassEarth(
    const std::vector<Impactor>& scenario, const Constants& constants) {
  double added_mass_earth = 0.0;
  for (const Impactor& impactor : scenario) {
    added_mass_earth += KgToEarthMass(impactor.mass_kg, constants);
  }
  return added_mass_earth;
}

[[nodiscard]] std::size_t CountPhase(const std::vector<Impactor>& scenario,
                                     GrowthPhase phase) {
  std::size_t count = 0;
  for (const Impactor& impactor : scenario) {
    if (impactor.phase == phase) {
      ++count;
    }
  }
  return count;
}

[[nodiscard]] std::string DescribeImpactor(const Constants& constants,
                                           const Impactor& impactor) {
  std::ostringstream stream;
  stream << std::fixed << std::setprecision(6);
  stream << impactor.label << " phase=" << ToString(impactor.phase)
         << " class=" << ToString(impactor.impactor_class)
         << " mass_earth=" << KgToEarthMass(impactor.mass_kg, constants)
         << " differentiated=" << (impactor.is_differentiated ? "yes" : "no")
         << " silicate_source=" << ToString(impactor.silicate_inventory_source)
         << " metal_source=" << ToString(impactor.metal_inventory_source);
  return stream.str();
}

void PrintInitialSummary(const Constants& constants, const PlanetState& planet) {
  std::cout << std::fixed << std::setprecision(6);
  const double radius_m = RadiusFromMassKg(planet.total_mass_kg(), constants);
  const double gravity_m_s2 = Gravity(planet.total_mass_kg(), radius_m, constants);
  const double escape_velocity_m_s =
      EscapeVelocity(planet.total_mass_kg(), radius_m, constants);
  const double mantle_thickness_m = MantleThickness(planet, constants);
  std::cout << "Scaffold ready\n";
  std::cout << "initial_planet_mass_earth="
            << KgToEarthMass(planet.total_mass_kg(), constants) << '\n';
  std::cout << "mantle_mass_kg=" << planet.mantle.mass_kg << '\n';
  std::cout << "core_mass_kg=" << planet.core.mass_kg << '\n';
  std::cout << "planet_radius_m=" << radius_m << '\n';
  std::cout << "surface_gravity_m_s2=" << gravity_m_s2 << '\n';
  std::cout << "escape_velocity_m_s=" << escape_velocity_m_s << '\n';
  std::cout << "mantle_thickness_m=" << mantle_thickness_m << '\n';
  std::cout << "mantle_C_ppm=" << planet.mantle.volatiles_ppm.carbon_ppm << '\n';
  std::cout << "mantle_H_ppm=" << planet.mantle.volatiles_ppm.hydrogen_ppm
            << '\n';
  std::cout << "mantle_N_ppm=" << planet.mantle.volatiles_ppm.nitrogen_ppm
            << '\n';
  std::cout << "mantle_delta15N_permil=" << planet.mantle.delta15n_permil
            << '\n';
  std::cout << "representative_alloy_xFe="
            << constants.representative_ec_alloy.x_fe() << '\n';
}

void PrintScenarioSummary(const Constants& constants, const PlanetState& initial_planet,
                          const std::vector<Impactor>& scenario) {
  const std::size_t precompute_count =
      CountPhase(scenario, GrowthPhase::kPrecomputeAccretion);
  const std::size_t gi_count = CountPhase(scenario, GrowthPhase::kGiantImpact);
  const std::size_t late_veneer_count =
      CountPhase(scenario, GrowthPhase::kLateVeneer);
  const double added_mass_earth = ComputeScenarioAddedMassEarth(scenario, constants);
  const double final_mass_earth =
      KgToEarthMass(initial_planet.total_mass_kg(), constants) + added_mass_earth;

  std::cout << "scenario_total_steps=" << scenario.size() << '\n';
  std::cout << "scenario_precompute_steps=" << precompute_count << '\n';
  std::cout << "scenario_gi_steps=" << gi_count << '\n';
  std::cout << "scenario_late_veneer_steps=" << late_veneer_count << '\n';
  std::cout << "scenario_added_mass_earth=" << added_mass_earth << '\n';
  std::cout << "scenario_final_mass_earth=" << final_mass_earth << '\n';

  if (scenario.empty()) {
    return;
  }

  std::cout << "scenario_first_precompute="
            << DescribeImpactor(constants, scenario.front()) << '\n';
  std::cout << "scenario_last_precompute="
            << DescribeImpactor(constants, scenario.at(precompute_count - 1)) << '\n';

  for (std::size_t index = precompute_count; index < scenario.size(); ++index) {
    std::cout << "scenario_major_step="
              << DescribeImpactor(constants, scenario.at(index)) << '\n';
  }
}

void PrintMeltPoolSummary(const Impactor& impactor, const MeltPoolState& melt_pool) {
  std::cout << "melt_pool_demo_impactor=" << impactor.label << '\n';
  std::cout << "melt_pool_iterations=" << melt_pool.iteration_count << '\n';
  std::cout << "melt_pool_impactor_density_kg_m3="
            << melt_pool.impactor_density_kg_m3 << '\n';
  std::cout << "melt_pool_impactor_diameter_m="
            << melt_pool.impactor_diameter_m << '\n';
  std::cout << "melt_pool_impact_velocity_m_s="
            << melt_pool.impact_velocity_m_s << '\n';
  std::cout << "melt_pool_target_volume_m3="
            << melt_pool.target_melt_volume_m3 << '\n';
  std::cout << "melt_pool_target_mass_kg="
            << melt_pool.target_melt_mass_kg << '\n';
  std::cout << "melt_pool_raw_depth_m=" << melt_pool.raw_depth_m << '\n';
  std::cout << "melt_pool_clamped_depth_m=" << melt_pool.depth_m << '\n';
  std::cout << "melt_pool_pressure_gpa=" << melt_pool.pressure_gpa << '\n';
  std::cout << "melt_pool_temperature_k=" << melt_pool.temperature_k << '\n';
  std::cout << "melt_pool_equilibrating_silicate_mass_kg="
            << melt_pool.equilibrating_silicate_mass_kg << '\n';
}

void PrintSimplifiedRunSummary(const Constants& constants,
                               const SimplifiedRunResult& run) {
  std::cout << "simple_run_completed_steps=" << run.completed_steps << '\n';
  std::cout << "simple_run_has_precompute_snapshot="
            << (run.has_precompute_snapshot ? "yes" : "no") << '\n';
  std::cout << "simple_run_final_planet_mass_earth="
            << KgToEarthMass(run.final_planet.total_mass_kg(), constants) << '\n';
  std::cout << "simple_run_final_mantle_mass_earth="
            << KgToEarthMass(run.final_planet.mantle.mass_kg, constants) << '\n';
  std::cout << "simple_run_final_core_mass_earth="
            << KgToEarthMass(run.final_planet.core.mass_kg, constants) << '\n';
  std::cout << "simple_run_final_mantle_C_ppm="
            << run.final_planet.mantle.volatiles_ppm.carbon_ppm << '\n';
  std::cout << "simple_run_final_mantle_H_ppm="
            << run.final_planet.mantle.volatiles_ppm.hydrogen_ppm << '\n';
  std::cout << "simple_run_final_mantle_N_ppm="
            << run.final_planet.mantle.volatiles_ppm.nitrogen_ppm << '\n';
  std::cout << "simple_run_final_core_C_ppm="
            << run.final_planet.core.volatiles_ppm.carbon_ppm << '\n';
  std::cout << "simple_run_final_core_H_ppm="
            << run.final_planet.core.volatiles_ppm.hydrogen_ppm << '\n';
  std::cout << "simple_run_final_core_N_ppm="
            << run.final_planet.core.volatiles_ppm.nitrogen_ppm << '\n';
  std::cout << "simple_run_final_atm_C_kg="
            << run.final_planet.atmosphere.elements_kg.carbon_kg << '\n';
  std::cout << "simple_run_final_atm_H_kg="
            << run.final_planet.atmosphere.elements_kg.hydrogen_kg << '\n';
  std::cout << "simple_run_final_atm_N_kg="
            << run.final_planet.atmosphere.elements_kg.nitrogen_kg << '\n';
  std::cout << "simple_run_final_atm_species_C="
            << ToString(run.final_planet.atmosphere.carbon_species) << '\n';
  std::cout << "simple_run_final_atm_species_H="
            << ToString(run.final_planet.atmosphere.hydrogen_species) << '\n';
  std::cout << "simple_run_final_atm_species_N="
            << ToString(run.final_planet.atmosphere.nitrogen_species) << '\n';
  std::cout << "simple_run_cumulative_eroded_C_kg="
            << run.cumulative_eroded_kg.carbon_kg << '\n';
  std::cout << "simple_run_cumulative_eroded_H_kg="
            << run.cumulative_eroded_kg.hydrogen_kg << '\n';
  std::cout << "simple_run_cumulative_eroded_N_kg="
            << run.cumulative_eroded_kg.nitrogen_kg << '\n';
  std::cout << "simple_run_max_abs_conservation_error_C_kg="
            << run.max_abs_conservation_error_kg.carbon_kg << '\n';
  std::cout << "simple_run_max_abs_conservation_error_H_kg="
            << run.max_abs_conservation_error_kg.hydrogen_kg << '\n';
  std::cout << "simple_run_max_abs_conservation_error_N_kg="
            << run.max_abs_conservation_error_kg.nitrogen_kg << '\n';
  std::cout << "simple_run_first_step_error_max_kg="
            << MaxComponentAbs(run.first_step.conservation_error_kg) << '\n';
  std::cout << "simple_run_last_step_error_max_kg="
            << MaxComponentAbs(run.last_step.conservation_error_kg) << '\n';
}

}  // namespace self_oxidation

int main() {
  try {
    const self_oxidation::Constants constants{};
    constants.representative_ec_alloy.Validate();

    const self_oxidation::PlanetState initial_planet =
        self_oxidation::MakeInitialPlanetState(constants);
    const std::vector<self_oxidation::Impactor> scenario =
        self_oxidation::BuildScenario(constants);
    const self_oxidation::SimplifiedRunResult simple_run =
        self_oxidation::RunScenarioSimplified(initial_planet, scenario, constants);

    self_oxidation::PrintInitialSummary(constants, initial_planet);
    self_oxidation::PrintScenarioSummary(constants, initial_planet, scenario);
    self_oxidation::PrintMeltPoolSummary(simple_run.first_step.resolved_impactor,
                                         simple_run.first_step.melt_pool);
    self_oxidation::PrintSimplifiedRunSummary(constants, simple_run);
    return 0;
  } catch (const std::exception& exception) {
    std::cerr << "error: " << exception.what() << '\n';
    return 1;
  }
}
