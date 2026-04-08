#include <algorithm>
#include <cmath>
#include <fstream>
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

enum class RedoxSource {
  kImpactorValue,
  kPrecomputeOutputAt0p10EarthMass,
  kCurrentTargetMantle,
};

enum class ErosionSizeDistribution {
  kDMinus2,
  kDMinus3,
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

  void Validate() const {
    if (fe3_over_fe_total < 0.0 || fe3_over_fe_total >= 1.0) {
      throw std::runtime_error("fe3_over_fe_total must be in [0, 1). ");
    }
    if (x_feo_silicate < 0.0 || x_feo_silicate > 1.0) {
      throw std::runtime_error("x_feo_silicate must be in [0, 1].");
    }
  }
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
    redox.Validate();
  }
};

struct Impactor {
  std::string label;
  GrowthPhase phase = GrowthPhase::kPrecomputeAccretion;
  ImpactorClass impactor_class = ImpactorClass::kECLikeUndifferentiated;
  double mass_kg = 0.0;
  double metallic_fe_fraction = 0.325;
  bool is_differentiated = false;
  InventorySource bulk_inventory_source = InventorySource::kNone;
  InventorySource silicate_inventory_source = InventorySource::kNone;
  InventorySource metal_inventory_source = InventorySource::kNone;
  InventorySource alloy_source = InventorySource::kNone;
  ElementPpm bulk_volatiles_ppm{};
  ElementPpm silicate_volatiles_ppm{};
  ElementPpm metal_volatiles_ppm{};
  double delta15n_silicate_permil = 0.0;
  double delta15n_metal_permil = 0.0;
  AlloyComposition alloy{};
  RedoxSource redox_source = RedoxSource::kImpactorValue;
  RedoxState redox{};

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
    const bool has_bulk_inventory =
        bulk_inventory_source != InventorySource::kNone;
    const bool has_phase_resolved_inventory =
        silicate_inventory_source != InventorySource::kNone ||
        metal_inventory_source != InventorySource::kNone;
    if (has_bulk_inventory && has_phase_resolved_inventory) {
      throw std::runtime_error(
          "Impactor cannot mix bulk and phase-resolved volatile inventories.");
    }
    alloy.Validate();
    redox.Validate();
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

struct HistoryRecord {
  std::size_t step_index = 0;
  std::string label;
  GrowthPhase phase = GrowthPhase::kPrecomputeAccretion;
  ImpactorClass impactor_class = ImpactorClass::kECLikeUndifferentiated;
  double planet_mass_before_earth = 0.0;
  double planet_mass_after_earth = 0.0;
  double mantle_mass_after_earth = 0.0;
  double core_mass_after_earth = 0.0;
  double bse_carbon_ppm = 0.0;
  double bse_hydrogen_ppm = 0.0;
  double bse_nitrogen_ppm = 0.0;
  double mantle_carbon_ppm = 0.0;
  double mantle_hydrogen_ppm = 0.0;
  double mantle_nitrogen_ppm = 0.0;
  double core_carbon_ppm = 0.0;
  double core_hydrogen_ppm = 0.0;
  double core_nitrogen_ppm = 0.0;
  double atmosphere_carbon_ppm_silicate_earth = 0.0;
  double atmosphere_hydrogen_ppm_silicate_earth = 0.0;
  double atmosphere_nitrogen_ppm_silicate_earth = 0.0;
  double mantle_delta15n_permil = 0.0;
  double core_delta15n_permil = 0.0;
  double atmosphere_delta15n_permil = 0.0;
};

struct DifferentiatedImpactorSnapshot {
  ElementPpm silicate_volatiles_ppm{};
  ElementPpm metal_volatiles_ppm{};
  double delta15n_silicate_permil = 0.0;
  double delta15n_metal_permil = 0.0;
  AlloyComposition alloy{};
  RedoxState redox{};
};

struct SimplifiedRunResult {
  PlanetState final_planet{};
  std::size_t completed_steps = 0;
  bool has_precompute_snapshot = false;
  DifferentiatedImpactorSnapshot precompute_snapshot{};
  std::vector<HistoryRecord> history{};
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

struct ImpactorVolatileMassesKg {
  ElementMassKg silicate_kg{};
  ElementMassKg metal_kg{};
};

struct Step2EquilibriumState {
  ElementPartitionKg carbon{};
  ElementPartitionKg hydrogen{};
  ElementPartitionKg nitrogen{};
  double mean_molar_mass_kg_per_mol = 0.0;
  double carbon_partial_pressure_pa = 0.0;
  double hydrogen_partial_pressure_pa = 0.0;
  double nitrogen_partial_pressure_pa = 0.0;
  int outer_iteration_count = 0;
};

struct ErosionEfficiencies {
  double direct_atmosphere_efficiency = 0.0;
  double impactor_vapor_escape_efficiency = 0.0;
};

struct AtmosphericErosionOutcome {
  double erosion_fraction = 0.0;
  double impactor_vapor_escape_efficiency = 0.0;
  ElementMassKg eroded_elements_kg{};
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

  double fixed_nbo_over_t = 2.5;

  double step2_outer_rel_tolerance = 1.0e-4;
  int step2_outer_max_iterations = 32;
  double step2_root_rel_tolerance = 1.0e-10;
  int step2_root_max_iterations = 128;

  double carbon_solubility_oxidized_ppm = 1.6;
  double nitrogen_solubility_oxidized_ppm = 1.0;
  double carbon_solubility_intermediate_ppm = 0.55;
  double nitrogen_solubility_intermediate_ppm = 5.0;
  double carbon_solubility_reduced_ppm = 0.22;
  double nitrogen_solubility_reduced_ppm = 50.0;
  double hydrogen_h2_solubility_ppm_sqrt_mpa = 5.0;

  double moore_a = 2565.0;
  double moore_b_al = -1.997;
  double moore_b_fe = -0.9275;
  double moore_b_na = 2.736;
  double moore_x_al = 0.137;
  double moore_x_fe = 0.124;
  double moore_x_na = 0.0268;
  double moore_c = 1.171;
  double moore_d = -14.21;
  double moore_r_c_to_x = 0.309;

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

  double oxygen_log10_kd_intercept = 0.6;
  double oxygen_log10_kd_temperature_coeff = -3800.0;
  double oxygen_log10_kd_pressure_coeff = 0.0;

  double dc_intercept = 1.49;
  double dc_temperature_coeff = 3000.0;
  double dc_pressure_coeff = -235.0;
  double dc_xs_coeff = 9.6;
  double dc_xo_coeff = -19.5;
  double dc_nbo_coeff = -0.118;
  double dc_delta_iw_coeff = -0.238;
  double erosion_magma_ocean_temperature_k = 1500.0;
  double erosion_late_veneer_temperature_k = 288.0;
  double erosion_target_surface_density_kg_m3 = 2630.0;
  double erosion_impactor_density_kg_m3 = 3320.0;
  double erosion_integration_rel_tolerance = 1.0e-2;
  int erosion_integration_max_refinements = 10;

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

[[nodiscard]] ElementMassKg ComputeMantleVolatileMassKg(const PlanetState& planet) {
  return PpmToElementMassKg(planet.mantle.volatiles_ppm, planet.mantle.mass_kg);
}

[[nodiscard]] ElementMassKg ComputeCoreVolatileMassKg(const PlanetState& planet) {
  return PpmToElementMassKg(planet.core.volatiles_ppm, planet.core.mass_kg);
}

[[nodiscard]] ElementMassKg ComputeBseVolatileMassKg(const PlanetState& planet) {
  const ElementMassKg mantle = ComputeMantleVolatileMassKg(planet);
  return {
      mantle.carbon_kg + planet.atmosphere.elements_kg.carbon_kg,
      mantle.hydrogen_kg + planet.atmosphere.elements_kg.hydrogen_kg,
      mantle.nitrogen_kg + planet.atmosphere.elements_kg.nitrogen_kg,
  };
}

[[nodiscard]] ElementPpm ComputeBseVolatilesPpmPerPlanetMass(
    const PlanetState& planet) {
  return ElementMassKgToPpm(ComputeBseVolatileMassKg(planet), planet.total_mass_kg());
}

[[nodiscard]] ElementPpm ComputeAtmosphereVolatilesPpmPerSilicateEarth(
    const PlanetState& planet) {
  return ElementMassKgToPpm(planet.atmosphere.elements_kg, planet.mantle.mass_kg);
}

[[nodiscard]] HistoryRecord CaptureHistoryRecord(const std::size_t index,
                                                 const PlanetState& planet_before,
                                                 const StepResult& step_result,
                                                 const Constants& constants) {
  const PlanetState& planet_after = step_result.planet_after;
  const ElementPpm bse_ppm = ComputeBseVolatilesPpmPerPlanetMass(planet_after);
  const ElementPpm atmosphere_ppm =
      ComputeAtmosphereVolatilesPpmPerSilicateEarth(planet_after);

  HistoryRecord record{};
  record.step_index = index;
  record.label = step_result.resolved_impactor.label;
  record.phase = step_result.resolved_impactor.phase;
  record.impactor_class = step_result.resolved_impactor.impactor_class;
  record.planet_mass_before_earth =
      KgToEarthMass(planet_before.total_mass_kg(), constants);
  record.planet_mass_after_earth =
      KgToEarthMass(planet_after.total_mass_kg(), constants);
  record.mantle_mass_after_earth =
      KgToEarthMass(planet_after.mantle.mass_kg, constants);
  record.core_mass_after_earth =
      KgToEarthMass(planet_after.core.mass_kg, constants);
  record.bse_carbon_ppm = bse_ppm.carbon_ppm;
  record.bse_hydrogen_ppm = bse_ppm.hydrogen_ppm;
  record.bse_nitrogen_ppm = bse_ppm.nitrogen_ppm;
  record.mantle_carbon_ppm = planet_after.mantle.volatiles_ppm.carbon_ppm;
  record.mantle_hydrogen_ppm = planet_after.mantle.volatiles_ppm.hydrogen_ppm;
  record.mantle_nitrogen_ppm = planet_after.mantle.volatiles_ppm.nitrogen_ppm;
  record.core_carbon_ppm = planet_after.core.volatiles_ppm.carbon_ppm;
  record.core_hydrogen_ppm = planet_after.core.volatiles_ppm.hydrogen_ppm;
  record.core_nitrogen_ppm = planet_after.core.volatiles_ppm.nitrogen_ppm;
  record.atmosphere_carbon_ppm_silicate_earth = atmosphere_ppm.carbon_ppm;
  record.atmosphere_hydrogen_ppm_silicate_earth = atmosphere_ppm.hydrogen_ppm;
  record.atmosphere_nitrogen_ppm_silicate_earth = atmosphere_ppm.nitrogen_ppm;
  record.mantle_delta15n_permil = planet_after.mantle.delta15n_permil;
  record.core_delta15n_permil = planet_after.core.delta15n_permil;
  record.atmosphere_delta15n_permil = planet_after.atmosphere.delta15n_permil;
  return record;
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

  const double log10_dn =
      constants.dn_intercept + constants.dn_temperature_coeff / temperature_k +
      constants.dn_pressure_coeff * pressure_gpa / temperature_k +
      constants.dn_nbo_coeff * constants.fixed_nbo_over_t +
      constants.dn_delta_iw_coeff * delta_iw_eq +
      constants.dn_xs_coeff *
          SafeLog10(SafeOneMinus(alloy.x_s, "x_s^metal"), "1 - x_s^metal") +
      constants.dn_xc_coeff *
          SafeLog10(SafeOneMinus(alloy.x_c, "x_c^metal"), "1 - x_c^metal") +
      constants.dn_xni_coeff *
          SafeLog10(SafeOneMinus(alloy.x_ni, "x_ni^metal"), "1 - x_ni^metal") +
      constants.dn_xsi_coeff *
          SafeLog10(SafeOneMinus(alloy.x_si, "x_si^metal"), "1 - x_si^metal");
  return std::pow(10.0, log10_dn);
}

[[nodiscard]] double ComputeNitrogenIsotopeFractionationPermil(double temperature_k,
                                                               double delta_iw_eq,
                                                               const Constants& constants) {
  if (temperature_k <= 0.0) {
    throw std::runtime_error(
        "ComputeNitrogenIsotopeFractionationPermil requires T > 0.");
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

  const double log10_dh =
      constants.dh_intercept + constants.dh_temperature_coeff / temperature_k +
      constants.dh_pressure_coeff * pressure_gpa / temperature_k +
      constants.dh_xc_coeff *
          SafeLog10(SafeOneMinus(alloy.x_c, "x_c^metal"), "1 - x_c^metal") -
      0.25 * delta_iw_eq;
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
    throw std::runtime_error(
        "ComputeOxygenMoleFractionInMetal requires X_FeO >= 0.");
  }

  const double x_fe = alloy.x_fe();
  const double log10_kd_o =
      constants.oxygen_log10_kd_intercept +
      constants.oxygen_log10_kd_temperature_coeff / temperature_k +
      constants.oxygen_log10_kd_pressure_coeff * pressure_gpa / temperature_k;
  const double kd_o = std::pow(10.0, log10_kd_o);
  const double raw_x_o =
      kd_o * x_feo_silicate / std::max(x_fe, 1.0e-12);
  if (!std::isfinite(raw_x_o)) {
    throw std::runtime_error(
        "ComputeOxygenMoleFractionInMetal produced a non-finite X_O.");
  }

  if (!(raw_x_o < 1.0)) {
    std::ostringstream message;
    message << "ComputeOxygenMoleFractionInMetal produced X_O >= 1.0"
            << "; X_O=" << raw_x_o
            << ", P=" << pressure_gpa
            << " GPa, T=" << temperature_k
            << " K, X_FeO^sil=" << x_feo_silicate
            << ", X_Fe^metal=" << x_fe;
    throw std::runtime_error(message.str());
  }
  return raw_x_o;
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
    throw std::runtime_error(
        "ComputeCarbonPartitionCoefficient requires non-negative X_S and X_O.");
  }

  const double log10_dc =
      constants.dc_intercept + constants.dc_temperature_coeff / temperature_k +
      constants.dc_pressure_coeff * pressure_gpa / temperature_k +
      constants.dc_xs_coeff *
          SafeLog10(SafeOneMinus(x_s_metal, "X_S^metal"), "1 - X_S^metal") +
      constants.dc_xo_coeff *
          SafeLog10(SafeOneMinus(x_o_metal, "X_O^metal"), "1 - X_O^metal") +
      constants.dc_nbo_coeff * constants.fixed_nbo_over_t +
      constants.dc_delta_iw_coeff * delta_iw_eq;
  return std::pow(10.0, log10_dc);
}

[[nodiscard]] AlloyComposition UpdatePrecomputeAlloyFromCoreCarbonPpm(
    double core_carbon_ppm, const AlloyComposition& reference_alloy) {
  reference_alloy.Validate();
  if (core_carbon_ppm < 0.0) {
    throw std::runtime_error(
        "UpdatePrecomputeAlloyFromCoreCarbonPpm requires non-negative C ppm.");
  }

  constexpr double kAtomicMassCarbon = 12.011;
  constexpr double kAtomicMassSulfur = 32.06;
  constexpr double kAtomicMassNickel = 58.6934;
  constexpr double kAtomicMassSilicon = 28.085;
  constexpr double kAtomicMassIron = 55.845;

  const double fixed_minor_sum =
      reference_alloy.x_s + reference_alloy.x_ni + reference_alloy.x_si;
  const double max_x_c = std::max(0.0, 1.0 - fixed_minor_sum - 1.0e-12);
  const double carbon_mass_fraction = PpmToMassFraction(core_carbon_ppm);
  if (carbon_mass_fraction < 0.0 || carbon_mass_fraction > 1.0) {
    throw std::runtime_error(
        "Precompute core carbon ppm implies an invalid carbon mass fraction.");
  }

  const double non_carbon_molar_mass =
      reference_alloy.x_s * kAtomicMassSulfur +
      reference_alloy.x_ni * kAtomicMassNickel +
      reference_alloy.x_si * kAtomicMassSilicon +
      (1.0 - fixed_minor_sum) * kAtomicMassIron;
  const double denominator =
      kAtomicMassCarbon -
      carbon_mass_fraction * (kAtomicMassCarbon - kAtomicMassIron);
  if (!(denominator > 0.0)) {
    throw std::runtime_error(
        "Carbon mole-fraction denominator became non-positive.");
  }

  const double raw_x_c =
      carbon_mass_fraction * non_carbon_molar_mass / denominator;
  if (raw_x_c < 0.0 || raw_x_c > max_x_c) {
    throw std::runtime_error(
        "Precompute core carbon ppm implies an invalid x_C^metal.");
  }

  AlloyComposition updated = reference_alloy;
  updated.x_c = raw_x_c;
  updated.Validate();
  return updated;
}

struct SilicateIronInventory {
  double non_fe_cation_mol = 0.0;
  double feo_mol = 0.0;
  double feo1p5_mol = 0.0;

  [[nodiscard]] double total_cation_mol() const {
    return non_fe_cation_mol + feo_mol + feo1p5_mol;
  }

  [[nodiscard]] double total_fe_mol() const {
    return feo_mol + feo1p5_mol;
  }

  [[nodiscard]] double ferric_fraction() const {
    const double total_fe = total_fe_mol();
    if (total_fe <= 0.0) {
      return 0.0;
    }
    return feo1p5_mol / total_fe;
  }
};

struct OxideCellEosParams {
  double v0_a3 = 0.0;
  double k0_gpa = 0.0;
  double k0_prime = 0.0;
  double k0_second_gpa_inv = 0.0;
  double thermal_a = 0.0;
  double thermal_b = 0.0;
  double thermal_c = 0.0;
};

struct IwPolynomial {
  double m0 = 0.0;
  double m1 = 0.0;
  double m2 = 0.0;
  double m3 = 0.0;
  double m4 = 0.0;
};

constexpr double kUniversalGasConstantJPerMolK = 8.31446261815324;
constexpr double kLog10Factor = 2.302585092994046;
constexpr double kReferencePressureGpa = 1.0e-4;
constexpr double kFerricReferenceTemperatureK = 1673.15;
constexpr double kFerricDeltaCpJPerMolK = 33.25;
constexpr double kFerricFitA = 0.1917;
constexpr double kFerricFitB = -1.961;
constexpr double kFerricFitC = 4158.1;
constexpr double kFerricY1 = -520.46;
constexpr double kFerricY2 = -185.37;
constexpr double kFerricY3 = 494.39;
constexpr double kFerricY4 = 1838.34;
constexpr double kFerricY5 = 2888.48;
constexpr double kFerricY6 = 3473.68;
constexpr double kFerricY7 = -4473.6;
constexpr double kFerricY8 = -1245.09;
constexpr double kFerricY9 = -1156.86;
constexpr double kSilicateXSio2 = 0.380;
constexpr double kSilicateXAlO1p5 = 0.045;
constexpr double kSilicateXMgO = 0.479;
constexpr double kSilicateXCaO = 0.031;
constexpr double kSilicateXFeOInitial = 0.056;
constexpr double kSilicateXTiO2 = 0.001;
constexpr double kSilicateXNaO0p5 = 0.007;
constexpr double kSilicateXKO0p5 = 0.0;
constexpr double kSilicateXPO2p5 = 0.0;
constexpr double kSilicateReferenceNonFeCationFractionSum =
    kSilicateXSio2 + kSilicateXAlO1p5 + kSilicateXMgO + kSilicateXCaO +
    kSilicateXTiO2 + kSilicateXNaO0p5 + kSilicateXKO0p5 + kSilicateXPO2p5;
constexpr double kMolarMassSiO2KgPerMol = 60.0843e-3;
constexpr double kMolarMassAlO1p5KgPerMol = 50.9800385e-3;
constexpr double kMolarMassMgOKgPerMol = 40.3044e-3;
constexpr double kMolarMassCaOKgPerMol = 56.0774e-3;
constexpr double kMolarMassFeOKgPerMol = 71.844e-3;
constexpr double kMolarMassFeO1p5KgPerMol = 79.8435e-3;
constexpr double kMolarMassTiO2KgPerMol = 79.8658e-3;
constexpr double kMolarMassNaO0p5KgPerMol = 30.98976928e-3;
constexpr double kMolarMassKO0p5KgPerMol = 47.0978e-3;
constexpr double kMolarMassPO2p5KgPerMol = 70.971761998e-3;
constexpr double kMolarMassFeKgPerMol = 55.845e-3;
constexpr double kMolarMassCarbonElementKgPerMol = 12.011e-3;
constexpr double kMolarMassHydrogenElementKgPerMol = 2.0 * 1.00794e-3;
constexpr double kMolarMassNitrogenElementKgPerMol = 28.0134e-3;
constexpr double kMolarMassH2OKgPerMol = 18.01528e-3;
constexpr double kAngstrom3ToCm3PerMol = 0.602214076;
constexpr double kSurfaceAdiabaticGradientKPerGpa = 15.0;
constexpr int kPressureIntegralSteps = 48;
constexpr int kVolumeSolveMaxIterations = 80;
constexpr double kVolumeSolveRelativeTolerance = 1.0e-10;
constexpr double kIwBoundaryA = -18.64;
constexpr double kIwBoundaryB = 0.04359;
constexpr double kIwBoundaryC = -5.069e-6;
constexpr IwPolynomial kIwCoeffLowA{6.844864, 1.175691e-1, 1.143873e-3, 0.0, 0.0};
constexpr IwPolynomial kIwCoeffLowB{5.791364e-4, -2.891434e-4, -2.737171e-7, 0.0, 0.0};
constexpr IwPolynomial kIwCoeffLowC{-7.971469e-5, 3.198005e-5, 0.0, 1.059554e-10, 2.014461e-7};
constexpr IwPolynomial kIwCoeffLowD{-2.769002e4, 5.285977e2, -2.919275, 0.0, 0.0};
constexpr IwPolynomial kIwCoeffHighA{8.463095, -3.000307e-3, 7.213445e-5, 0.0, 0.0};
constexpr IwPolynomial kIwCoeffHighB{1.148738e-3, -9.352312e-5, 5.161592e-7, 0.0, 0.0};
constexpr IwPolynomial kIwCoeffHighC{-7.448624e-4, -6.329325e-6, 0.0, -1.407339e-10, 1.830014e-4};
constexpr IwPolynomial kIwCoeffHighD{-2.782082e4, 5.285977e2, -8.473231e-1, 0.0, 0.0};
constexpr OxideCellEosParams kFerricCellEos{1204.69, 21.35, 3.626, 0.009, 34.53, 68.64, 35.27};
constexpr OxideCellEosParams kFerrousCellEos{1180.10, 24.22, 3.382, 0.012, 35.70, 71.10, 36.60};

[[nodiscard]] double ComputeFerricComponentMoleFraction(const RedoxState& redox) {
  redox.Validate();
  if (redox.fe3_over_fe_total <= 0.0) {
    return 0.0;
  }
  return redox.x_feo_silicate * redox.fe3_over_fe_total /
         SafeOneMinus(redox.fe3_over_fe_total, "fe3_over_fe_total");
}

[[nodiscard]] double ComputeTotalIronCationFraction(const RedoxState& redox) {
  return redox.x_feo_silicate + ComputeFerricComponentMoleFraction(redox);
}

[[nodiscard]] double ComputeNonFeScaleFactor(const RedoxState& redox) {
  const double remaining_cation_fraction =
      SafeOneMinus(ComputeTotalIronCationFraction(redox),
                   "X_Fe,total^silicate");
  return remaining_cation_fraction / kSilicateReferenceNonFeCationFractionSum;
}

[[nodiscard]] double ComputeSilicateMeanCationMolarMassKgPerMol(
    const RedoxState& redox) {
  const double ferric_fraction = ComputeFerricComponentMoleFraction(redox);
  const double non_fe_scale = ComputeNonFeScaleFactor(redox);
  return non_fe_scale *
             (kSilicateXSio2 * kMolarMassSiO2KgPerMol +
              kSilicateXAlO1p5 * kMolarMassAlO1p5KgPerMol +
              kSilicateXMgO * kMolarMassMgOKgPerMol +
              kSilicateXCaO * kMolarMassCaOKgPerMol +
              kSilicateXTiO2 * kMolarMassTiO2KgPerMol +
              kSilicateXNaO0p5 * kMolarMassNaO0p5KgPerMol +
              kSilicateXKO0p5 * kMolarMassKO0p5KgPerMol +
              kSilicateXPO2p5 * kMolarMassPO2p5KgPerMol) +
         redox.x_feo_silicate * kMolarMassFeOKgPerMol +
         ferric_fraction * kMolarMassFeO1p5KgPerMol;
}

[[nodiscard]] SilicateIronInventory ComputeSilicateIronInventory(
    double silicate_mass_kg, const RedoxState& redox) {
  if (silicate_mass_kg < 0.0) {
    throw std::runtime_error(
        "ComputeSilicateIronInventory requires non-negative mass.");
  }
  if (silicate_mass_kg == 0.0) {
    return {};
  }

  const double ferric_fraction = ComputeFerricComponentMoleFraction(redox);
  const double mean_molar_mass_kg_per_mol =
      ComputeSilicateMeanCationMolarMassKgPerMol(redox);
  const double total_cation_mol = silicate_mass_kg / mean_molar_mass_kg_per_mol;
  SilicateIronInventory inventory{};
  inventory.feo_mol = redox.x_feo_silicate * total_cation_mol;
  inventory.feo1p5_mol = ferric_fraction * total_cation_mol;
  inventory.non_fe_cation_mol =
      total_cation_mol - inventory.feo_mol - inventory.feo1p5_mol;
  return inventory;
}

[[nodiscard]] SilicateIronInventory AddSilicateIronInventory(
    const SilicateIronInventory& left, const SilicateIronInventory& right) {
  return {
      left.non_fe_cation_mol + right.non_fe_cation_mol,
      left.feo_mol + right.feo_mol,
      left.feo1p5_mol + right.feo1p5_mol,
  };
}

[[nodiscard]] RedoxState MakeRedoxStateFromIronInventory(
    const SilicateIronInventory& inventory, double delta_iw_eq,
    double delta_iw_surface, double oxidized_fe0_mass_kg) {
  RedoxState redox{};
  const double total_cation_mol = inventory.total_cation_mol();
  const double total_fe_mol = inventory.total_fe_mol();
  redox.fe3_over_fe_total =
      total_fe_mol > 0.0 ? inventory.feo1p5_mol / total_fe_mol : 0.0;
  redox.x_feo_silicate =
      total_cation_mol > 0.0 ? inventory.feo_mol / total_cation_mol : 0.0;
  redox.delta_iw_eq = delta_iw_eq;
  redox.delta_iw_surface = delta_iw_surface;
  redox.oxidized_fe0_mass_kg = oxidized_fe0_mass_kg;
  redox.Validate();
  return redox;
}

[[nodiscard]] double ComputeFerricHeatCapacityTermLog10(double temperature_k) {
  if (temperature_k <= 0.0) {
    throw std::runtime_error("ComputeFerricHeatCapacityTermLog10 requires T > 0.");
  }
  return kFerricDeltaCpJPerMolK / (kUniversalGasConstantJPerMolK * kLog10Factor) *
         (1.0 - kFerricReferenceTemperatureK / temperature_k -
          std::log(temperature_k / kFerricReferenceTemperatureK));
}

[[nodiscard]] double ComputeFerricCompositionTermLog10(double temperature_k) {
  if (temperature_k <= 0.0) {
    throw std::runtime_error("ComputeFerricCompositionTermLog10 requires T > 0.");
  }
  const double composition_sum =
      kFerricY1 * kSilicateXSio2 + kFerricY2 * kSilicateXTiO2 +
      kFerricY3 * kSilicateXMgO + kFerricY4 * kSilicateXCaO +
      kFerricY5 * kSilicateXNaO0p5 + kFerricY6 * kSilicateXKO0p5 +
      kFerricY7 * kSilicateXPO2p5 +
      kFerricY8 * kSilicateXSio2 * kSilicateXAlO1p5 +
      kFerricY9 * kSilicateXSio2 * kSilicateXMgO;
  return composition_sum / temperature_k;
}

[[nodiscard]] double EvaluateIwPolynomial(const IwPolynomial& polynomial,
                                          double pressure_gpa) {
  const double pressure_sqrt = std::sqrt(std::max(pressure_gpa, 0.0));
  return polynomial.m0 + polynomial.m1 * pressure_gpa +
         polynomial.m2 * pressure_gpa * pressure_gpa +
         polynomial.m3 * pressure_gpa * pressure_gpa * pressure_gpa +
         polynomial.m4 * pressure_sqrt;
}

[[nodiscard]] double ComputeIwBoundaryPressureGpa(double temperature_k) {
  return kIwBoundaryA + kIwBoundaryB * temperature_k +
         kIwBoundaryC * temperature_k * temperature_k;
}

[[nodiscard]] double ComputeLog10FugacityIw(double temperature_k,
                                            double pressure_gpa) {
  if (temperature_k <= 0.0) {
    throw std::runtime_error("ComputeLog10FugacityIw requires T > 0.");
  }
  if (pressure_gpa < 0.0) {
    throw std::runtime_error("ComputeLog10FugacityIw requires P >= 0.");
  }

  const bool use_hcp = pressure_gpa > ComputeIwBoundaryPressureGpa(temperature_k);
  const IwPolynomial& coeff_a = use_hcp ? kIwCoeffHighA : kIwCoeffLowA;
  const IwPolynomial& coeff_b = use_hcp ? kIwCoeffHighB : kIwCoeffLowB;
  const IwPolynomial& coeff_c = use_hcp ? kIwCoeffHighC : kIwCoeffLowC;
  const IwPolynomial& coeff_d = use_hcp ? kIwCoeffHighD : kIwCoeffLowD;
  const double a_value = EvaluateIwPolynomial(coeff_a, pressure_gpa);
  const double b_value = EvaluateIwPolynomial(coeff_b, pressure_gpa);
  const double c_value = EvaluateIwPolynomial(coeff_c, pressure_gpa);
  const double d_value = EvaluateIwPolynomial(coeff_d, pressure_gpa);
  return a_value + b_value * temperature_k +
         c_value * temperature_k * std::log(temperature_k) +
         d_value / temperature_k;
}

[[nodiscard]] double ComputeLog10FugacityFromDeltaIw(double delta_iw,
                                                     double temperature_k,
                                                     double pressure_gpa) {
  return delta_iw + ComputeLog10FugacityIw(temperature_k, pressure_gpa);
}

[[nodiscard]] double ComputeBirchMurnaghanPressureAtReferenceGpa(
    double volume_a3, const OxideCellEosParams& params) {
  if (volume_a3 <= 0.0) {
    throw std::runtime_error(
        "ComputeBirchMurnaghanPressureAtReferenceGpa requires V > 0.");
  }

  const double eta2 = std::pow(params.v0_a3 / volume_a3, 2.0 / 3.0);
  const double finite_strain = 0.5 * (eta2 - 1.0);
  const double correction =
      1.0 + 1.5 * (params.k0_prime - 4.0) * finite_strain +
      1.5 * (params.k0_gpa * params.k0_second_gpa_inv +
             (params.k0_prime - 4.0) * (params.k0_prime - 3.0) + 35.0 / 9.0) *
          finite_strain * finite_strain;
  return 3.0 * params.k0_gpa * finite_strain *
         std::pow(1.0 + 2.0 * finite_strain, 2.5) * correction;
}

[[nodiscard]] double ComputeThermalPressureCoefficientGpaPerK(
    double volume_a3, const OxideCellEosParams& params) {
  if (volume_a3 <= 0.0) {
    throw std::runtime_error(
        "ComputeThermalPressureCoefficientGpaPerK requires V > 0.");
  }
  const double ratio = volume_a3 / params.v0_a3;
  return (params.thermal_a - params.thermal_b * ratio +
          params.thermal_c * ratio * ratio) /
         1000.0;
}

[[nodiscard]] double ComputeOxideCellPressureGpa(double volume_a3,
                                                 double temperature_k,
                                                 const OxideCellEosParams& params) {
  return ComputeBirchMurnaghanPressureAtReferenceGpa(volume_a3, params) +
         ComputeThermalPressureCoefficientGpaPerK(volume_a3, params) *
             (temperature_k - 3000.0);
}

[[nodiscard]] double SolveOxideCellVolumeA3(double target_pressure_gpa,
                                            double temperature_k,
                                            const OxideCellEosParams& params) {
  if (target_pressure_gpa < 0.0) {
    throw std::runtime_error("SolveOxideCellVolumeA3 requires P >= 0.");
  }

  constexpr int kBracketScanSteps = 600;
  const double minimum_multiplier = 0.35;
  const double maximum_multiplier = 3.0;
  double previous_multiplier = minimum_multiplier;
  double previous_pressure_gpa = ComputeOxideCellPressureGpa(
      previous_multiplier * params.v0_a3, temperature_k, params);

  for (int step = 1; step <= kBracketScanSteps; ++step) {
    const double current_multiplier =
        minimum_multiplier +
        (maximum_multiplier - minimum_multiplier) *
            static_cast<double>(step) / static_cast<double>(kBracketScanSteps);
    const double current_pressure_gpa = ComputeOxideCellPressureGpa(
        current_multiplier * params.v0_a3, temperature_k, params);

    const double previous_residual = previous_pressure_gpa - target_pressure_gpa;
    const double current_residual = current_pressure_gpa - target_pressure_gpa;
    if (previous_residual == 0.0) {
      return previous_multiplier * params.v0_a3;
    }
    if (previous_residual * current_residual <= 0.0) {
      double lower_volume_a3 = previous_multiplier * params.v0_a3;
      double upper_volume_a3 = current_multiplier * params.v0_a3;
      double lower_residual = previous_residual;
      for (int iteration = 0; iteration < kVolumeSolveMaxIterations; ++iteration) {
        const double midpoint_volume_a3 = 0.5 * (lower_volume_a3 + upper_volume_a3);
        const double midpoint_residual =
            ComputeOxideCellPressureGpa(midpoint_volume_a3, temperature_k, params) -
            target_pressure_gpa;
        if (std::fabs(upper_volume_a3 - lower_volume_a3) <=
            kVolumeSolveRelativeTolerance * params.v0_a3) {
          return midpoint_volume_a3;
        }
        if (lower_residual * midpoint_residual <= 0.0) {
          upper_volume_a3 = midpoint_volume_a3;
        } else {
          lower_volume_a3 = midpoint_volume_a3;
          lower_residual = midpoint_residual;
        }
      }
      return 0.5 * (lower_volume_a3 + upper_volume_a3);
    }

    previous_multiplier = current_multiplier;
    previous_pressure_gpa = current_pressure_gpa;
  }

  throw std::runtime_error("Failed to bracket oxide volume for EOS solve.");
}

[[nodiscard]] double ComputeReactionDeltaVCm3PerMol(double pressure_gpa,
                                                    double temperature_k) {
  const double ferric_cell_volume_a3 =
      SolveOxideCellVolumeA3(pressure_gpa, temperature_k, kFerricCellEos);
  const double ferrous_cell_volume_a3 =
      SolveOxideCellVolumeA3(pressure_gpa, temperature_k, kFerrousCellEos);
  return 0.5 * (ferric_cell_volume_a3 - ferrous_cell_volume_a3) *
         kAngstrom3ToCm3PerMol;
}

[[nodiscard]] double ComputePressureIntegralTermLog10(double pressure_gpa,
                                                      double temperature_k) {
  if (pressure_gpa <= kReferencePressureGpa) {
    return 0.0;
  }
  const double pressure_step_gpa =
      (pressure_gpa - kReferencePressureGpa) / static_cast<double>(kPressureIntegralSteps);
  double integral_gpa_cm3_per_mol = 0.0;
  double previous_delta_v_cm3_per_mol =
      ComputeReactionDeltaVCm3PerMol(kReferencePressureGpa, temperature_k);
  for (int step = 1; step <= kPressureIntegralSteps; ++step) {
    const double current_pressure_gpa =
        kReferencePressureGpa + pressure_step_gpa * static_cast<double>(step);
    const double current_delta_v_cm3_per_mol =
        ComputeReactionDeltaVCm3PerMol(current_pressure_gpa, temperature_k);
    integral_gpa_cm3_per_mol +=
        0.5 * (previous_delta_v_cm3_per_mol + current_delta_v_cm3_per_mol) *
        pressure_step_gpa;
    previous_delta_v_cm3_per_mol = current_delta_v_cm3_per_mol;
  }
  const double integral_j_per_mol = integral_gpa_cm3_per_mol * 1000.0;
  return integral_j_per_mol /
         (kUniversalGasConstantJPerMolK * temperature_k * kLog10Factor);
}

[[nodiscard]] double ComputeLog10FerricRatio(double temperature_k,
                                             double pressure_gpa,
                                             double log10_fugacity_o2) {
  return kFerricFitA * log10_fugacity_o2 + kFerricFitB +
         kFerricFitC / temperature_k -
         ComputeFerricHeatCapacityTermLog10(temperature_k) -
         ComputePressureIntegralTermLog10(pressure_gpa, temperature_k) +
         ComputeFerricCompositionTermLog10(temperature_k);
}

[[nodiscard]] double ComputeFerricFractionAtEquilibrium(double temperature_k,
                                                        double pressure_gpa,
                                                        double delta_iw_eq) {
  const double log10_fugacity_o2 =
      ComputeLog10FugacityFromDeltaIw(delta_iw_eq, temperature_k, pressure_gpa);
  const double log10_r_ferric =
      ComputeLog10FerricRatio(temperature_k, pressure_gpa, log10_fugacity_o2);
  const double r_ferric = std::pow(10.0, log10_r_ferric);
  return r_ferric / (1.0 + r_ferric);
}

[[nodiscard]] double ComputeSurfaceTemperatureK(double equilibration_temperature_k,
                                                double equilibration_pressure_gpa,
                                                const Constants& constants) {
  return std::max(equilibration_temperature_k -
                      kSurfaceAdiabaticGradientKPerGpa * equilibration_pressure_gpa,
                  constants.liquidus_surface_temperature_k);
}

[[nodiscard]] double ComputeDeltaIwSurfaceFromFerricFraction(
    double ferric_fraction, double surface_temperature_k) {
  const double clamped_ferric_fraction =
      std::clamp(ferric_fraction, 1.0e-12, 1.0 - 1.0e-12);
  const double ferric_ratio =
      clamped_ferric_fraction / (1.0 - clamped_ferric_fraction);
  const double log10_fugacity_o2 =
      (std::log10(ferric_ratio) - kFerricFitB -
       kFerricFitC / surface_temperature_k +
       ComputeFerricHeatCapacityTermLog10(surface_temperature_k) -
       ComputeFerricCompositionTermLog10(surface_temperature_k)) /
      kFerricFitA;
  return log10_fugacity_o2 -
         ComputeLog10FugacityIw(surface_temperature_k, kReferencePressureGpa);
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
  if (pressure_gpa < constants.liquidus_slope_break_gpa) {
    return constants.liquidus_surface_temperature_k +
           constants.liquidus_slope_low_pressure_k_gpa * pressure_gpa;
  }

  return constants.liquidus_surface_temperature_k +
         constants.liquidus_slope_low_pressure_k_gpa *
             constants.liquidus_slope_break_gpa +
         constants.liquidus_slope_high_pressure_k_gpa *
             (pressure_gpa - constants.liquidus_slope_break_gpa);
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

void ApplyAtmosphereSpeciationFromDeltaIw(AtmosphereState* atmosphere,
                                        double delta_iw_surface) {
  if (atmosphere == nullptr) {
    throw std::runtime_error(
        "ApplyAtmosphereSpeciationFromDeltaIw received null atmosphere.");
  }

  if (delta_iw_surface > -1.0) {
    atmosphere->carbon_species = CarbonGasSpecies::kCO2;
    atmosphere->hydrogen_species = HydrogenGasSpecies::kH2O;
    atmosphere->nitrogen_species = NitrogenGasSpecies::kN2;
  } else if (delta_iw_surface > -3.0) {
    atmosphere->carbon_species = CarbonGasSpecies::kCO;
    atmosphere->hydrogen_species = HydrogenGasSpecies::kH2;
    atmosphere->nitrogen_species = NitrogenGasSpecies::kN2;
  } else {
    atmosphere->carbon_species = CarbonGasSpecies::kCH4;
    atmosphere->hydrogen_species = HydrogenGasSpecies::kH2;
    atmosphere->nitrogen_species = NitrogenGasSpecies::kN2;
  }
}

[[nodiscard]] ErosionEfficiencies AddErosionEfficiencies(
    const ErosionEfficiencies& left, const ErosionEfficiencies& right) {
  return {
      left.direct_atmosphere_efficiency + right.direct_atmosphere_efficiency,
      left.impactor_vapor_escape_efficiency + right.impactor_vapor_escape_efficiency,
  };
}

[[nodiscard]] ErosionEfficiencies ScaleErosionEfficiencies(
    const ErosionEfficiencies& efficiencies, double factor) {
  return {
      efficiencies.direct_atmosphere_efficiency * factor,
      efficiencies.impactor_vapor_escape_efficiency * factor,
  };
}

template <typename Func>
[[nodiscard]] double SimpsonIntegrateScalar(double x_min, double x_max,
                                            double relative_tolerance,
                                            int max_refinements, Func&& func) {
  if (x_max <= x_min) {
    return 0.0;
  }

  int interval_count = 2;
  auto integrate = [&](int num_intervals) {
    const double dx = (x_max - x_min) / static_cast<double>(num_intervals);
    double weighted_sum = func(x_min) + func(x_max);
    for (int index = 1; index < num_intervals; ++index) {
      const double x = x_min + dx * static_cast<double>(index);
      weighted_sum += (index % 2 == 0 ? 2.0 : 4.0) * func(x);
    }
    return weighted_sum * dx / 3.0;
  };

  double previous_value = integrate(interval_count);
  for (int refinement = 0; refinement < max_refinements; ++refinement) {
    interval_count *= 2;
    const double current_value = integrate(interval_count);
    const double scale = std::max(std::fabs(previous_value), 1.0e-30);
    if (std::fabs(current_value - previous_value) / scale < relative_tolerance) {
      return current_value;
    }
    previous_value = current_value;
  }
  return previous_value;
}

template <typename Func>
[[nodiscard]] ErosionEfficiencies SimpsonIntegrateEfficiencies(
    double x_min, double x_max, double relative_tolerance,
    int max_refinements, Func&& func) {
  if (x_max <= x_min) {
    return {};
  }

  int interval_count = 2;
  auto integrate = [&](int num_intervals) {
    const double dx = (x_max - x_min) / static_cast<double>(num_intervals);
    ErosionEfficiencies weighted_sum =
        AddErosionEfficiencies(func(x_min), func(x_max));
    for (int index = 1; index < num_intervals; ++index) {
      const double x = x_min + dx * static_cast<double>(index);
      const double weight = index % 2 == 0 ? 2.0 : 4.0;
      weighted_sum = AddErosionEfficiencies(
          weighted_sum, ScaleErosionEfficiencies(func(x), weight));
    }
    return ScaleErosionEfficiencies(weighted_sum, dx / 3.0);
  };

  ErosionEfficiencies previous_value = integrate(interval_count);
  for (int refinement = 0; refinement < max_refinements; ++refinement) {
    interval_count *= 2;
    const ErosionEfficiencies current_value = integrate(interval_count);
    const double direct_scale =
        std::max(std::fabs(previous_value.direct_atmosphere_efficiency), 1.0e-30);
    const double vapor_scale =
        std::max(std::fabs(previous_value.impactor_vapor_escape_efficiency), 1.0e-30);
    const double direct_change =
        std::fabs(current_value.direct_atmosphere_efficiency -
                  previous_value.direct_atmosphere_efficiency) /
        direct_scale;
    const double vapor_change =
        std::fabs(current_value.impactor_vapor_escape_efficiency -
                  previous_value.impactor_vapor_escape_efficiency) /
        vapor_scale;
    if (direct_change < relative_tolerance && vapor_change < relative_tolerance) {
      return current_value;
    }
    previous_value = current_value;
  }
  return previous_value;
}

[[nodiscard]] double CarbonSpeciesMolarMassKgPerMol(CarbonGasSpecies species) {
  switch (species) {
    case CarbonGasSpecies::kCO2:
      return 44.0095e-3;
    case CarbonGasSpecies::kCO:
      return 28.0101e-3;
    case CarbonGasSpecies::kCH4:
      return 16.0425e-3;
  }
  throw std::runtime_error("Unknown carbon gas species.");
}

[[nodiscard]] double HydrogenSpeciesMolarMassKgPerMol(HydrogenGasSpecies species) {
  switch (species) {
    case HydrogenGasSpecies::kH2O:
      return 18.01528e-3;
    case HydrogenGasSpecies::kH2:
      return 2.01588e-3;
  }
  throw std::runtime_error("Unknown hydrogen gas species.");
}

[[nodiscard]] double NitrogenSpeciesMolarMassKgPerMol(NitrogenGasSpecies species) {
  switch (species) {
    case NitrogenGasSpecies::kN2:
      return 28.0134e-3;
  }
  throw std::runtime_error("Unknown nitrogen gas species.");
}

[[nodiscard]] double ComputeAtmosphereTotalMolecularMassKg(
    const AtmosphereState& atmosphere) {
  return atmosphere.elements_kg.carbon_kg *
             (CarbonSpeciesMolarMassKgPerMol(atmosphere.carbon_species) /
              12.011e-3) +
         atmosphere.elements_kg.hydrogen_kg *
             (HydrogenSpeciesMolarMassKgPerMol(atmosphere.hydrogen_species) /
              (2.0 * 1.00794e-3)) +
         atmosphere.elements_kg.nitrogen_kg *
             (NitrogenSpeciesMolarMassKgPerMol(atmosphere.nitrogen_species) /
              28.0134e-3);
}

[[nodiscard]] double ComputeAtmosphereTotalMoles(const AtmosphereState& atmosphere) {
  return atmosphere.elements_kg.carbon_kg / 12.011e-3 +
         atmosphere.elements_kg.hydrogen_kg / (2.0 * 1.00794e-3) +
         atmosphere.elements_kg.nitrogen_kg / 28.0134e-3;
}

[[nodiscard]] double ComputeAtmosphereMeanMolecularMolarMassKgPerMol(
    const AtmosphereState& atmosphere) {
  const double total_moles = ComputeAtmosphereTotalMoles(atmosphere);
  if (!(total_moles > 0.0)) {
    return 0.0;
  }
  return ComputeAtmosphereTotalMolecularMassKg(atmosphere) / total_moles;
}

[[nodiscard]] ElementMassKg ComputeImpactorBulkElementMassKg(
    const Impactor& impactor) {
  impactor.Validate();
  if (impactor.bulk_inventory_source != InventorySource::kNone) {
    return PpmToElementMassKg(impactor.bulk_volatiles_ppm, impactor.mass_kg);
  }

  const double impactor_silicate_mass_kg =
      (1.0 - impactor.metallic_fe_fraction) * impactor.mass_kg;
  const double impactor_metal_mass_kg =
      impactor.metallic_fe_fraction * impactor.mass_kg;
  return AddElementMassKg(
      PpmToElementMassKg(impactor.silicate_volatiles_ppm, impactor_silicate_mass_kg),
      PpmToElementMassKg(impactor.metal_volatiles_ppm, impactor_metal_mass_kg));
}

[[nodiscard]] ImpactorVolatileMassesKg ComputeImpactorEquilibratingVolatileMassesKg(
    const Impactor& impactor, double impactor_silicate_mass_kg,
    double impactor_metal_mass_kg) {
  impactor.Validate();
  if (impactor.bulk_inventory_source != InventorySource::kNone) {
    return {
        PpmToElementMassKg(impactor.bulk_volatiles_ppm, impactor.mass_kg),
        {},
    };
  }
  return {
      PpmToElementMassKg(impactor.silicate_volatiles_ppm, impactor_silicate_mass_kg),
      PpmToElementMassKg(impactor.metal_volatiles_ppm, impactor_metal_mass_kg),
  };
}

[[nodiscard]] ElementMassKg ComputeImpactorBulkElementMassFractions(
    const Impactor& impactor) {
  const ElementMassKg impactor_bulk_elements_kg =
      ComputeImpactorBulkElementMassKg(impactor);
  return {
      impactor_bulk_elements_kg.carbon_kg / impactor.mass_kg,
      impactor_bulk_elements_kg.hydrogen_kg / impactor.mass_kg,
      impactor_bulk_elements_kg.nitrogen_kg / impactor.mass_kg,
  };
}

[[nodiscard]] double ComputeMassDistributionWeight(
    double log10_diameter_cm, ErosionSizeDistribution distribution) {
  constexpr double kLog10DiameterMinCm = 3.5;
  constexpr double kLog10DiameterMaxCm = 8.0;
  static const double dminus2_normalization =
      (std::pow(10.0, 2.0 * kLog10DiameterMaxCm) -
       std::pow(10.0, 2.0 * kLog10DiameterMinCm)) /
      (2.0 * kLog10Factor);
  static const double dminus3_normalization =
      (std::pow(10.0, kLog10DiameterMaxCm) -
       std::pow(10.0, kLog10DiameterMinCm)) /
      kLog10Factor;

  switch (distribution) {
    case ErosionSizeDistribution::kDMinus2:
      return std::pow(10.0, 2.0 * log10_diameter_cm) / dminus2_normalization;
    case ErosionSizeDistribution::kDMinus3:
      return std::pow(10.0, log10_diameter_cm) / dminus3_normalization;
  }
  throw std::runtime_error("Unknown erosion size distribution.");
}

[[nodiscard]] double ComputeVelocityDistributionWeight(double velocity_factor) {
  constexpr double kVelocityFactorMin = 1.0;
  constexpr double kVelocityFactorMax = 4.324;
  static const double normalization = SimpsonIntegrateScalar(
      kVelocityFactorMin, kVelocityFactorMax, 1.0e-10, 14,
      [](double factor) {
        const double excess_velocity = factor * factor - 1.0;
        return excess_velocity <= 0.0
                   ? 0.0
                   : std::sqrt(excess_velocity) * std::exp(-0.5 * excess_velocity);
      });

  if (velocity_factor < kVelocityFactorMin || velocity_factor > kVelocityFactorMax) {
    return 0.0;
  }
  const double excess_velocity = velocity_factor * velocity_factor - 1.0;
  return excess_velocity <= 0.0
             ? 0.0
             : std::sqrt(excess_velocity) * std::exp(-0.5 * excess_velocity) /
                   normalization;
}

[[nodiscard]] ErosionEfficiencies ComputeSingleImpactErosionEfficiencies(
    double impactor_diameter_m, double impact_velocity_m_s,
    double escape_velocity_m_s, double atmospheric_scale_height_m,
    double atmospheric_density_kg_m3, const Constants& constants) {
  if (impactor_diameter_m <= 0.0 || impact_velocity_m_s <= 0.0 ||
      escape_velocity_m_s <= 0.0 || atmospheric_scale_height_m <= 0.0 ||
      atmospheric_density_kg_m3 <= 0.0) {
    return {};
  }

  const double diameter_term =
      std::pow(impactor_diameter_m / atmospheric_scale_height_m, 3.0);
  const double density_term =
      constants.erosion_impactor_density_kg_m3 *
      constants.erosion_target_surface_density_kg_m3 /
      (atmospheric_density_kg_m3 *
       (constants.erosion_target_surface_density_kg_m3 +
        constants.erosion_impactor_density_kg_m3));
  const double velocity_term =
      std::pow(impact_velocity_m_s / escape_velocity_m_s, 2.0) - 1.0;
  const double u_term =
      constants.erosion_target_surface_density_kg_m3 /
      constants.erosion_impactor_density_kg_m3 *
      impact_velocity_m_s / escape_velocity_m_s;
  const double impactor_vapor_candidate_2 = 0.07 * u_term;
  const double full_integral = 256.0 / 693.0;
  const double alpha_factor = 5.0 + Pi() - 4.0 / 3.0;

  double direct_efficiency = 0.0;
  if (velocity_term >= 0.0) {
    const double attenuation = atmospheric_density_kg_m3 /
                               constants.erosion_impactor_density_kg_m3 *
                               ((atmospheric_scale_height_m / impactor_diameter_m) +
                                (2.0 * atmospheric_scale_height_m * atmospheric_scale_height_m *
                                 std::sqrt(atmospheric_density_kg_m3 /
                                           constants.erosion_impactor_density_kg_m3)) /
                                    (0.75 * impactor_diameter_m * impactor_diameter_m) +
                                (std::pow(atmospheric_scale_height_m, 3.0) *
                                 atmospheric_density_kg_m3) /
                                    ((std::pow(impactor_diameter_m, 3.0) / 8.0) *
                                     constants.erosion_impactor_density_kg_m3));
    const double surface_impact_velocity_m_s =
        impact_velocity_m_s * std::exp(-attenuation);
    const double plume_velocity_m_s = std::sqrt(26.0) * surface_impact_velocity_m_s;
    if (plume_velocity_m_s >= escape_velocity_m_s) {
      const double velocity_ratio = escape_velocity_m_s / plume_velocity_m_s;
      const double plume_integral =
          -std::pow(velocity_ratio, 11.0) / 11.0 +
          5.0 * std::pow(velocity_ratio, 9.0) / 9.0 -
          10.0 * std::pow(velocity_ratio, 7.0) / 7.0 +
          2.0 * std::pow(velocity_ratio, 5.0) -
          5.0 * std::pow(velocity_ratio, 3.0) / 3.0 + velocity_ratio;
      const double integral_fraction = (full_integral - plume_integral) / full_integral;
      direct_efficiency =
          0.75 * atmospheric_density_kg_m3 /
          constants.erosion_impactor_density_kg_m3 *
          (2.0 * atmospheric_scale_height_m / impactor_diameter_m +
           16.0 * atmospheric_scale_height_m * atmospheric_scale_height_m *
               std::sqrt(atmospheric_density_kg_m3 /
                         constants.erosion_impactor_density_kg_m3) /
               (3.0 * impactor_diameter_m * impactor_diameter_m) +
           16.0 * std::pow(atmospheric_scale_height_m, 3.0) *
               atmospheric_density_kg_m3 /
               (std::pow(impactor_diameter_m, 3.0) *
                constants.erosion_impactor_density_kg_m3)) *
          integral_fraction * alpha_factor;
    }
  }
  direct_efficiency = std::max(0.0, direct_efficiency);

  double impactor_vapor_efficiency = 0.0;
  if (velocity_term > 0.0) {
    const double x_parameter = diameter_term * density_term * velocity_term;
    if (x_parameter > 0.0) {
      const double impactor_vapor_candidate_1 =
          0.035 * u_term * (std::log10(x_parameter) - 1.0);
      impactor_vapor_efficiency =
          std::clamp(std::min(impactor_vapor_candidate_1,
                              impactor_vapor_candidate_2),
                     0.0, 1.0);
    }
  }

  return {direct_efficiency, impactor_vapor_efficiency};
}

[[nodiscard]] ErosionEfficiencies ComputeDiameterAveragedErosionEfficiencies(
    double impact_velocity_m_s, double escape_velocity_m_s,
    double atmospheric_scale_height_m, double atmospheric_density_kg_m3,
    ErosionSizeDistribution distribution, const Constants& constants) {
  constexpr double kLog10DiameterMinCm = 3.5;
  constexpr double kLog10DiameterSplitCm = 3.6874278;
  constexpr double kLog10DiameterMaxCm = 8.0;

  auto integrate_range = [&](double lower, double upper) {
    return SimpsonIntegrateEfficiencies(
        lower, upper, constants.erosion_integration_rel_tolerance,
        constants.erosion_integration_max_refinements,
        [&](double log10_diameter_cm) {
          const double impactor_diameter_m = std::pow(10.0, log10_diameter_cm - 2.0);
          return ScaleErosionEfficiencies(
              ComputeSingleImpactErosionEfficiencies(
                  impactor_diameter_m, impact_velocity_m_s, escape_velocity_m_s,
                  atmospheric_scale_height_m, atmospheric_density_kg_m3,
                  constants),
              ComputeMassDistributionWeight(log10_diameter_cm, distribution));
        });
  };

  return AddErosionEfficiencies(
      integrate_range(kLog10DiameterMinCm, kLog10DiameterSplitCm),
      integrate_range(kLog10DiameterSplitCm, kLog10DiameterMaxCm));
}

[[nodiscard]] ErosionEfficiencies ComputeAverageAtmosphericErosionEfficiencies(
    double escape_velocity_m_s, double atmospheric_scale_height_m,
    double atmospheric_density_kg_m3, ErosionSizeDistribution distribution,
    const Constants& constants) {
  constexpr double kVelocityFactorMin = 1.0;
  constexpr double kVelocityFactorMax = 4.324;
  if (escape_velocity_m_s <= 0.0 || atmospheric_scale_height_m <= 0.0 ||
      atmospheric_density_kg_m3 <= 0.0) {
    return {};
  }

  return SimpsonIntegrateEfficiencies(
      kVelocityFactorMin, kVelocityFactorMax,
      constants.erosion_integration_rel_tolerance,
      constants.erosion_integration_max_refinements,
      [&](double velocity_factor) {
        return ScaleErosionEfficiencies(
            ComputeDiameterAveragedErosionEfficiencies(
                velocity_factor * escape_velocity_m_s, escape_velocity_m_s,
                atmospheric_scale_height_m, atmospheric_density_kg_m3,
                distribution, constants),
            ComputeVelocityDistributionWeight(velocity_factor));
      });
}

[[nodiscard]] AtmosphericErosionOutcome ComputeAccretionPhaseAtmosphericErosion(
    const PlanetState& planet_before, const AtmosphereState& atmosphere_before_erosion,
    const Impactor& impactor, GrowthPhase phase, const Constants& constants) {
  planet_before.Validate();
  impactor.Validate();
  AtmosphericErosionOutcome outcome{};

  const double total_molecular_mass_kg =
      ComputeAtmosphereTotalMolecularMassKg(atmosphere_before_erosion);
  const double mean_molecular_mass_kg_per_mol =
      ComputeAtmosphereMeanMolecularMolarMassKgPerMol(atmosphere_before_erosion);
  if (!(total_molecular_mass_kg > 0.0) ||
      !(mean_molecular_mass_kg_per_mol > 0.0)) {
    return outcome;
  }

  const double planet_radius_m =
      RadiusFromMassKg(planet_before.total_mass_kg(), constants);
  const double surface_gravity_m_s2 =
      Gravity(planet_before.total_mass_kg(), planet_radius_m, constants);
  const double escape_velocity_m_s =
      EscapeVelocity(planet_before.total_mass_kg(), planet_radius_m, constants);
  const double atmosphere_temperature_k =
      phase == GrowthPhase::kLateVeneer
          ? constants.erosion_late_veneer_temperature_k
          : constants.erosion_magma_ocean_temperature_k;
  const double atmospheric_scale_height_m =
      kUniversalGasConstantJPerMolK * atmosphere_temperature_k /
      (mean_molecular_mass_kg_per_mol * surface_gravity_m_s2);
  const double atmospheric_density_kg_m3 =
      total_molecular_mass_kg /
      (4.0 * Pi() * planet_radius_m * planet_radius_m * atmospheric_scale_height_m);

  const ErosionSizeDistribution distribution =
      phase == GrowthPhase::kLateVeneer ? ErosionSizeDistribution::kDMinus3
                                        : ErosionSizeDistribution::kDMinus2;
  const ErosionEfficiencies efficiencies =
      ComputeAverageAtmosphericErosionEfficiencies(
          escape_velocity_m_s, atmospheric_scale_height_m,
          atmospheric_density_kg_m3, distribution, constants);
  outcome.impactor_vapor_escape_efficiency =
      efficiencies.impactor_vapor_escape_efficiency;

  const ElementMassKg impactor_bulk_mass_fractions =
      ComputeImpactorBulkElementMassFractions(impactor);
  auto compute_loss = [&](double atmospheric_element_mass_kg,
                          double impactor_bulk_fraction) {
    const double raw_loss_kg =
        (efficiencies.impactor_vapor_escape_efficiency * impactor_bulk_fraction +
         efficiencies.direct_atmosphere_efficiency * atmospheric_element_mass_kg /
             total_molecular_mass_kg) *
        impactor.mass_kg;
    return std::clamp(raw_loss_kg, 0.0, atmospheric_element_mass_kg);
  };

  outcome.eroded_elements_kg = {
      compute_loss(atmosphere_before_erosion.elements_kg.carbon_kg,
                   impactor_bulk_mass_fractions.carbon_kg),
      compute_loss(atmosphere_before_erosion.elements_kg.hydrogen_kg,
                   impactor_bulk_mass_fractions.hydrogen_kg),
      compute_loss(atmosphere_before_erosion.elements_kg.nitrogen_kg,
                   impactor_bulk_mass_fractions.nitrogen_kg),
  };
  const double eroded_total_molecular_mass_kg =
      outcome.eroded_elements_kg.carbon_kg *
          (CarbonSpeciesMolarMassKgPerMol(atmosphere_before_erosion.carbon_species) /
           12.011e-3) +
      outcome.eroded_elements_kg.hydrogen_kg *
          (HydrogenSpeciesMolarMassKgPerMol(atmosphere_before_erosion.hydrogen_species) /
           (2.0 * 1.00794e-3)) +
      outcome.eroded_elements_kg.nitrogen_kg *
          (NitrogenSpeciesMolarMassKgPerMol(atmosphere_before_erosion.nitrogen_species) /
           28.0134e-3);
  outcome.erosion_fraction = total_molecular_mass_kg > 0.0
                                 ? eroded_total_molecular_mass_kg /
                                       total_molecular_mass_kg
                                 : 0.0;
  return outcome;
}

struct AtmosphereHostGeometry {
  double host_mass_kg = 0.0;
  double radius_m = 0.0;
  double surface_gravity_m_s2 = 0.0;
};

[[nodiscard]] AtmosphereHostGeometry ComputeAtmosphereHostGeometry(
    const PlanetState& planet_before, const Impactor& impactor,
    const Constants& constants) {
  const double host_mass_kg = impactor.phase == GrowthPhase::kGiantImpact
                                  ? planet_before.total_mass_kg() + impactor.mass_kg
                                  : planet_before.total_mass_kg();
  if (!(host_mass_kg > 0.0)) {
    return {};
  }
  const double radius_m = RadiusFromMassKg(host_mass_kg, constants);
  return {
      host_mass_kg,
      radius_m,
      Gravity(host_mass_kg, radius_m, constants),
  };
}

[[nodiscard]] AtmosphericErosionOutcome ComputeGiantImpactAtmosphericErosion(
    const PlanetState& planet_before, const AtmosphereState& atmosphere_before_erosion,
    const Impactor& impactor, const Constants& constants) {
  planet_before.Validate();
  impactor.Validate();
  AtmosphericErosionOutcome outcome{};
  const AtmosphereHostGeometry geometry =
      ComputeAtmosphereHostGeometry(planet_before, impactor, constants);
  if (!(geometry.host_mass_kg > 0.0)) {
    return outcome;
  }

  const double mass_ratio = impactor.mass_kg / geometry.host_mass_kg;
  const double erosion_fraction = std::clamp(
      0.4 * mass_ratio + 1.4 * mass_ratio * mass_ratio -
          0.8 * mass_ratio * mass_ratio * mass_ratio,
      0.0, 1.0);
  outcome.erosion_fraction = erosion_fraction;
  outcome.eroded_elements_kg = {
      atmosphere_before_erosion.elements_kg.carbon_kg * erosion_fraction,
      atmosphere_before_erosion.elements_kg.hydrogen_kg * erosion_fraction,
      atmosphere_before_erosion.elements_kg.nitrogen_kg * erosion_fraction,
  };
  return outcome;
}

[[nodiscard]] Impactor ResolveImpactorForStep(
    const Impactor& impactor,
    const std::optional<DifferentiatedImpactorSnapshot>& precompute_snapshot,
    const RedoxState& target_redox) {
  Impactor resolved = impactor;

  switch (resolved.bulk_inventory_source) {
    case InventorySource::kNone:
      resolved.bulk_volatiles_ppm = {};
      break;
    case InventorySource::kEcBulk:
    case InventorySource::kCiBulk:
      break;
    case InventorySource::kPrecomputeOutputAt0p10EarthMass:
    case InventorySource::kPrecomputeCoreAt0p10EarthMass:
      throw std::runtime_error(
          "Bulk inventory source cannot use precompute-resolved phase snapshots.");
  }

  switch (resolved.silicate_inventory_source) {
    case InventorySource::kNone:
      resolved.silicate_volatiles_ppm = {};
      if (resolved.bulk_inventory_source == InventorySource::kNone) {
        resolved.delta15n_silicate_permil = 0.0;
      }
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
      break;
    case InventorySource::kPrecomputeOutputAt0p10EarthMass:
      throw std::runtime_error(
          "Metal inventory source cannot be precompute silicate snapshot.");
  }

  switch (resolved.alloy_source) {
    case InventorySource::kNone:
      break;
    case InventorySource::kPrecomputeCoreAt0p10EarthMass:
      if (!precompute_snapshot.has_value()) {
        throw std::runtime_error("Impactor requires precompute alloy snapshot.");
      }
      resolved.alloy = precompute_snapshot->alloy;
      break;
    case InventorySource::kEcBulk:
    case InventorySource::kCiBulk:
    case InventorySource::kPrecomputeOutputAt0p10EarthMass:
      throw std::runtime_error(
          "Alloy source cannot be bulk inventory or precompute silicate snapshot.");
  }

  switch (resolved.redox_source) {
    case RedoxSource::kImpactorValue:
      break;
    case RedoxSource::kPrecomputeOutputAt0p10EarthMass:
      if (!precompute_snapshot.has_value()) {
        throw std::runtime_error("Impactor requires precompute redox snapshot.");
      }
      resolved.redox = precompute_snapshot->redox;
      break;
    case RedoxSource::kCurrentTargetMantle:
      resolved.redox = target_redox;
      break;
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
      planet.redox,
  };
}

template <typename Func>
[[nodiscard]] double SolveBisection(double lower_bound, double upper_bound,
                                    double relative_tolerance,
                                    int max_iterations, Func&& residual,
                                    const char* label) {
  if (!(upper_bound >= lower_bound)) {
    throw std::runtime_error(std::string(label) + " requires upper_bound >= lower_bound.");
  }
  if (upper_bound == lower_bound) {
    return lower_bound;
  }

  double lower = lower_bound;
  double upper = upper_bound;
  double f_lower = residual(lower);
  double f_upper = residual(upper);
  if (!std::isfinite(f_lower) || !std::isfinite(f_upper)) {
    throw std::runtime_error(std::string(label) + " residual is not finite at the bracket endpoints.");
  }
  if (f_lower == 0.0) {
    return lower;
  }
  if (f_upper == 0.0) {
    return upper;
  }
  if (f_lower * f_upper > 0.0) {
    throw std::runtime_error(std::string(label) + " root is not bracketed.");
  }

  const double scale = std::max(std::fabs(upper_bound), 1.0);
  for (int iteration = 0; iteration < max_iterations; ++iteration) {
    const double midpoint = 0.5 * (lower + upper);
    const double f_midpoint = residual(midpoint);
    if (!std::isfinite(f_midpoint)) {
      throw std::runtime_error(std::string(label) + " residual became non-finite during bisection.");
    }

    if (f_lower * f_midpoint <= 0.0) {
      upper = midpoint;
      f_upper = f_midpoint;
    } else {
      lower = midpoint;
      f_lower = f_midpoint;
    }

    if (std::fabs(upper - lower) / scale < relative_tolerance) {
      return 0.5 * (lower + upper);
    }
  }

  throw std::runtime_error(std::string(label) + " bisection did not converge.");
}

[[nodiscard]] double ComputeAtmosphericPartialPressurePa(
    double atmospheric_element_mass_kg, double element_molar_mass_kg_per_mol,
    double mean_molar_mass_kg_per_mol, double surface_gravity_m_s2,
    double planet_radius_m) {
  if (atmospheric_element_mass_kg <= 0.0 || element_molar_mass_kg_per_mol <= 0.0 ||
      mean_molar_mass_kg_per_mol <= 0.0 || surface_gravity_m_s2 <= 0.0 ||
      planet_radius_m <= 0.0) {
    return 0.0;
  }

  const double molecule_moles =
      atmospheric_element_mass_kg / element_molar_mass_kg_per_mol;
  return surface_gravity_m_s2 * mean_molar_mass_kg_per_mol * molecule_moles /
         (4.0 * Pi() * planet_radius_m * planet_radius_m);
}

[[nodiscard]] double ComputeHenrySilicateConcentrationPpm(double partial_pressure_pa,
                                                          double solubility_ppm,
                                                          double stoichiometric_x) {
  if (partial_pressure_pa <= 0.0) {
    return 0.0;
  }
  if (solubility_ppm < 0.0 || stoichiometric_x <= 0.0) {
    throw std::runtime_error("Henry solubility parameters are invalid.");
  }
  return solubility_ppm *
         std::pow(partial_pressure_pa / 1.0e6, 1.0 / stoichiometric_x);
}

[[nodiscard]] double SelectCarbonSolubilityPpm(CarbonGasSpecies species,
                                               const Constants& constants) {
  switch (species) {
    case CarbonGasSpecies::kCO2:
      return constants.carbon_solubility_oxidized_ppm;
    case CarbonGasSpecies::kCO:
      return constants.carbon_solubility_intermediate_ppm;
    case CarbonGasSpecies::kCH4:
      return constants.carbon_solubility_reduced_ppm;
  }
  throw std::runtime_error("Unknown carbon gas species.");
}

[[nodiscard]] double SelectNitrogenSolubilityPpm(CarbonGasSpecies carbon_species,
                                                 const Constants& constants) {
  switch (carbon_species) {
    case CarbonGasSpecies::kCO2:
      return constants.nitrogen_solubility_oxidized_ppm;
    case CarbonGasSpecies::kCO:
      return constants.nitrogen_solubility_intermediate_ppm;
    case CarbonGasSpecies::kCH4:
      return constants.nitrogen_solubility_reduced_ppm;
  }
  throw std::runtime_error("Unknown carbon gas species for nitrogen solubility.");
}

[[nodiscard]] double SelectHydrogenHenrySolubilityPpmSqrtMpa(
    HydrogenGasSpecies species, const Constants& constants) {
  switch (species) {
    case HydrogenGasSpecies::kH2:
      return constants.hydrogen_h2_solubility_ppm_sqrt_mpa;
    case HydrogenGasSpecies::kH2O:
      throw std::runtime_error(
          "H2O solubility must be solved with the Moore (1998) model.");
  }
  throw std::runtime_error("Unknown hydrogen gas species.");
}

[[nodiscard]] ElementPartitionKg SolveHenryAtmosphereSilicateMetalPartition(
    double total_element_kg, double silicate_mass_kg, double metal_mass_kg,
    double partition_coefficient, double solubility_ppm, double stoichiometric_x,
    double element_molar_mass_kg_per_mol, double mean_molar_mass_kg_per_mol,
    double surface_gravity_m_s2, double planet_radius_m,
    const Constants& constants, const char* label) {
  if (total_element_kg < 0.0 || silicate_mass_kg < 0.0 || metal_mass_kg < 0.0 ||
      partition_coefficient < 0.0) {
    throw std::runtime_error(std::string(label) +
                             " received a negative mass or partition coefficient.");
  }
  if (total_element_kg == 0.0) {
    return {};
  }
  if (silicate_mass_kg <= 0.0) {
    return {total_element_kg, 0.0, 0.0};
  }

  auto residual = [&](double atmosphere_element_mass_kg) {
    const double partial_pressure_pa = ComputeAtmosphericPartialPressurePa(
        atmosphere_element_mass_kg, element_molar_mass_kg_per_mol,
        mean_molar_mass_kg_per_mol, surface_gravity_m_s2, planet_radius_m);
    const double silicate_concentration_ppm = ComputeHenrySilicateConcentrationPpm(
        partial_pressure_pa, solubility_ppm, stoichiometric_x);
    const double silicate_mass_fraction =
        PpmToMassFraction(silicate_concentration_ppm);
    if (silicate_mass_fraction < 0.0 || silicate_mass_fraction > 1.0) {
      throw std::runtime_error(std::string(label) +
                               " silicate mass fraction left [0, 1].");
    }
    return total_element_kg - atmosphere_element_mass_kg -
           silicate_mass_fraction * silicate_mass_kg -
           partition_coefficient * silicate_mass_fraction * metal_mass_kg;
  };

  const double atmosphere_element_mass_kg = SolveBisection(
      0.0, total_element_kg, constants.step2_root_rel_tolerance,
      constants.step2_root_max_iterations, residual, label);
  const double partial_pressure_pa = ComputeAtmosphericPartialPressurePa(
      atmosphere_element_mass_kg, element_molar_mass_kg_per_mol,
      mean_molar_mass_kg_per_mol, surface_gravity_m_s2, planet_radius_m);
  const double silicate_concentration_ppm = ComputeHenrySilicateConcentrationPpm(
      partial_pressure_pa, solubility_ppm, stoichiometric_x);
  const double silicate_mass_fraction = PpmToMassFraction(silicate_concentration_ppm);
  const double silicate_element_mass_kg = silicate_mass_fraction * silicate_mass_kg;
  const double metal_element_mass_kg =
      partition_coefficient * silicate_mass_fraction * metal_mass_kg;
  const double atmosphere_mass_kg =
      total_element_kg - silicate_element_mass_kg - metal_element_mass_kg;
  if (atmosphere_mass_kg < -1.0e-9 * std::max(total_element_kg, 1.0)) {
    throw std::runtime_error(std::string(label) +
                             " atmosphere mass became negative after Henry solve.");
  }

  return {
      std::max(0.0, atmosphere_mass_kg),
      silicate_element_mass_kg,
      metal_element_mass_kg,
  };
}

[[nodiscard]] ElementPartitionKg SolveMooreHydrogenAtmosphereSilicateMetalPartition(
    double total_hydrogen_kg, double silicate_mass_kg, double metal_mass_kg,
    double partition_coefficient, double temperature_k,
    double mean_molar_mass_kg_per_mol, double surface_gravity_m_s2,
    double planet_radius_m, double carbon_partial_pressure_pa,
    double nitrogen_partial_pressure_pa, const Constants& constants) {
  if (total_hydrogen_kg < 0.0 || silicate_mass_kg < 0.0 || metal_mass_kg < 0.0 ||
      partition_coefficient < 0.0) {
    throw std::runtime_error(
        "SolveMooreHydrogenAtmosphereSilicateMetalPartition received invalid inputs.");
  }
  if (total_hydrogen_kg == 0.0) {
    return {};
  }
  if (silicate_mass_kg <= 0.0) {
    return {total_hydrogen_kg, 0.0, 0.0};
  }

  const double metal_to_silicate_ratio = metal_mass_kg / silicate_mass_kg;
  const double upper_silicate_hydrogen_kg =
      total_hydrogen_kg /
      (1.0 + partition_coefficient * metal_to_silicate_ratio);
  const double epsilon = std::max(1.0e-30, total_hydrogen_kg * 1.0e-12);
  double lower_bound = std::min(epsilon, 0.5 * upper_silicate_hydrogen_kg);
  double upper_bound = upper_silicate_hydrogen_kg - lower_bound;
  if (!(upper_bound > lower_bound)) {
    lower_bound = 0.25 * upper_silicate_hydrogen_kg;
    upper_bound = 0.75 * upper_silicate_hydrogen_kg;
  }
  if (!(upper_bound > lower_bound)) {
    throw std::runtime_error("Moore H solver could not define a valid bracket.");
  }

  const double moore_b_sum =
      constants.moore_b_al * constants.moore_x_al +
      constants.moore_b_fe * constants.moore_x_fe +
      constants.moore_b_na * constants.moore_x_na;

  auto residual = [&](double silicate_hydrogen_kg) {
    const double metal_hydrogen_kg =
        partition_coefficient * metal_to_silicate_ratio * silicate_hydrogen_kg;
    const double atmosphere_hydrogen_kg =
        total_hydrogen_kg - silicate_hydrogen_kg - metal_hydrogen_kg;
    const double water_partial_pressure_pa = ComputeAtmosphericPartialPressurePa(
        atmosphere_hydrogen_kg, kMolarMassHydrogenElementKgPerMol,
        mean_molar_mass_kg_per_mol, surface_gravity_m_s2, planet_radius_m);
    if (!(water_partial_pressure_pa > 0.0) || !(silicate_hydrogen_kg > 0.0)) {
      return std::numeric_limits<double>::quiet_NaN();
    }

    return constants.moore_c * std::log(water_partial_pressure_pa / 1.0e6) +
           moore_b_sum *
               (water_partial_pressure_pa + carbon_partial_pressure_pa +
                nitrogen_partial_pressure_pa) /
               (temperature_k * 1.0e6) -
           2.0 * std::log((kMolarMassH2OKgPerMol * silicate_hydrogen_kg) /
                          (constants.moore_r_c_to_x * silicate_mass_kg *
                           kMolarMassHydrogenElementKgPerMol)) +
           constants.moore_a / temperature_k + constants.moore_d;
  };

  const double silicate_hydrogen_kg = SolveBisection(
      lower_bound, upper_bound, constants.step2_root_rel_tolerance,
      constants.step2_root_max_iterations, residual, "Moore hydrogen partition");
  const double metal_hydrogen_kg =
      partition_coefficient * metal_to_silicate_ratio * silicate_hydrogen_kg;
  const double atmosphere_hydrogen_kg =
      total_hydrogen_kg - silicate_hydrogen_kg - metal_hydrogen_kg;
  if (atmosphere_hydrogen_kg < -1.0e-9 * std::max(total_hydrogen_kg, 1.0)) {
    throw std::runtime_error(
        "Atmosphere hydrogen mass became negative after Moore solve.");
  }

  return {
      std::max(0.0, atmosphere_hydrogen_kg),
      silicate_hydrogen_kg,
      metal_hydrogen_kg,
  };
}

[[nodiscard]] Step2EquilibriumState SolveThreePhaseVolatileEquilibrium(
    const ElementMassKg& total_equilibrating_volatiles_kg,
    const AtmosphereState& initial_mean_guess_state, double silicate_mass_kg,
    double metal_mass_kg, double carbon_partition_coefficient,
    double hydrogen_partition_coefficient, double nitrogen_partition_coefficient,
    double temperature_k, double surface_gravity_m_s2, double planet_radius_m,
    const Constants& constants) {
  Step2EquilibriumState state{};
  if (silicate_mass_kg < 0.0 || metal_mass_kg < 0.0) {
    throw std::runtime_error("SolveThreePhaseVolatileEquilibrium requires non-negative host masses.");
  }

  AtmosphereState mean_guess_state = initial_mean_guess_state;
  double mean_molar_mass_kg_per_mol =
      ComputeAtmosphereMeanMolecularMolarMassKgPerMol(mean_guess_state);
  if (!(mean_molar_mass_kg_per_mol > 0.0)) {
    mean_guess_state.elements_kg = total_equilibrating_volatiles_kg;
    mean_molar_mass_kg_per_mol =
        ComputeAtmosphereMeanMolecularMolarMassKgPerMol(mean_guess_state);
  }
  if (!(mean_molar_mass_kg_per_mol > 0.0)) {
    state.mean_molar_mass_kg_per_mol = 0.0;
    return state;
  }

  const double carbon_solubility_ppm =
      SelectCarbonSolubilityPpm(mean_guess_state.carbon_species, constants);
  const double nitrogen_solubility_ppm =
      SelectNitrogenSolubilityPpm(mean_guess_state.carbon_species, constants);

  for (int iteration = 0; iteration < constants.step2_outer_max_iterations;
       ++iteration) {
    state.carbon = SolveHenryAtmosphereSilicateMetalPartition(
        total_equilibrating_volatiles_kg.carbon_kg, silicate_mass_kg,
        metal_mass_kg, carbon_partition_coefficient, carbon_solubility_ppm, 1.0,
        kMolarMassCarbonElementKgPerMol, mean_molar_mass_kg_per_mol,
        surface_gravity_m_s2, planet_radius_m, constants, "carbon partition");
    state.nitrogen = SolveHenryAtmosphereSilicateMetalPartition(
        total_equilibrating_volatiles_kg.nitrogen_kg, silicate_mass_kg,
        metal_mass_kg, nitrogen_partition_coefficient, nitrogen_solubility_ppm,
        1.0, kMolarMassNitrogenElementKgPerMol,
        mean_molar_mass_kg_per_mol, surface_gravity_m_s2, planet_radius_m,
        constants, "nitrogen partition");

    state.carbon_partial_pressure_pa = ComputeAtmosphericPartialPressurePa(
        state.carbon.atmosphere_kg, kMolarMassCarbonElementKgPerMol,
        mean_molar_mass_kg_per_mol, surface_gravity_m_s2, planet_radius_m);
    state.nitrogen_partial_pressure_pa = ComputeAtmosphericPartialPressurePa(
        state.nitrogen.atmosphere_kg, kMolarMassNitrogenElementKgPerMol,
        mean_molar_mass_kg_per_mol, surface_gravity_m_s2, planet_radius_m);

    switch (mean_guess_state.hydrogen_species) {
      case HydrogenGasSpecies::kH2O:
        state.hydrogen = SolveMooreHydrogenAtmosphereSilicateMetalPartition(
            total_equilibrating_volatiles_kg.hydrogen_kg, silicate_mass_kg,
            metal_mass_kg, hydrogen_partition_coefficient, temperature_k,
            mean_molar_mass_kg_per_mol, surface_gravity_m_s2, planet_radius_m,
            state.carbon_partial_pressure_pa, state.nitrogen_partial_pressure_pa,
            constants);
        break;
      case HydrogenGasSpecies::kH2:
        state.hydrogen = SolveHenryAtmosphereSilicateMetalPartition(
            total_equilibrating_volatiles_kg.hydrogen_kg, silicate_mass_kg,
            metal_mass_kg, hydrogen_partition_coefficient,
            SelectHydrogenHenrySolubilityPpmSqrtMpa(
                mean_guess_state.hydrogen_species, constants),
            2.0, kMolarMassHydrogenElementKgPerMol,
            mean_molar_mass_kg_per_mol, surface_gravity_m_s2, planet_radius_m,
            constants, "hydrogen partition");
        break;
    }

    state.hydrogen_partial_pressure_pa = ComputeAtmosphericPartialPressurePa(
        state.hydrogen.atmosphere_kg, kMolarMassHydrogenElementKgPerMol,
        mean_molar_mass_kg_per_mol, surface_gravity_m_s2, planet_radius_m);

    AtmosphereState updated_atmosphere = mean_guess_state;
    updated_atmosphere.elements_kg = {
        state.carbon.atmosphere_kg,
        state.hydrogen.atmosphere_kg,
        state.nitrogen.atmosphere_kg,
    };
    const double new_mean_molar_mass_kg_per_mol =
        ComputeAtmosphereMeanMolecularMolarMassKgPerMol(updated_atmosphere);
    state.outer_iteration_count = iteration + 1;
    state.mean_molar_mass_kg_per_mol = new_mean_molar_mass_kg_per_mol;

    if (!(new_mean_molar_mass_kg_per_mol > 0.0)) {
      return state;
    }
    if (std::fabs(new_mean_molar_mass_kg_per_mol - mean_molar_mass_kg_per_mol) /
            mean_molar_mass_kg_per_mol <
        constants.step2_outer_rel_tolerance) {
      return state;
    }

    mean_molar_mass_kg_per_mol = new_mean_molar_mass_kg_per_mol;
  }

  throw std::runtime_error(
      "SolveThreePhaseVolatileEquilibrium outer iteration did not converge.");
}

[[nodiscard]] StepResult RunOneImpactStep(std::size_t index, const PlanetState& planet_before,
                                          const Impactor& resolved_impactor,
                                          const Constants& constants) {
  planet_before.Validate();
  resolved_impactor.Validate();

  StepResult result{};
  result.index = index;
  result.resolved_impactor = resolved_impactor;
  result.melt_pool =
      ComputeMeltPoolState(planet_before, resolved_impactor, constants);

  const double target_unmelted_mantle_mass_kg =
      planet_before.mantle.mass_kg - result.melt_pool.target_melt_mass_kg;
  if (target_unmelted_mantle_mass_kg < -1.0e-6) {
    throw std::runtime_error("Melt pool consumed more mantle mass than available.");
  }

  const SilicateIronInventory target_melt_iron = ComputeSilicateIronInventory(
      result.melt_pool.target_melt_mass_kg, planet_before.redox);
  const SilicateIronInventory impactor_silicate_iron =
      ComputeSilicateIronInventory(result.melt_pool.impactor_silicate_mass_kg,
                                   resolved_impactor.redox);
  const SilicateIronInventory mixed_preoxidation_iron =
      AddSilicateIronInventory(target_melt_iron, impactor_silicate_iron);
  const double ferric_fraction_old = mixed_preoxidation_iron.ferric_fraction();
  const double ferric_fraction_equilibrium = std::max(
      ferric_fraction_old,
      ComputeFerricFractionAtEquilibrium(result.melt_pool.temperature_k,
                                         result.melt_pool.pressure_gpa,
                                         planet_before.redox.delta_iw_eq));
  const double max_reaction_extent_mol = mixed_preoxidation_iron.feo_mol / 3.0;
  const double raw_reaction_extent_mol =
      mixed_preoxidation_iron.total_fe_mol() > 0.0
          ? mixed_preoxidation_iron.total_fe_mol() *
                (ferric_fraction_equilibrium - ferric_fraction_old) /
                (2.0 + ferric_fraction_equilibrium)
          : 0.0;
  const double reaction_extent_mol =
      std::clamp(raw_reaction_extent_mol, 0.0, std::max(0.0, max_reaction_extent_mol));

  SilicateIronInventory equilibrated_melt_iron = mixed_preoxidation_iron;
  equilibrated_melt_iron.feo_mol -= 3.0 * reaction_extent_mol;
  equilibrated_melt_iron.feo1p5_mol += 2.0 * reaction_extent_mol;
  const double oxidized_fe0_mass_kg = reaction_extent_mol * kMolarMassFeKgPerMol;
  const double surface_temperature_k = ComputeSurfaceTemperatureK(
      result.melt_pool.temperature_k, result.melt_pool.pressure_gpa, constants);
  const double delta_iw_surface = ComputeDeltaIwSurfaceFromFerricFraction(
      equilibrated_melt_iron.ferric_fraction(), surface_temperature_k);
  const RedoxState equilibrated_melt_redox = MakeRedoxStateFromIronInventory(
      equilibrated_melt_iron, planet_before.redox.delta_iw_eq, delta_iw_surface,
      oxidized_fe0_mass_kg);

  result.dn_metal_silicate = ComputeNitrogenPartitionCoefficient(
      result.melt_pool.temperature_k, result.melt_pool.pressure_gpa,
      planet_before.redox.delta_iw_eq, resolved_impactor.alloy, constants);
  result.dh_metal_silicate = ComputeHydrogenPartitionCoefficient(
      result.melt_pool.temperature_k, result.melt_pool.pressure_gpa,
      planet_before.redox.delta_iw_eq, resolved_impactor.alloy, constants);
  result.xo_metal = ComputeOxygenMoleFractionInMetal(
      result.melt_pool.temperature_k, result.melt_pool.pressure_gpa,
      equilibrated_melt_redox.x_feo_silicate, resolved_impactor.alloy, constants);
  result.dc_metal_silicate = ComputeCarbonPartitionCoefficient(
      result.melt_pool.temperature_k, result.melt_pool.pressure_gpa,
      planet_before.redox.delta_iw_eq, resolved_impactor.alloy.x_s, result.xo_metal,
      constants);

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
  const ImpactorVolatileMassesKg impactor_equilibrating_volatiles_kg =
      ComputeImpactorEquilibratingVolatileMassesKg(
          resolved_impactor, result.melt_pool.impactor_silicate_mass_kg,
          result.melt_pool.impactor_metal_mass_kg);
  const ElementMassKg impactor_silicate_volatiles_kg =
      impactor_equilibrating_volatiles_kg.silicate_kg;
  const ElementMassKg impactor_metal_volatiles_kg =
      impactor_equilibrating_volatiles_kg.metal_kg;

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

  AtmosphereState equilibrium_species_state = planet_before.atmosphere;
  ApplyAtmosphereSpeciationFromDeltaIw(&equilibrium_species_state, delta_iw_surface);
  equilibrium_species_state.elements_kg = old_atmosphere_volatiles_kg;

  const ElementMassKg total_equilibrating_volatiles_kg{
      old_atmosphere_volatiles_kg.carbon_kg + mixed_silicate_volatiles_kg.carbon_kg +
          impactor_metal_volatiles_kg.carbon_kg,
      old_atmosphere_volatiles_kg.hydrogen_kg +
          mixed_silicate_volatiles_kg.hydrogen_kg +
          impactor_metal_volatiles_kg.hydrogen_kg,
      old_atmosphere_volatiles_kg.nitrogen_kg +
          mixed_silicate_volatiles_kg.nitrogen_kg +
          impactor_metal_volatiles_kg.nitrogen_kg,
  };
  const AtmosphereHostGeometry partition_geometry =
      ComputeAtmosphereHostGeometry(planet_before, resolved_impactor, constants);
  const Step2EquilibriumState step2_equilibrium =
      SolveThreePhaseVolatileEquilibrium(
          total_equilibrating_volatiles_kg, equilibrium_species_state,
          result.melt_pool.equilibrating_silicate_mass_kg,
          result.melt_pool.impactor_metal_mass_kg, result.dc_metal_silicate,
          result.dh_metal_silicate, result.dn_metal_silicate,
          result.melt_pool.temperature_k, partition_geometry.surface_gravity_m_s2,
          partition_geometry.radius_m, constants);

  const ElementPartitionKg carbon_partition = step2_equilibrium.carbon;
  const ElementPartitionKg hydrogen_partition = step2_equilibrium.hydrogen;
  const ElementPartitionKg nitrogen_partition = step2_equilibrium.nitrogen;

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
    const double delta15n_metal_silicate =
        ComputeNitrogenIsotopeFractionationPermil(
            result.melt_pool.temperature_k, planet_before.redox.delta_iw_eq,
            constants);
    result.delta15n_silicate_permil =
        delta15n_bulk_permil -
        (metal_after_equilibrium_kg.nitrogen_kg / total_n_before_equilibrium_kg) *
            delta15n_metal_silicate;
    result.delta15n_metal_permil =
        result.delta15n_silicate_permil + delta15n_metal_silicate;
  }

  PlanetState planet_after = planet_before;
  planet_after.core.mass_kg = planet_before.core.mass_kg +
                              result.melt_pool.impactor_metal_mass_kg +
                              oxidized_fe0_mass_kg;
  planet_after.mantle.mass_kg = planet_before.mantle.mass_kg +
                                result.melt_pool.impactor_silicate_mass_kg -
                                oxidized_fe0_mass_kg;

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
  const SilicateIronInventory unmelted_mantle_iron = ComputeSilicateIronInventory(
      std::max(0.0, target_unmelted_mantle_mass_kg), planet_before.redox);
  planet_after.redox = MakeRedoxStateFromIronInventory(
      AddSilicateIronInventory(unmelted_mantle_iron, equilibrated_melt_iron),
      planet_before.redox.delta_iw_eq, delta_iw_surface, oxidized_fe0_mass_kg);
  ApplyAtmosphereSpeciationFromDeltaIw(&planet_after.atmosphere,
                                       planet_after.redox.delta_iw_surface);

  AtmosphericErosionOutcome erosion_outcome{};
  switch (resolved_impactor.phase) {
    case GrowthPhase::kPrecomputeAccretion:
      erosion_outcome = ComputeAccretionPhaseAtmosphericErosion(
          planet_before, planet_after.atmosphere, resolved_impactor,
          GrowthPhase::kPrecomputeAccretion, constants);
      break;
    case GrowthPhase::kGiantImpact:
      erosion_outcome = ComputeGiantImpactAtmosphericErosion(
          planet_before, planet_after.atmosphere, resolved_impactor, constants);
      break;
    case GrowthPhase::kLateVeneer:
      erosion_outcome = ComputeAccretionPhaseAtmosphericErosion(
          planet_before, planet_after.atmosphere, resolved_impactor,
          GrowthPhase::kLateVeneer, constants);
      break;
  }
  result.erosion_fraction = erosion_outcome.erosion_fraction;
  result.eroded_volatiles_kg = erosion_outcome.eroded_elements_kg;
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
    Impactor resolved_impactor = ResolveImpactorForStep(
        scenario.at(index), precompute_snapshot, current_planet.redox);
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
    run.history.push_back(
        CaptureHistoryRecord(index, current_planet, step_result, constants));
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
      precompute_snapshot = CaptureDifferentiatedImpactorSnapshot(
          current_planet, current_precompute_alloy);
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

[[nodiscard]] std::string ToString(RedoxSource redox_source) {
  switch (redox_source) {
    case RedoxSource::kImpactorValue:
      return "impactor_value";
    case RedoxSource::kPrecomputeOutputAt0p10EarthMass:
      return "precompute_output_at_0p10";
    case RedoxSource::kCurrentTargetMantle:
      return "current_target_mantle";
  }
  throw std::runtime_error("Unknown RedoxSource.");
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
  impactor.bulk_inventory_source = InventorySource::kEcBulk;
  impactor.silicate_inventory_source = InventorySource::kNone;
  impactor.metal_inventory_source = InventorySource::kNone;
  impactor.alloy_source = InventorySource::kNone;
  impactor.bulk_volatiles_ppm = constants.ec_bulk_volatiles_ppm;
  impactor.silicate_volatiles_ppm = {};
  impactor.metal_volatiles_ppm = {};
  impactor.delta15n_silicate_permil = constants.ec_delta15n_permil;
  impactor.delta15n_metal_permil = std::numeric_limits<double>::quiet_NaN();
  impactor.alloy = constants.representative_ec_alloy;
  impactor.redox_source = RedoxSource::kImpactorValue;
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
  impactor.bulk_inventory_source = InventorySource::kNone;
  impactor.silicate_inventory_source =
      InventorySource::kPrecomputeOutputAt0p10EarthMass;
  impactor.metal_inventory_source = InventorySource::kPrecomputeCoreAt0p10EarthMass;
  impactor.alloy_source = InventorySource::kPrecomputeCoreAt0p10EarthMass;
  impactor.bulk_volatiles_ppm = {};
  impactor.silicate_volatiles_ppm = {};
  impactor.metal_volatiles_ppm = {};
  impactor.delta15n_silicate_permil = constants.ec_delta15n_permil;
  impactor.delta15n_metal_permil = std::numeric_limits<double>::quiet_NaN();
  impactor.alloy = constants.representative_ec_alloy;
  impactor.redox_source = RedoxSource::kPrecomputeOutputAt0p10EarthMass;
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
  impactor.bulk_inventory_source = InventorySource::kCiBulk;
  impactor.silicate_inventory_source = InventorySource::kNone;
  impactor.metal_inventory_source = InventorySource::kNone;
  impactor.alloy_source = InventorySource::kPrecomputeCoreAt0p10EarthMass;
  impactor.bulk_volatiles_ppm = constants.ci_bulk_volatiles_ppm;
  impactor.silicate_volatiles_ppm = {};
  impactor.metal_volatiles_ppm = {};
  impactor.delta15n_silicate_permil = constants.ci_delta15n_permil;
  impactor.delta15n_metal_permil = std::numeric_limits<double>::quiet_NaN();
  impactor.alloy = constants.representative_ec_alloy;
  impactor.redox_source = RedoxSource::kCurrentTargetMantle;
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
  impactor.bulk_inventory_source = InventorySource::kEcBulk;
  impactor.silicate_inventory_source = InventorySource::kNone;
  impactor.metal_inventory_source = InventorySource::kNone;
  impactor.alloy_source = InventorySource::kPrecomputeCoreAt0p10EarthMass;
  impactor.bulk_volatiles_ppm = constants.ec_bulk_volatiles_ppm;
  impactor.silicate_volatiles_ppm = {};
  impactor.metal_volatiles_ppm = {};
  impactor.delta15n_silicate_permil = constants.ec_delta15n_permil;
  impactor.delta15n_metal_permil = std::numeric_limits<double>::quiet_NaN();
  impactor.alloy = constants.representative_ec_alloy;
  impactor.redox_source = RedoxSource::kImpactorValue;
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
         << " bulk_source=" << ToString(impactor.bulk_inventory_source)
         << " silicate_source=" << ToString(impactor.silicate_inventory_source)
         << " metal_source=" << ToString(impactor.metal_inventory_source)
         << " alloy_source=" << ToString(impactor.alloy_source)
         << " redox_source=" << ToString(impactor.redox_source);
  return stream.str();
}

struct CommandLineOptions {
  bool quiet = false;
  std::string history_csv_path;
};

void WriteHistoryCsv(const std::string& path, const SimplifiedRunResult& run) {
  std::ofstream output(path);
  if (!output) {
    throw std::runtime_error("Failed to open history CSV for writing: " + path);
  }

  output << std::setprecision(15);
  output << "step_index,label,phase,impactor_class,planet_mass_before_earth,"
         << "planet_mass_after_earth,mantle_mass_after_earth,core_mass_after_earth,"
         << "bse_carbon_ppm,bse_hydrogen_ppm,bse_nitrogen_ppm,"
         << "mantle_carbon_ppm,mantle_hydrogen_ppm,mantle_nitrogen_ppm,"
         << "core_carbon_ppm,core_hydrogen_ppm,core_nitrogen_ppm,"
         << "atmosphere_carbon_ppm_silicate_earth,"
         << "atmosphere_hydrogen_ppm_silicate_earth,"
         << "atmosphere_nitrogen_ppm_silicate_earth,"
         << "mantle_delta15n_permil,core_delta15n_permil,"
         << "atmosphere_delta15n_permil\n";

  for (const HistoryRecord& row : run.history) {
    output << row.step_index << ','
           << row.label << ','
           << ToString(row.phase) << ','
           << ToString(row.impactor_class) << ','
           << row.planet_mass_before_earth << ','
           << row.planet_mass_after_earth << ','
           << row.mantle_mass_after_earth << ','
           << row.core_mass_after_earth << ','
           << row.bse_carbon_ppm << ','
           << row.bse_hydrogen_ppm << ','
           << row.bse_nitrogen_ppm << ','
           << row.mantle_carbon_ppm << ','
           << row.mantle_hydrogen_ppm << ','
           << row.mantle_nitrogen_ppm << ','
           << row.core_carbon_ppm << ','
           << row.core_hydrogen_ppm << ','
           << row.core_nitrogen_ppm << ','
           << row.atmosphere_carbon_ppm_silicate_earth << ','
           << row.atmosphere_hydrogen_ppm_silicate_earth << ','
           << row.atmosphere_nitrogen_ppm_silicate_earth << ','
           << row.mantle_delta15n_permil << ','
           << row.core_delta15n_permil << ','
           << row.atmosphere_delta15n_permil << '\n';
  }
}

void PrintUsage(const char* argv0) {
  std::cout << "Usage: " << argv0
            << " [--quiet] [--history-csv PATH]\n";
}

[[nodiscard]] CommandLineOptions ParseCommandLine(int argc, char** argv) {
  CommandLineOptions options{};
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--quiet") {
      options.quiet = true;
      continue;
    }
    if (arg == "--help" || arg == "-h") {
      PrintUsage(argv[0]);
      std::exit(0);
    }
    if (arg == "--history-csv") {
      if (i + 1 >= argc) {
        throw std::runtime_error("--history-csv requires a path argument.");
      }
      options.history_csv_path = argv[++i];
      continue;
    }
    const std::string prefix = "--history-csv=";
    if (arg.rfind(prefix, 0) == 0) {
      options.history_csv_path = arg.substr(prefix.size());
      continue;
    }
    throw std::runtime_error("Unknown command-line option: " + arg);
  }
  return options;
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
  std::cout << "simple_run_completed_steps=" << run.completed_steps << std::endl;
  std::cout << "simple_run_has_precompute_snapshot="
            << (run.has_precompute_snapshot ? "yes" : "no") << std::endl;
  std::cout << "simple_run_final_planet_mass_earth="
            << KgToEarthMass(run.final_planet.total_mass_kg(), constants) << std::endl;
  std::cout << "simple_run_final_mantle_mass_earth="
            << KgToEarthMass(run.final_planet.mantle.mass_kg, constants) << std::endl;
  std::cout << "simple_run_final_core_mass_earth="
            << KgToEarthMass(run.final_planet.core.mass_kg, constants) << std::endl;
  std::cout << "simple_run_final_mantle_C_ppm="
            << run.final_planet.mantle.volatiles_ppm.carbon_ppm << std::endl;
  std::cout << "simple_run_final_mantle_H_ppm="
            << run.final_planet.mantle.volatiles_ppm.hydrogen_ppm << std::endl;
  std::cout << "simple_run_final_mantle_N_ppm="
            << run.final_planet.mantle.volatiles_ppm.nitrogen_ppm << std::endl;
  std::cout << "simple_run_final_core_C_ppm="
            << run.final_planet.core.volatiles_ppm.carbon_ppm << std::endl;
  std::cout << "simple_run_final_core_H_ppm="
            << run.final_planet.core.volatiles_ppm.hydrogen_ppm << std::endl;
  std::cout << "simple_run_final_core_N_ppm="
            << run.final_planet.core.volatiles_ppm.nitrogen_ppm << std::endl;
  std::cout << "simple_run_final_mantle_fe3_over_fe_total="
            << run.final_planet.redox.fe3_over_fe_total << std::endl;
  std::cout << "simple_run_final_mantle_x_feo_silicate="
            << run.final_planet.redox.x_feo_silicate << std::endl;
  std::cout << "simple_run_final_delta_iw_surface="
            << run.final_planet.redox.delta_iw_surface << std::endl;
  std::cout << "simple_run_final_oxidized_fe0_mass_kg="
            << run.final_planet.redox.oxidized_fe0_mass_kg << std::endl;
  std::cout << "simple_run_final_atm_C_kg="
            << run.final_planet.atmosphere.elements_kg.carbon_kg << std::endl;
  std::cout << "simple_run_final_atm_H_kg="
            << run.final_planet.atmosphere.elements_kg.hydrogen_kg << std::endl;
  std::cout << "simple_run_final_atm_N_kg="
            << run.final_planet.atmosphere.elements_kg.nitrogen_kg << std::endl;
  std::cout << "simple_run_final_atm_species_C="
            << ToString(run.final_planet.atmosphere.carbon_species) << std::endl;
  std::cout << "simple_run_final_atm_species_H="
            << ToString(run.final_planet.atmosphere.hydrogen_species) << std::endl;
  std::cout << "simple_run_final_atm_species_N="
            << ToString(run.final_planet.atmosphere.nitrogen_species) << std::endl;
  std::cout << "simple_run_cumulative_eroded_C_kg="
            << run.cumulative_eroded_kg.carbon_kg << std::endl;
  std::cout << "simple_run_cumulative_eroded_H_kg="
            << run.cumulative_eroded_kg.hydrogen_kg << std::endl;
  std::cout << "simple_run_cumulative_eroded_N_kg="
            << run.cumulative_eroded_kg.nitrogen_kg << std::endl;
  std::cout << "simple_run_max_abs_conservation_error_C_kg="
            << run.max_abs_conservation_error_kg.carbon_kg << std::endl;
  std::cout << "simple_run_max_abs_conservation_error_H_kg="
            << run.max_abs_conservation_error_kg.hydrogen_kg << std::endl;
  std::cout << "simple_run_max_abs_conservation_error_N_kg="
            << run.max_abs_conservation_error_kg.nitrogen_kg << std::endl;
  std::cout << "simple_run_first_step_error_max_kg="
            << MaxComponentAbs(run.first_step.conservation_error_kg) << std::endl;
  std::cout << "simple_run_last_step_error_max_kg="
            << MaxComponentAbs(run.last_step.conservation_error_kg) << std::endl;
}

}  // namespace self_oxidation

int main(int argc, char** argv) {
  try {
    const self_oxidation::CommandLineOptions options =
        self_oxidation::ParseCommandLine(argc, argv);

    const self_oxidation::Constants constants{};
    constants.representative_ec_alloy.Validate();

    const self_oxidation::PlanetState initial_planet =
        self_oxidation::MakeInitialPlanetState(constants);
    const std::vector<self_oxidation::Impactor> scenario =
        self_oxidation::BuildScenario(constants);
    const self_oxidation::SimplifiedRunResult simple_run =
        self_oxidation::RunScenarioSimplified(initial_planet, scenario, constants);

    if (!options.history_csv_path.empty()) {
      self_oxidation::WriteHistoryCsv(options.history_csv_path, simple_run);
    }

    if (!options.quiet) {
      self_oxidation::PrintInitialSummary(constants, initial_planet);
      self_oxidation::PrintScenarioSummary(constants, initial_planet, scenario);
      self_oxidation::PrintMeltPoolSummary(simple_run.first_step.resolved_impactor,
                                           simple_run.first_step.melt_pool);
      self_oxidation::PrintSimplifiedRunSummary(constants, simple_run);
      if (!options.history_csv_path.empty()) {
        std::cout << "history_csv=" << options.history_csv_path << std::endl;
      }
    }
    return 0;
  } catch (const std::exception& exception) {
    std::cerr << "error: " << exception.what() << '\n';
    return 1;
  }
}
