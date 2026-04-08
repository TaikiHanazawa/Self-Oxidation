#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
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
  double target_melt_volume_m3 = 0.0;
  double target_melt_mass_kg = 0.0;
  double impactor_silicate_mass_kg = 0.0;
  double impactor_metal_mass_kg = 0.0;
  double equilibrating_silicate_mass_kg = 0.0;
  double depth_m = 0.0;
  double pressure_gpa = 0.0;
  double temperature_k = 0.0;
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
  PlanetState planet_after{};
  MeltPoolState melt_pool{};
  double dn_metal_silicate = std::numeric_limits<double>::quiet_NaN();
  double dh_metal_silicate = std::numeric_limits<double>::quiet_NaN();
  double dc_metal_silicate = std::numeric_limits<double>::quiet_NaN();
};

struct Constants {
  double earth_mass_kg = 5.9722e24;
  double earth_radius_m = 6.371e6;
  double gravitational_constant_si = 6.67430e-11;
  double mantle_density_kg_m3 = 3300.0;
  double core_density_kg_m3 = 8000.0;

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
  std::cout << "Scaffold ready\n";
  std::cout << "initial_planet_mass_earth="
            << KgToEarthMass(planet.total_mass_kg(), constants) << '\n';
  std::cout << "mantle_mass_kg=" << planet.mantle.mass_kg << '\n';
  std::cout << "core_mass_kg=" << planet.core.mass_kg << '\n';
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

}  // namespace self_oxidation

int main() {
  try {
    const self_oxidation::Constants constants{};
    constants.representative_ec_alloy.Validate();

    const self_oxidation::PlanetState initial_planet =
        self_oxidation::MakeInitialPlanetState(constants);
    const std::vector<self_oxidation::Impactor> scenario =
        self_oxidation::BuildScenario(constants);

    self_oxidation::PrintInitialSummary(constants, initial_planet);
    self_oxidation::PrintScenarioSummary(constants, initial_planet, scenario);
    return 0;
  } catch (const std::exception& exception) {
    std::cerr << "error: " << exception.what() << '\n';
    return 1;
  }
}
