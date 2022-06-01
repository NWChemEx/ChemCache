/*
 * This file is autogenerated by: generate_basis.py
 *
 * NOTE: Any modifications made in this file will be lost next time
 *       generate_basis.py is run.
 */

#include "../basis_set_list.hpp"
#include <chemist/chemist.hpp>

namespace chemcache::basis_sets {

void load_aug_dash_cc_dash_pvdz_dash_pp_dash_optri(
  chemist::BasisSetManager& bsm) {
    chemist::BasisSetManager::ao_basis_map basis_map;

    // Z = 29
    basis_map.emplace(29, chemist::AtomicBasisSet<double>(
                            "aug-cc-pvdz-pp-optri", 29, 0.0, 0.0, 0.0));

    basis_map.at(29).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.3615150000e+01});
    basis_map.at(29).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{7.4436000000e+00});
    basis_map.at(29).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{4.6986000000e+00});
    basis_map.at(29).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{2.3742140000e+00});
    basis_map.at(29).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.4483190000e+00});
    basis_map.at(29).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{3.4686000000e-01});
    basis_map.at(29).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.1092380000e+01});
    basis_map.at(29).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{7.3506280000e+00});
    basis_map.at(29).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{4.9315000000e+00});
    basis_map.at(29).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{2.0277600000e+00});
    basis_map.at(29).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{7.8788100000e-01});
    basis_map.at(29).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{4.3122300000e-01});
    basis_map.at(29).add_shell(chemist::ShellType::pure, 2,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.5688530000e+01});
    basis_map.at(29).add_shell(chemist::ShellType::pure, 2,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{6.4603340000e+00});
    basis_map.at(29).add_shell(chemist::ShellType::pure, 2,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{4.1723270000e+00});
    basis_map.at(29).add_shell(chemist::ShellType::pure, 2,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{2.7920000000e+00});
    basis_map.at(29).add_shell(chemist::ShellType::pure, 2,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.0688050000e+00});
    basis_map.at(29).add_shell(chemist::ShellType::pure, 2,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{7.1144700000e-01});
    basis_map.at(29).add_shell(chemist::ShellType::pure, 3,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{2.1890890000e+01});
    basis_map.at(29).add_shell(chemist::ShellType::pure, 3,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{8.3307000000e+00});
    basis_map.at(29).add_shell(chemist::ShellType::pure, 3,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{5.5693780000e+00});
    basis_map.at(29).add_shell(chemist::ShellType::pure, 3,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.3523440000e+00});
    basis_map.at(29).add_shell(chemist::ShellType::pure, 3,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{3.9762600000e-01});
    basis_map.at(29).add_shell(chemist::ShellType::pure, 4,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.4168310000e+01});
    basis_map.at(29).add_shell(chemist::ShellType::pure, 4,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{5.3687560000e+00});
    basis_map.at(29).add_shell(chemist::ShellType::pure, 4,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{2.0874870000e+00});
    basis_map.at(29).add_shell(chemist::ShellType::pure, 4,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{7.5582100000e-01});
    basis_map.at(29).add_shell(chemist::ShellType::pure, 5,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{8.0312760000e+00});
    basis_map.at(29).add_shell(chemist::ShellType::pure, 5,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{2.6281350000e+00});
    basis_map.at(29).add_shell(chemist::ShellType::pure, 6,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{5.4472130000e+00});

    // Z = 30
    basis_map.emplace(30, chemist::AtomicBasisSet<double>(
                            "aug-cc-pvdz-pp-optri", 30, 0.0, 0.0, 0.0));

    basis_map.at(30).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.3700190000e+01});
    basis_map.at(30).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{9.1444660000e+00});
    basis_map.at(30).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{5.0733000000e+00});
    basis_map.at(30).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{2.2237420000e+00});
    basis_map.at(30).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.4850770000e+00});
    basis_map.at(30).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{4.1816000000e-01});
    basis_map.at(30).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{2.3262550000e+01});
    basis_map.at(30).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.5673670000e+01});
    basis_map.at(30).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{5.3705000000e+00});
    basis_map.at(30).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.9247800000e+00});
    basis_map.at(30).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{8.4626400000e-01});
    basis_map.at(30).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{5.3571500000e-01});
    basis_map.at(30).add_shell(chemist::ShellType::pure, 2,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.8613380000e+01});
    basis_map.at(30).add_shell(chemist::ShellType::pure, 2,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{7.5009300000e+00});
    basis_map.at(30).add_shell(chemist::ShellType::pure, 2,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{4.9823310000e+00});
    basis_map.at(30).add_shell(chemist::ShellType::pure, 2,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{3.3188000000e+00});
    basis_map.at(30).add_shell(chemist::ShellType::pure, 2,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.6151700000e+00});
    basis_map.at(30).add_shell(chemist::ShellType::pure, 2,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.0748830000e+00});
    basis_map.at(30).add_shell(chemist::ShellType::pure, 3,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{2.5132730000e+01});
    basis_map.at(30).add_shell(chemist::ShellType::pure, 3,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{9.5199000000e+00});
    basis_map.at(30).add_shell(chemist::ShellType::pure, 3,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{6.3711850000e+00});
    basis_map.at(30).add_shell(chemist::ShellType::pure, 3,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.6204940000e+00});
    basis_map.at(30).add_shell(chemist::ShellType::pure, 3,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{4.5101300000e-01});
    basis_map.at(30).add_shell(chemist::ShellType::pure, 4,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.7163670000e+01});
    basis_map.at(30).add_shell(chemist::ShellType::pure, 4,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{6.4469590000e+00});
    basis_map.at(30).add_shell(chemist::ShellType::pure, 4,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{2.6776790000e+00});
    basis_map.at(30).add_shell(chemist::ShellType::pure, 4,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.0182460000e+00});
    basis_map.at(30).add_shell(chemist::ShellType::pure, 5,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{9.1262860000e+00});
    basis_map.at(30).add_shell(chemist::ShellType::pure, 5,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{3.0970450000e+00});
    basis_map.at(30).add_shell(chemist::ShellType::pure, 6,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{6.4053200000e+00});

    // Z = 47
    basis_map.emplace(47, chemist::AtomicBasisSet<double>(
                            "aug-cc-pvdz-pp-optri", 47, 0.0, 0.0, 0.0));

    basis_map.at(47).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.2639210000e+01});
    basis_map.at(47).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{8.5539920000e+00});
    basis_map.at(47).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{2.3985000000e+00});
    basis_map.at(47).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.6093170000e+00});
    basis_map.at(47).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{8.9113700000e-01});
    basis_map.at(47).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{5.9950200000e-01});
    basis_map.at(47).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.1469780000e+01});
    basis_map.at(47).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{7.7147390000e+00});
    basis_map.at(47).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{5.1799100000e+00});
    basis_map.at(47).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{2.5852000000e+00});
    basis_map.at(47).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.7497360000e+00});
    basis_map.at(47).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{5.2696900000e-01});
    basis_map.at(47).add_shell(chemist::ShellType::pure, 2,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{2.1001460000e+01});
    basis_map.at(47).add_shell(chemist::ShellType::pure, 2,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{6.3441670000e+00});
    basis_map.at(47).add_shell(chemist::ShellType::pure, 2,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{4.2406120000e+00});
    basis_map.at(47).add_shell(chemist::ShellType::pure, 2,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{2.8365000000e+00});
    basis_map.at(47).add_shell(chemist::ShellType::pure, 2,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{8.3016400000e-01});
    basis_map.at(47).add_shell(chemist::ShellType::pure, 2,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{5.5366500000e-01});
    basis_map.at(47).add_shell(chemist::ShellType::pure, 3,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.1659440000e+01});
    basis_map.at(47).add_shell(chemist::ShellType::pure, 3,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{7.3985980000e+00});
    basis_map.at(47).add_shell(chemist::ShellType::pure, 3,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{3.9751000000e+00});
    basis_map.at(47).add_shell(chemist::ShellType::pure, 3,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{2.6532840000e+00});
    basis_map.at(47).add_shell(chemist::ShellType::pure, 3,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{6.6316400000e-01});
    basis_map.at(47).add_shell(chemist::ShellType::pure, 4,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{7.0187400000e+00});
    basis_map.at(47).add_shell(chemist::ShellType::pure, 4,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{3.0883250000e+00});
    basis_map.at(47).add_shell(chemist::ShellType::pure, 4,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.3945630000e+00});
    basis_map.at(47).add_shell(chemist::ShellType::pure, 4,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{5.8335200000e-01});
    basis_map.at(47).add_shell(chemist::ShellType::pure, 5,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{3.1015640000e+00});
    basis_map.at(47).add_shell(chemist::ShellType::pure, 5,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.3194290000e+00});
    basis_map.at(47).add_shell(chemist::ShellType::pure, 6,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{2.2667700000e+00});

    // Z = 48
    basis_map.emplace(48, chemist::AtomicBasisSet<double>(
                            "aug-cc-pvdz-pp-optri", 48, 0.0, 0.0, 0.0));

    basis_map.at(48).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.2869320000e+01});
    basis_map.at(48).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{6.3258540000e+00});
    basis_map.at(48).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{2.6738000000e+00});
    basis_map.at(48).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.5123050000e+00});
    basis_map.at(48).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{6.9460800000e-01});
    basis_map.at(48).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{3.4276000000e-01});
    basis_map.at(48).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.2571740000e+01});
    basis_map.at(48).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{8.5174520000e+00});
    basis_map.at(48).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{5.7467910000e+00});
    basis_map.at(48).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{2.7157000000e+00});
    basis_map.at(48).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.8560780000e+00});
    basis_map.at(48).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{5.2131100000e-01});
    basis_map.at(48).add_shell(chemist::ShellType::pure, 2,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{2.1832180000e+01});
    basis_map.at(48).add_shell(chemist::ShellType::pure, 2,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{7.3798510000e+00});
    basis_map.at(48).add_shell(chemist::ShellType::pure, 2,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{4.9417370000e+00});
    basis_map.at(48).add_shell(chemist::ShellType::pure, 2,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{3.3125000000e+00});
    basis_map.at(48).add_shell(chemist::ShellType::pure, 2,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{2.0469140000e+00});
    basis_map.at(48).add_shell(chemist::ShellType::pure, 2,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{6.5284000000e-01});
    basis_map.at(48).add_shell(chemist::ShellType::pure, 3,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.5377210000e+01});
    basis_map.at(48).add_shell(chemist::ShellType::pure, 3,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{7.8405300000e+00});
    basis_map.at(48).add_shell(chemist::ShellType::pure, 3,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{4.5765000000e+00});
    basis_map.at(48).add_shell(chemist::ShellType::pure, 3,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{2.9481920000e+00});
    basis_map.at(48).add_shell(chemist::ShellType::pure, 3,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{7.7311600000e-01});
    basis_map.at(48).add_shell(chemist::ShellType::pure, 4,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{7.1320870000e+00});
    basis_map.at(48).add_shell(chemist::ShellType::pure, 4,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{3.4520730000e+00});
    basis_map.at(48).add_shell(chemist::ShellType::pure, 4,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.5944050000e+00});
    basis_map.at(48).add_shell(chemist::ShellType::pure, 4,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{6.8856900000e-01});
    basis_map.at(48).add_shell(chemist::ShellType::pure, 5,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{3.4499030000e+00});
    basis_map.at(48).add_shell(chemist::ShellType::pure, 5,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.5154460000e+00});
    basis_map.at(48).add_shell(chemist::ShellType::pure, 6,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{2.5575350000e+00});

    // Z = 79
    basis_map.emplace(79, chemist::AtomicBasisSet<double>(
                            "aug-cc-pvdz-pp-optri", 79, 0.0, 0.0, 0.0));

    basis_map.at(79).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.5054590000e+01});
    basis_map.at(79).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.0041750000e+01});
    basis_map.at(79).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{4.7854040000e+00});
    basis_map.at(79).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{2.0643000000e+00});
    basis_map.at(79).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.2661230000e+00});
    basis_map.at(79).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{8.4507100000e-01});
    basis_map.at(79).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{2.7283290000e+01});
    basis_map.at(79).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{9.0880680000e+00});
    basis_map.at(79).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{6.0522510000e+00});
    basis_map.at(79).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{2.2606000000e+00});
    basis_map.at(79).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.5092490000e+00});
    basis_map.at(79).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{7.7702100000e-01});
    basis_map.at(79).add_shell(chemist::ShellType::pure, 2,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.0167350000e+01});
    basis_map.at(79).add_shell(chemist::ShellType::pure, 2,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{6.8013230000e+00});
    basis_map.at(79).add_shell(chemist::ShellType::pure, 2,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{4.5419460000e+00});
    basis_map.at(79).add_shell(chemist::ShellType::pure, 2,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{2.0167000000e+00});
    basis_map.at(79).add_shell(chemist::ShellType::pure, 2,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.3564320000e+00});
    basis_map.at(79).add_shell(chemist::ShellType::pure, 2,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{5.0969300000e-01});
    basis_map.at(79).add_shell(chemist::ShellType::pure, 3,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{7.3994110000e+00});
    basis_map.at(79).add_shell(chemist::ShellType::pure, 3,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{4.6716510000e+00});
    basis_map.at(79).add_shell(chemist::ShellType::pure, 3,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{3.1253190000e+00});
    basis_map.at(79).add_shell(chemist::ShellType::pure, 3,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{2.0995000000e+00});
    basis_map.at(79).add_shell(chemist::ShellType::pure, 3,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{5.6652600000e-01});
    basis_map.at(79).add_shell(chemist::ShellType::pure, 4,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{4.4745740000e+00});
    basis_map.at(79).add_shell(chemist::ShellType::pure, 4,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{2.9860180000e+00});
    basis_map.at(79).add_shell(chemist::ShellType::pure, 4,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.3999910000e+00});
    basis_map.at(79).add_shell(chemist::ShellType::pure, 4,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{5.7447100000e-01});
    basis_map.at(79).add_shell(chemist::ShellType::pure, 5,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{2.9044630000e+00});
    basis_map.at(79).add_shell(chemist::ShellType::pure, 5,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.2824050000e+00});
    basis_map.at(79).add_shell(chemist::ShellType::pure, 6,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.7224630000e+00});

    // Z = 80
    basis_map.emplace(80, chemist::AtomicBasisSet<double>(
                            "aug-cc-pvdz-pp-optri", 80, 0.0, 0.0, 0.0));

    basis_map.at(80).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.3373940000e+01});
    basis_map.at(80).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{8.9239300000e+00});
    basis_map.at(80).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{5.7967510000e+00});
    basis_map.at(80).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{2.3359000000e+00});
    basis_map.at(80).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.4314090000e+00});
    basis_map.at(80).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{4.2733300000e-01});
    basis_map.at(80).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{2.2115300000e+01});
    basis_map.at(80).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.1483160000e+01});
    basis_map.at(80).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{7.6630510000e+00});
    basis_map.at(80).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{2.3585000000e+00});
    basis_map.at(80).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.4580420000e+00});
    basis_map.at(80).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{7.4860000000e-01});
    basis_map.at(80).add_shell(chemist::ShellType::pure, 2,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.0679550000e+01});
    basis_map.at(80).add_shell(chemist::ShellType::pure, 2,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{7.1339590000e+00});
    basis_map.at(80).add_shell(chemist::ShellType::pure, 2,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{4.7504670000e+00});
    basis_map.at(80).add_shell(chemist::ShellType::pure, 2,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{2.2228000000e+00});
    basis_map.at(80).add_shell(chemist::ShellType::pure, 2,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.4854680000e+00});
    basis_map.at(80).add_shell(chemist::ShellType::pure, 2,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{3.9450900000e-01});
    basis_map.at(80).add_shell(chemist::ShellType::pure, 3,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{8.3062810000e+00});
    basis_map.at(80).add_shell(chemist::ShellType::pure, 3,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{5.5539610000e+00});
    basis_map.at(80).add_shell(chemist::ShellType::pure, 3,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{3.7121290000e+00});
    basis_map.at(80).add_shell(chemist::ShellType::pure, 3,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{2.2727000000e+00});
    basis_map.at(80).add_shell(chemist::ShellType::pure, 3,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{6.5591400000e-01});
    basis_map.at(80).add_shell(chemist::ShellType::pure, 4,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{4.8929370000e+00});
    basis_map.at(80).add_shell(chemist::ShellType::pure, 4,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{3.2679300000e+00});
    basis_map.at(80).add_shell(chemist::ShellType::pure, 4,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.5760520000e+00});
    basis_map.at(80).add_shell(chemist::ShellType::pure, 4,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{6.6214800000e-01});
    basis_map.at(80).add_shell(chemist::ShellType::pure, 5,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{2.9810210000e+00});
    basis_map.at(80).add_shell(chemist::ShellType::pure, 5,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.2624060000e+00});
    basis_map.at(80).add_shell(chemist::ShellType::pure, 6,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.8808450000e+00});

    bsm.insert("aug-cc-pvdz-pp-optri", basis_map);
} // function load_aug_dash_cc_dash_pvdz_dash_pp_dash_optri

} // namespace chemcache::basis_sets
