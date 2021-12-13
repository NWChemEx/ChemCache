/*
 * This file is autogenerated by: generate_basis.py 
 * 
 * NOTE: Any modifications made in this file will be lost next time 
 *       generate_basis.py is run.
 */

#include "../basis_set_list.hpp"
# include <chemist/chemist.hpp>

namespace chemcache::basis_sets {

void load_aug_dash_cc_dash_pvtz_dash_optri(chemist::BasisSetManager& bsm) {
    chemist::BasisSetManager::ao_basis_map basis_map;

    // Z = 1
    basis_map.emplace(1, chemist::AtomicBasisSet<double>("aug-cc-pvtz-optri", 1, 0.0, 0.0, 0.0));

    basis_map.at(1).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.7470550000e+00});
    basis_map.at(1).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.8882500000e-01});
    basis_map.at(1).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.0767400000e-01});
    basis_map.at(1).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.8370960000e+00});
    basis_map.at(1).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.1141910000e+00});
    basis_map.at(1).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{7.0262600000e-01});
    basis_map.at(1).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.7849140000e+00});
    basis_map.at(1).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{6.1751300000e-01});
    basis_map.at(1).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.9187250000e+00});

    // Z = 2
    basis_map.emplace(2, chemist::AtomicBasisSet<double>("aug-cc-pvtz-optri", 2, 0.0, 0.0, 0.0));

    basis_map.at(2).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.7766510000e+00});
    basis_map.at(2).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.0273960000e+00});
    basis_map.at(2).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.1240600000e-01});
    basis_map.at(2).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{7.7011780000e+00});
    basis_map.at(2).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.0313920000e+00});
    basis_map.at(2).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.1906500000e-01});
    basis_map.at(2).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.5412270000e+00});
    basis_map.at(2).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.0478080000e+00});
    basis_map.at(2).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.4275120000e+00});

    // Z = 5
    basis_map.emplace(5, chemist::AtomicBasisSet<double>("aug-cc-pvtz-optri", 5, 0.0, 0.0, 0.0));

    basis_map.at(5).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.2671480000e+00});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{8.8007900000e-01});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.4478800000e-01});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.6164000000e-02});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.4841420000e+00});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{7.5285500000e-01});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.1907400000e-01});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.1560800000e-01});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.0116000000e-02});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.8473750000e+00});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.3974980000e+00});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.9366400000e-01});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{9.0864000000e-02});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.8840320000e+00});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{8.7878000000e-01});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.7007100000e-01});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 4,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.0223720000e+00});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 4,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.7533100000e-01});

    // Z = 6
    basis_map.emplace(6, chemist::AtomicBasisSet<double>("aug-cc-pvtz-optri", 6, 0.0, 0.0, 0.0));

    basis_map.at(6).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.9642050000e+00});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.3508990000e+00});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.4716000000e-01});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{7.1296000000e-02});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.3634370000e+00});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.0793800000e+00});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{5.7198100000e-01});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.8055100000e-01});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{5.3629000000e-02});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.9484350000e+00});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.0420780000e+00});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.9920500000e-01});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.4973800000e-01});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.4897200000e+00});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.3604080000e+00});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{5.0696900000e-01});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 4,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.6030840000e+00});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 4,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{5.8323100000e-01});

    // Z = 7
    basis_map.emplace(7, chemist::AtomicBasisSet<double>("aug-cc-pvtz-optri", 7, 0.0, 0.0, 0.0));

    basis_map.at(7).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{5.8477420000e+00});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.0001550000e+00});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.9413800000e-01});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{8.6385000000e-02});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.4463680000e+00});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.4548940000e+00});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{8.8508700000e-01});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.5803500000e-01});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{7.3452000000e-02});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{7.9836110000e+00});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.4268080000e+00});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{7.9084400000e-01});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.3959900000e-01});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{5.4040350000e+00});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.0241770000e+00});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{7.2907100000e-01});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 4,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.5286770000e+00});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 4,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{9.5053900000e-01});

    // Z = 8
    basis_map.emplace(8, chemist::AtomicBasisSet<double>("aug-cc-pvtz-optri", 8, 0.0, 0.0, 0.0));

    basis_map.at(8).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.1352581000e+01});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.3758930000e+00});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.5330200000e-01});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.1043700000e-01});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{5.4560910000e+00});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.3525840000e+00});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.0679520000e+00});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.3014500000e-01});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{8.9368000000e-02});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.2069849000e+01});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.2787000000e+00});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.3747810000e+00});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.1932700000e-01});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{6.9152510000e+00});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.2424930000e+00});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{7.4993500000e-01});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 4,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.1413190000e+00});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 4,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.0754710000e+00});

    // Z = 9
    basis_map.emplace(9, chemist::AtomicBasisSet<double>("aug-cc-pvtz-optri", 9, 0.0, 0.0, 0.0));

    basis_map.at(9).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.6249086000e+01});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.3814200000e+00});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.0996880000e+00});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.3730100000e-01});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{7.6657940000e+00});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.0196150000e+00});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.3629930000e+00});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.0423300000e-01});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.5171700000e-01});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.6217768000e+01});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{5.6999140000e+00});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.0227140000e+00});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.3819500000e-01});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{9.8128910000e+00});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.8118310000e+00});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.0807210000e+00});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 4,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.0534990000e+00});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 4,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.3579070000e+00});

    // Z = 10
    basis_map.emplace(10, chemist::AtomicBasisSet<double>("aug-cc-pvtz-optri", 10, 0.0, 0.0, 0.0));

    basis_map.at(10).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.1763868000e+01});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.4564600000e+00});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.4799290000e+00});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.2367400000e-01});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{8.8150200000e+00});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.5110110000e+00});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.7043880000e+00});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{5.4583500000e-01});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.1942800000e-01});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.6026006000e+01});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.1059234000e+01});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.7659490000e+00});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{6.3895300000e-01});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.3036287000e+01});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{7.6031790000e+00});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.6025530000e+00});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 4,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{5.1878530000e+00});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 4,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.7994760000e+00});

    // Z = 13
    basis_map.emplace(13, chemist::AtomicBasisSet<double>("aug-cc-pvtz-optri", 13, 0.0, 0.0, 0.0));

    basis_map.at(13).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.9522820000e+00});
    basis_map.at(13).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.9926600000e-01});
    basis_map.at(13).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.1361300000e-01});
    basis_map.at(13).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.4202000000e-02});
    basis_map.at(13).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.3226840000e+00});
    basis_map.at(13).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.9794900000e+00});
    basis_map.at(13).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.3229000000e+00});
    basis_map.at(13).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.3870200000e-01});
    basis_map.at(13).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.0589000000e-02});
    basis_map.at(13).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{9.2862280000e+00});
    basis_map.at(13).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.3523380000e+00});
    basis_map.at(13).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.9881800000e-01});
    basis_map.at(13).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.6639200000e-01});
    basis_map.at(13).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.3241760000e+00});
    basis_map.at(13).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.1097000000e-01});
    basis_map.at(13).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.5064000000e-01});
    basis_map.at(13).add_shell(chemist::ShellType::pure, 4,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.8185500000e-01});
    basis_map.at(13).add_shell(chemist::ShellType::pure, 4,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.2591700000e-01});

    // Z = 14
    basis_map.emplace(14, chemist::AtomicBasisSet<double>("aug-cc-pvtz-optri", 14, 0.0, 0.0, 0.0));

    basis_map.at(14).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{6.5782740000e+00});
    basis_map.at(14).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{7.0173700000e-01});
    basis_map.at(14).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.6770000000e-01});
    basis_map.at(14).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{5.1182000000e-02});
    basis_map.at(14).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{7.2071580000e+00});
    basis_map.at(14).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.1878490000e+00});
    basis_map.at(14).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.3682450000e+00});
    basis_map.at(14).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.8343600000e-01});
    basis_map.at(14).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.9935000000e-02});
    basis_map.at(14).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.0065409000e+01});
    basis_map.at(14).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.1416080000e+00});
    basis_map.at(14).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{9.1882100000e-01});
    basis_map.at(14).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.5590900000e-01});
    basis_map.at(14).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.9190560000e+00});
    basis_map.at(14).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{5.1551200000e-01});
    basis_map.at(14).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.2441700000e-01});
    basis_map.at(14).add_shell(chemist::ShellType::pure, 4,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{6.5190900000e-01});
    basis_map.at(14).add_shell(chemist::ShellType::pure, 4,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.1767000000e-01});

    // Z = 15
    basis_map.emplace(15, chemist::AtomicBasisSet<double>("aug-cc-pvtz-optri", 15, 0.0, 0.0, 0.0));

    basis_map.at(15).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{8.2324400000e+00});
    basis_map.at(15).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.4034300000e+00});
    basis_map.at(15).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.1018700000e-01});
    basis_map.at(15).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{6.6043000000e-02});
    basis_map.at(15).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{9.6778770000e+00});
    basis_map.at(15).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.5253700000e+00});
    basis_map.at(15).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.2013060000e+00});
    basis_map.at(15).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.9667000000e-01});
    basis_map.at(15).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.5997000000e-02});
    basis_map.at(15).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.3648659000e+01});
    basis_map.at(15).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.6678810000e+00});
    basis_map.at(15).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.0793380000e+00});
    basis_map.at(15).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.2669500000e-01});
    basis_map.at(15).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.2525360000e+00});
    basis_map.at(15).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{6.8876600000e-01});
    basis_map.at(15).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.0198000000e-01});
    basis_map.at(15).add_shell(chemist::ShellType::pure, 4,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{8.4846700000e-01});
    basis_map.at(15).add_shell(chemist::ShellType::pure, 4,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.0373300000e-01});

    // Z = 16
    basis_map.emplace(16, chemist::AtomicBasisSet<double>("aug-cc-pvtz-optri", 16, 0.0, 0.0, 0.0));

    basis_map.at(16).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{7.1544480000e+00});
    basis_map.at(16).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.7056790000e+00});
    basis_map.at(16).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.8796500000e-01});
    basis_map.at(16).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{7.9949000000e-02});
    basis_map.at(16).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.2343066000e+01});
    basis_map.at(16).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.4969470000e+00});
    basis_map.at(16).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.6303970000e+00});
    basis_map.at(16).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{5.3941300000e-01});
    basis_map.at(16).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{5.6470000000e-02});
    basis_map.at(16).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.9290141000e+01});
    basis_map.at(16).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{7.3717540000e+00});
    basis_map.at(16).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.8762860000e+00});
    basis_map.at(16).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.3652200000e-01});
    basis_map.at(16).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{5.4398970000e+00});
    basis_map.at(16).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{8.5729800000e-01});
    basis_map.at(16).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.7236000000e-01});
    basis_map.at(16).add_shell(chemist::ShellType::pure, 4,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.0512970000e+00});
    basis_map.at(16).add_shell(chemist::ShellType::pure, 4,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.8717400000e-01});

    // Z = 17
    basis_map.emplace(17, chemist::AtomicBasisSet<double>("aug-cc-pvtz-optri", 17, 0.0, 0.0, 0.0));

    basis_map.at(17).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{7.6159120000e+00});
    basis_map.at(17).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.7009310000e+00});
    basis_map.at(17).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.7369500000e-01});
    basis_map.at(17).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.0849800000e-01});
    basis_map.at(17).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.3838955000e+01});
    basis_map.at(17).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.1705200000e+00});
    basis_map.at(17).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.9624920000e+00});
    basis_map.at(17).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{5.5562800000e-01});
    basis_map.at(17).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{8.4739000000e-02});
    basis_map.at(17).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.4875052000e+01});
    basis_map.at(17).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.0080369000e+01});
    basis_map.at(17).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.0236740000e+00});
    basis_map.at(17).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{5.0847400000e-01});
    basis_map.at(17).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{6.8531320000e+00});
    basis_map.at(17).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.0625280000e+00});
    basis_map.at(17).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.7185200000e-01});
    basis_map.at(17).add_shell(chemist::ShellType::pure, 4,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.2891180000e+00});
    basis_map.at(17).add_shell(chemist::ShellType::pure, 4,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{5.9232000000e-01});

    // Z = 18
    basis_map.emplace(18, chemist::AtomicBasisSet<double>("aug-cc-pvtz-optri", 18, 0.0, 0.0, 0.0));

    basis_map.at(18).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{8.6144650000e+00});
    basis_map.at(18).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.8927430000e+00});
    basis_map.at(18).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{5.6983700000e-01});
    basis_map.at(18).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.3065600000e-01});
    basis_map.at(18).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.6824368000e+01});
    basis_map.at(18).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{6.3006820000e+00});
    basis_map.at(18).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.8069490000e+00});
    basis_map.at(18).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{5.3505100000e-01});
    basis_map.at(18).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{8.2595000000e-02});
    basis_map.at(18).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.1255380000e+01});
    basis_map.at(18).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.2449987000e+01});
    basis_map.at(18).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.3563880000e+00});
    basis_map.at(18).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{6.4288400000e-01});
    basis_map.at(18).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{8.1895310000e+00});
    basis_map.at(18).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.3323360000e+00});
    basis_map.at(18).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{6.0262800000e-01});
    basis_map.at(18).add_shell(chemist::ShellType::pure, 4,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.5693630000e+00});
    basis_map.at(18).add_shell(chemist::ShellType::pure, 4,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{7.4704700000e-01});


    bsm.insert("aug-cc-pvtz-optri", basis_map);
} // function load_aug_dash_cc_dash_pvtz_dash_optri

} // namespace chemcache::basis_sets
