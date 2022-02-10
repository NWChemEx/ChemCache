/*
 * This file is autogenerated by: generate_basis.py 
 * 
 * NOTE: Any modifications made in this file will be lost next time 
 *       generate_basis.py is run.
 */

#include "../basis_set_list.hpp"
# include <chemist/chemist.hpp>

namespace chemcache::basis_sets {

void load_d_dash_aug_dash_cc_dash_pvtz(chemist::BasisSetManager& bsm) {
    chemist::BasisSetManager::ao_basis_map basis_map;

    // Z = 1
    basis_map.emplace(1, chemist::AtomicBasisSet<double>("d-aug-cc-pvtz", 1, 0.0, 0.0, 0.0));

    basis_map.at(1).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.2580000000e-01});
    basis_map.at(1).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{6.0680000000e-03,4.5308000000e-02,2.0282200000e-01,5.0390300000e-01,3.8342100000e-01},
        std::vector<double>{3.3870000000e+01,5.0950000000e+00,1.1590000000e+00,3.2580000000e-01,1.0270000000e-01});
    basis_map.at(1).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.0270000000e-01});
    basis_map.at(1).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.5260000000e-02});
    basis_map.at(1).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{6.2100000000e-03});
    basis_map.at(1).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.4070000000e+00});
    basis_map.at(1).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.8800000000e-01});
    basis_map.at(1).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.0200000000e-01});
    basis_map.at(1).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.6800000000e-02});
    basis_map.at(1).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.0570000000e+00});
    basis_map.at(1).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.4700000000e-01});
    basis_map.at(1).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{5.7700000000e-02});

    // Z = 2
    basis_map.emplace(2, chemist::AtomicBasisSet<double>("d-aug-cc-pvtz", 2, 0.0, 0.0, 0.0));

    basis_map.at(2).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{6.6690000000e-01});
    basis_map.at(2).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{2.5870000000e-03,1.9533000000e-02,9.0998000000e-02,2.7205000000e-01,4.7806500000e-01,3.0773700000e-01},
        std::vector<double>{2.3400000000e+02,3.5160000000e+01,7.9890000000e+00,2.2120000000e+00,6.6690000000e-01,2.0890000000e-01});
    basis_map.at(2).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.0890000000e-01});
    basis_map.at(2).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{5.1380000000e-02});
    basis_map.at(2).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.2600000000e-02});
    basis_map.at(2).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.0440000000e+00});
    basis_map.at(2).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{7.5800000000e-01});
    basis_map.at(2).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.9930000000e-01});
    basis_map.at(2).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{5.2400000000e-02});
    basis_map.at(2).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.9650000000e+00});
    basis_map.at(2).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.5920000000e-01});
    basis_map.at(2).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.0700000000e-01});

    // Z = 5
    basis_map.emplace(5, chemist::AtomicBasisSet<double>("d-aug-cc-pvtz", 5, 0.0, 0.0, 0.0));

    basis_map.at(5).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{5.5500000000e-04,4.2910000000e-03,2.1949000000e-02,8.4441000000e-02,2.3855700000e-01,4.3507200000e-01,3.4195500000e-01,3.6856000000e-02,-9.5450000000e-03,2.3680000000e-03},
        std::vector<double>{5.4730000000e+03,8.2090000000e+02,1.8680000000e+02,5.2830000000e+01,1.7080000000e+01,5.9990000000e+00,2.2080000000e+00,5.8790000000e-01,2.4150000000e-01,8.6100000000e-02});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{5.8790000000e-01});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{-1.1200000000e-04,-8.6800000000e-04,-4.4840000000e-03,-1.7683000000e-02,-5.3639000000e-02,-1.1900500000e-01,-1.6582400000e-01,1.2010700000e-01,5.9598100000e-01,4.1102100000e-01},
        std::vector<double>{5.4730000000e+03,8.2090000000e+02,1.8680000000e+02,5.2830000000e+01,1.7080000000e+01,5.9990000000e+00,2.2080000000e+00,5.8790000000e-01,2.4150000000e-01,8.6100000000e-02});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{8.6100000000e-02});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.9140000000e-02});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{9.8600000000e-03});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.3850000000e-01});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.3118000000e-02,7.9896000000e-02,2.7727500000e-01,5.0427000000e-01,3.5368000000e-01},
        std::vector<double>{1.2050000000e+01,2.6130000000e+00,7.4750000000e-01,2.3850000000e-01,7.6980000000e-02});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{7.6980000000e-02});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.0960000000e-02});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{5.7100000000e-03});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{6.6100000000e-01});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.9900000000e-01});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{6.0400000000e-02});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.8300000000e-02});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.9000000000e-01});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.6300000000e-01});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{5.4200000000e-02});

    // Z = 6
    basis_map.emplace(6, chemist::AtomicBasisSet<double>("d-aug-cc-pvtz", 6, 0.0, 0.0, 0.0));

    basis_map.at(6).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{5.3100000000e-04,4.1080000000e-03,2.1087000000e-02,8.1853000000e-02,2.3481700000e-01,4.3440100000e-01,3.4612900000e-01,3.9378000000e-02,-8.9830000000e-03,2.3850000000e-03},
        std::vector<double>{8.2360000000e+03,1.2350000000e+03,2.8080000000e+02,7.9270000000e+01,2.5590000000e+01,8.9970000000e+00,3.3190000000e+00,9.0590000000e-01,3.6430000000e-01,1.2850000000e-01});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{9.0590000000e-01});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{-1.1300000000e-04,-8.7800000000e-04,-4.5400000000e-03,-1.8133000000e-02,-5.5760000000e-02,-1.2689500000e-01,-1.7035200000e-01,1.4038200000e-01,5.9868400000e-01,3.9538900000e-01},
        std::vector<double>{8.2360000000e+03,1.2350000000e+03,2.8080000000e+02,7.9270000000e+01,2.5590000000e+01,8.9970000000e+00,3.3190000000e+00,9.0590000000e-01,3.6430000000e-01,1.2850000000e-01});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.2850000000e-01});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.4020000000e-02});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.5100000000e-02});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.8270000000e-01});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.4031000000e-02,8.6866000000e-02,2.9021600000e-01,5.0100800000e-01,3.4340600000e-01},
        std::vector<double>{1.8710000000e+01,4.1330000000e+00,1.2000000000e+00,3.8270000000e-01,1.2090000000e-01});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.2090000000e-01});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.5690000000e-02});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.0500000000e-02});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.0970000000e+00});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.1800000000e-01});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.0000000000e-01});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.1400000000e-02});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{7.6100000000e-01});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.6800000000e-01});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{9.4400000000e-02});

    // Z = 7
    basis_map.emplace(7, chemist::AtomicBasisSet<double>("d-aug-cc-pvtz", 7, 0.0, 0.0, 0.0));

    basis_map.at(7).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{5.2300000000e-04,4.0450000000e-03,2.0775000000e-02,8.0727000000e-02,2.3307400000e-01,4.3350100000e-01,3.4747200000e-01,4.1262000000e-02,-8.5080000000e-03,2.3840000000e-03},
        std::vector<double>{1.1420000000e+04,1.7120000000e+03,3.8930000000e+02,1.1000000000e+02,3.5570000000e+01,1.2540000000e+01,4.6440000000e+00,1.2930000000e+00,5.1180000000e-01,1.7870000000e-01});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.2930000000e+00});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{-1.1500000000e-04,-8.9500000000e-04,-4.6240000000e-03,-1.8528000000e-02,-5.7339000000e-02,-1.3207600000e-01,-1.7251000000e-01,1.5181400000e-01,5.9994400000e-01,3.8746200000e-01},
        std::vector<double>{1.1420000000e+04,1.7120000000e+03,3.8930000000e+02,1.1000000000e+02,3.5570000000e+01,1.2540000000e+01,4.6440000000e+00,1.2930000000e+00,5.1180000000e-01,1.7870000000e-01});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.7870000000e-01});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{5.7600000000e-02});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.8600000000e-02});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{5.5500000000e-01});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.4670000000e-02,9.1764000000e-02,2.9868300000e-01,4.9848700000e-01,3.3702300000e-01},
        std::vector<double>{2.6630000000e+01,5.9480000000e+00,1.7420000000e+00,5.5500000000e-01,1.7250000000e-01});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.7250000000e-01});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.9100000000e-02});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.4000000000e-02});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.6540000000e+00});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.6900000000e-01});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.5100000000e-01});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.8600000000e-02});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.0930000000e+00});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.6400000000e-01});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.2100000000e-01});

    // Z = 8
    basis_map.emplace(8, chemist::AtomicBasisSet<double>("d-aug-cc-pvtz", 8, 0.0, 0.0, 0.0));

    basis_map.at(8).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{5.0800000000e-04,3.9290000000e-03,2.0243000000e-02,7.9181000000e-02,2.3068700000e-01,4.3311800000e-01,3.5026000000e-01,4.2728000000e-02,-8.1540000000e-03,2.3810000000e-03},
        std::vector<double>{1.5330000000e+04,2.2990000000e+03,5.2240000000e+02,1.4730000000e+02,4.7550000000e+01,1.6760000000e+01,6.2070000000e+00,1.7520000000e+00,6.8820000000e-01,2.3840000000e-01});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.7520000000e+00});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{-1.1500000000e-04,-8.9500000000e-04,-4.6360000000e-03,-1.8724000000e-02,-5.8463000000e-02,-1.3646300000e-01,-1.7574000000e-01,1.6093400000e-01,6.0341800000e-01,3.7876500000e-01},
        std::vector<double>{1.5330000000e+04,2.2990000000e+03,5.2240000000e+02,1.4730000000e+02,4.7550000000e+01,1.6760000000e+01,6.2070000000e+00,1.7520000000e+00,6.8820000000e-01,2.3840000000e-01});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.3840000000e-01});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{7.3760000000e-02});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.2800000000e-02});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{7.1560000000e-01});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.5928000000e-02,9.9740000000e-02,3.1049200000e-01,4.9102600000e-01,3.3633700000e-01},
        std::vector<double>{3.4460000000e+01,7.7490000000e+00,2.2800000000e+00,7.1560000000e-01,2.1400000000e-01});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.1400000000e-01});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{5.9740000000e-02});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.6700000000e-02});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.3140000000e+00});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{6.4500000000e-01});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.1400000000e-01});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{7.1000000000e-02});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.4280000000e+00});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{5.0000000000e-01});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.7500000000e-01});

    // Z = 9
    basis_map.emplace(9, chemist::AtomicBasisSet<double>("d-aug-cc-pvtz", 9, 0.0, 0.0, 0.0));

    basis_map.at(9).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{5.0700000000e-04,3.9230000000e-03,2.0200000000e-02,7.9010000000e-02,2.3043900000e-01,4.3287200000e-01,3.4996400000e-01,4.3233000000e-02,-7.8920000000e-03,2.3840000000e-03},
        std::vector<double>{1.9500000000e+04,2.9230000000e+03,6.6450000000e+02,1.8750000000e+02,6.0620000000e+01,2.1420000000e+01,7.9500000000e+00,2.2570000000e+00,8.8150000000e-01,3.0410000000e-01});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.2570000000e+00});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{-1.1700000000e-04,-9.1200000000e-04,-4.7170000000e-03,-1.9086000000e-02,-5.9655000000e-02,-1.4001000000e-01,-1.7678200000e-01,1.7162500000e-01,6.0504300000e-01,3.6951200000e-01},
        std::vector<double>{1.9500000000e+04,2.9230000000e+03,6.6450000000e+02,1.8750000000e+02,6.0620000000e+01,2.1420000000e+01,7.9500000000e+00,2.2570000000e+00,8.8150000000e-01,3.0410000000e-01});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.0410000000e-01});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{9.1580000000e-02});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.7600000000e-02});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{9.1320000000e-01});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.6665000000e-02,1.0447200000e-01,3.1726000000e-01,4.8734300000e-01,3.3460400000e-01},
        std::vector<double>{4.3880000000e+01,9.9260000000e+00,2.9300000000e+00,9.1320000000e-01,2.6720000000e-01});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.6720000000e-01});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{7.3610000000e-02});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.0300000000e-02});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.1070000000e+00});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{8.5500000000e-01});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.9200000000e-01});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{9.9700000000e-02});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.9170000000e+00});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{7.2400000000e-01});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.7300000000e-01});

    // Z = 10
    basis_map.emplace(10, chemist::AtomicBasisSet<double>("d-aug-cc-pvtz", 10, 0.0, 0.0, 0.0));

    basis_map.at(10).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{5.0200000000e-04,3.8810000000e-03,1.9997000000e-02,7.8418000000e-02,2.2967600000e-01,4.3272200000e-01,3.5064200000e-01,4.3911000000e-02,-7.6450000000e-03,2.3750000000e-03},
        std::vector<double>{2.4350000000e+04,3.6500000000e+03,8.2960000000e+02,2.3400000000e+02,7.5610000000e+01,2.6730000000e+01,9.9270000000e+00,2.8360000000e+00,1.1020000000e+00,3.7820000000e-01});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.8360000000e+00});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{-1.1800000000e-04,-9.1500000000e-04,-4.7370000000e-03,-1.9233000000e-02,-6.0369000000e-02,-1.4250800000e-01,-1.7771000000e-01,1.7735200000e-01,6.0583600000e-01,3.6510900000e-01},
        std::vector<double>{2.4350000000e+04,3.6500000000e+03,8.2960000000e+02,2.3400000000e+02,7.5610000000e+01,2.6730000000e+01,9.9270000000e+00,2.8360000000e+00,1.1020000000e+00,3.7820000000e-01});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.7820000000e-01});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.1330000000e-01});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.3900000000e-02});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.1430000000e+00});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.7151000000e-02,1.0765600000e-01,3.2168100000e-01,4.8523200000e-01,3.3258400000e-01},
        std::vector<double>{5.4700000000e+01,1.2430000000e+01,3.6790000000e+00,1.1430000000e+00,3.3000000000e-01});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.3000000000e-01});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{9.1750000000e-02});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.5500000000e-02});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.0140000000e+00});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.0960000000e+00});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.8600000000e-01});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.3600000000e-01});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.5440000000e+00});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.0840000000e+00});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.6200000000e-01});


    bsm.insert("d-aug-cc-pvtz", basis_map);
} // function load_d_dash_aug_dash_cc_dash_pvtz

} // namespace chemcache::basis_sets
