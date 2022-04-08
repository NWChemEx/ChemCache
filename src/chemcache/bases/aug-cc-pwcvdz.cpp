/*
 * This file is autogenerated by: generate_basis.py 
 * 
 * NOTE: Any modifications made in this file will be lost next time 
 *       generate_basis.py is run.
 */

#include "../basis_set_list.hpp"
# include <chemist/chemist.hpp>

namespace chemcache::basis_sets {

void load_aug_dash_cc_dash_pwcvdz(chemist::BasisSetManager& bsm) {
    chemist::BasisSetManager::ao_basis_map basis_map;

    // Z = 5
    basis_map.emplace(5, chemist::AtomicBasisSet<double>("aug-cc-pwcvdz", 5, 0.0, 0.0, 0.0));

    basis_map.at(5).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{6.9600000000e-04,5.3530000000e-03,2.7134000000e-02,1.0138000000e-01,2.7205500000e-01,4.4840300000e-01,2.9012300000e-01,1.4322000000e-02,-3.4860000000e-03},
        std::vector<double>{4.5700000000e+03,6.8590000000e+02,1.5650000000e+02,4.4470000000e+01,1.4480000000e+01,5.1310000000e+00,1.8980000000e+00,3.3290000000e-01,1.0430000000e-01});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{-1.3900000000e-04,-1.0970000000e-03,-5.4440000000e-03,-2.1916000000e-02,-5.9751000000e-02,-1.3873200000e-01,-1.3148200000e-01,5.3952600000e-01,5.8077400000e-01},
        std::vector<double>{4.5700000000e+03,6.8590000000e+02,1.5650000000e+02,4.4470000000e+01,1.4480000000e+01,5.1310000000e+00,1.8980000000e+00,3.3290000000e-01,1.0430000000e-01});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.0430000000e-01});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.4630000000e+00});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.1050000000e-02});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.5660000000e+00});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{3.5481000000e-02,1.9807200000e-01,5.0523000000e-01,4.7949900000e-01},
        std::vector<double>{6.0010000000e+00,1.2410000000e+00,3.3640000000e-01,9.5380000000e-02});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{9.5380000000e-02});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.3780000000e-02});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.4300000000e-01});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{9.0400000000e-02});

    // Z = 6
    basis_map.emplace(6, chemist::AtomicBasisSet<double>("aug-cc-pwcvdz", 6, 0.0, 0.0, 0.0));

    basis_map.at(6).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{6.9200000000e-04,5.3290000000e-03,2.7077000000e-02,1.0171800000e-01,2.7474000000e-01,4.4856400000e-01,2.8507400000e-01,1.5204000000e-02,-3.1910000000e-03},
        std::vector<double>{6.6650000000e+03,1.0000000000e+03,2.2800000000e+02,6.4710000000e+01,2.1060000000e+01,7.4950000000e+00,2.7970000000e+00,5.2150000000e-01,1.5960000000e-01});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{-1.4600000000e-04,-1.1540000000e-03,-5.7250000000e-03,-2.3312000000e-02,-6.3955000000e-02,-1.4998100000e-01,-1.2726200000e-01,5.4452900000e-01,5.8049600000e-01},
        std::vector<double>{6.6650000000e+03,1.0000000000e+03,2.2800000000e+02,6.4710000000e+01,2.1060000000e+01,7.4950000000e+00,2.7970000000e+00,5.2150000000e-01,1.5960000000e-01});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.5960000000e-01});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.3390000000e+00});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.6900000000e-02});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{7.2350000000e+00});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{3.8109000000e-02,2.0948000000e-01,5.0855700000e-01,4.6884200000e-01},
        std::vector<double>{9.4390000000e+00,2.0020000000e+00,5.4560000000e-01,1.5170000000e-01});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.5170000000e-01});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.0410000000e-02});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{5.5000000000e-01});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.5100000000e-01});

    // Z = 7
    basis_map.emplace(7, chemist::AtomicBasisSet<double>("aug-cc-pwcvdz", 7, 0.0, 0.0, 0.0));

    basis_map.at(7).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{7.0000000000e-04,5.3890000000e-03,2.7406000000e-02,1.0320700000e-01,2.7872300000e-01,4.4854000000e-01,2.7823800000e-01,1.5440000000e-02,-2.8640000000e-03},
        std::vector<double>{9.0460000000e+03,1.3570000000e+03,3.0930000000e+02,8.7730000000e+01,2.8560000000e+01,1.0210000000e+01,3.8380000000e+00,7.4660000000e-01,2.2480000000e-01});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{-1.5300000000e-04,-1.2080000000e-03,-5.9920000000e-03,-2.4544000000e-02,-6.7459000000e-02,-1.5807800000e-01,-1.2183100000e-01,5.4900300000e-01,5.7881500000e-01},
        std::vector<double>{9.0460000000e+03,1.3570000000e+03,3.0930000000e+02,8.7730000000e+01,2.8560000000e+01,1.0210000000e+01,3.8380000000e+00,7.4660000000e-01,2.2480000000e-01});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.2480000000e-01});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.3200000000e+00});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{6.1240000000e-02});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.0239000000e+01});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{3.9919000000e-02,2.1716900000e-01,5.1031900000e-01,4.6221400000e-01},
        std::vector<double>{1.3550000000e+01,2.9170000000e+00,7.9730000000e-01,2.1850000000e-01});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.1850000000e-01});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{5.6110000000e-02});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{8.1700000000e-01});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.3000000000e-01});

    // Z = 8
    basis_map.emplace(8, chemist::AtomicBasisSet<double>("aug-cc-pwcvdz", 8, 0.0, 0.0, 0.0));

    basis_map.at(8).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{7.1000000000e-04,5.4700000000e-03,2.7837000000e-02,1.0480000000e-01,2.8306200000e-01,4.4871900000e-01,2.7095200000e-01,1.5458000000e-02,-2.5850000000e-03},
        std::vector<double>{1.1720000000e+04,1.7590000000e+03,4.0080000000e+02,1.1370000000e+02,3.7030000000e+01,1.3270000000e+01,5.0250000000e+00,1.0130000000e+00,3.0230000000e-01});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{-1.6000000000e-04,-1.2630000000e-03,-6.2670000000e-03,-2.5716000000e-02,-7.0924000000e-02,-1.6541100000e-01,-1.1695500000e-01,5.5736800000e-01,5.7275900000e-01},
        std::vector<double>{1.1720000000e+04,1.7590000000e+03,4.0080000000e+02,1.1370000000e+02,3.7030000000e+01,1.3270000000e+01,5.0250000000e+00,1.0130000000e+00,3.0230000000e-01});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.0230000000e-01});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{5.4590000000e+00});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{7.8960000000e-02});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.3320000000e+01});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{4.3018000000e-02,2.2891300000e-01,5.0872800000e-01,4.6053100000e-01},
        std::vector<double>{1.7700000000e+01,3.8540000000e+00,1.0460000000e+00,2.7530000000e-01});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.7530000000e-01});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{6.8560000000e-02});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.1850000000e+00});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.3200000000e-01});

    // Z = 9
    basis_map.emplace(9, chemist::AtomicBasisSet<double>("aug-cc-pwcvdz", 9, 0.0, 0.0, 0.0));

    basis_map.at(9).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{7.2100000000e-04,5.5530000000e-03,2.8267000000e-02,1.0644400000e-01,2.8681400000e-01,4.4864100000e-01,2.6476100000e-01,1.5333000000e-02,-2.3320000000e-03},
        std::vector<double>{1.4710000000e+04,2.2070000000e+03,5.0280000000e+02,1.4260000000e+02,4.6470000000e+01,1.6700000000e+01,6.3560000000e+00,1.3160000000e+00,3.8970000000e-01});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{-1.6500000000e-04,-1.3080000000e-03,-6.4950000000e-03,-2.6691000000e-02,-7.3690000000e-02,-1.7077600000e-01,-1.1232700000e-01,5.6281400000e-01,5.6877800000e-01},
        std::vector<double>{1.4710000000e+04,2.2070000000e+03,5.0280000000e+02,1.4260000000e+02,4.6470000000e+01,1.6700000000e+01,6.3560000000e+00,1.3160000000e+00,3.8970000000e-01});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.8970000000e-01});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{6.7240000000e+00});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{9.8630000000e-02});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.6803000000e+01});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{4.4878000000e-02,2.3571800000e-01,5.0852100000e-01,4.5812000000e-01},
        std::vector<double>{2.2670000000e+01,4.9770000000e+00,1.3470000000e+00,3.4710000000e-01});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.4710000000e-01});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{8.5020000000e-02});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.6400000000e+00});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.6400000000e-01});

    // Z = 10
    basis_map.emplace(10, chemist::AtomicBasisSet<double>("aug-cc-pwcvdz", 10, 0.0, 0.0, 0.0));

    basis_map.at(10).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{7.3800000000e-04,5.6770000000e-03,2.8883000000e-02,1.0854000000e-01,2.9090700000e-01,4.4832400000e-01,2.5802600000e-01,1.5063000000e-02,-2.1000000000e-03},
        std::vector<double>{1.7880000000e+04,2.6830000000e+03,6.1150000000e+02,1.7350000000e+02,5.6640000000e+01,2.0420000000e+01,7.8100000000e+00,1.6530000000e+00,4.8690000000e-01});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{-1.7200000000e-04,-1.3570000000e-03,-6.7370000000e-03,-2.7663000000e-02,-7.6208000000e-02,-1.7522700000e-01,-1.0703800000e-01,5.6705000000e-01,5.6521600000e-01},
        std::vector<double>{1.7880000000e+04,2.6830000000e+03,6.1150000000e+02,1.7350000000e+02,5.6640000000e+01,2.0420000000e+01,7.8100000000e+00,1.6530000000e+00,4.8690000000e-01});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.8690000000e-01});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{8.1160000000e+00});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.2300000000e-01});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.0658000000e+01});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{4.6087000000e-02,2.4018100000e-01,5.0874400000e-01,4.5566000000e-01},
        std::vector<double>{2.8390000000e+01,6.2700000000e+00,1.6950000000e+00,4.3170000000e-01});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.3170000000e-01});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.0640000000e-01});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.2020000000e+00});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{6.3100000000e-01});

    // Z = 13
    basis_map.emplace(13, chemist::AtomicBasisSet<double>("aug-cc-pwcvdz", 13, 0.0, 0.0, 0.0));

    basis_map.at(13).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{2.9025000000e-04,2.2506400000e-03,1.1645900000e-02,4.6737700000e-02,1.4629900000e-01,3.3028300000e-01,4.1586100000e-01,1.8925300000e-01,1.1588900000e-02,-1.2838500000e-03,4.2588300000e-04,-1.9928000000e-04},
        std::vector<double>{6.4150000000e+04,9.6170000000e+03,2.1890000000e+03,6.2050000000e+02,2.0270000000e+02,7.3150000000e+01,2.8550000000e+01,1.1770000000e+01,3.3000000000e+00,1.1730000000e+00,1.7520000000e-01,6.4730000000e-02});
    basis_map.at(13).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{-7.5804800000e-05,-5.8179100000e-04,-3.0811300000e-03,-1.2311200000e-02,-4.1978100000e-02,-1.0337100000e-01,-1.9630800000e-01,-8.3000200000e-02,5.4104000000e-01,5.7879600000e-01,2.8814700000e-02,-9.5379500000e-03},
        std::vector<double>{6.4150000000e+04,9.6170000000e+03,2.1890000000e+03,6.2050000000e+02,2.0270000000e+02,7.3150000000e+01,2.8550000000e+01,1.1770000000e+01,3.3000000000e+00,1.1730000000e+00,1.7520000000e-01,6.4730000000e-02});
    basis_map.at(13).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.7507800000e-05,1.3420800000e-04,7.1244200000e-04,2.8433000000e-03,9.7684200000e-03,2.4185000000e-02,4.7499300000e-02,2.0362100000e-02,-1.5878800000e-01,-3.1169400000e-01,6.2014700000e-01,5.2094300000e-01},
        std::vector<double>{6.4150000000e+04,9.6170000000e+03,2.1890000000e+03,6.2050000000e+02,2.0270000000e+02,7.3150000000e+01,2.8550000000e+01,1.1770000000e+01,3.3000000000e+00,1.1730000000e+00,1.7520000000e-01,6.4730000000e-02});
    basis_map.at(13).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{6.4730000000e-02});
    basis_map.at(13).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.5990000000e+00});
    basis_map.at(13).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.3100000000e-02});
    basis_map.at(13).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{4.0684700000e-03,3.0681500000e-02,1.2914900000e-01,3.2083100000e-01,4.5381500000e-01,2.7506600000e-01,1.9080700000e-02,-3.1284800000e-03},
        std::vector<double>{2.5880000000e+02,6.0890000000e+01,1.9140000000e+01,6.8810000000e+00,2.5740000000e+00,9.5720000000e-01,2.0990000000e-01,5.9860000000e-02});
    basis_map.at(13).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{-7.4805300000e-04,-5.4579600000e-03,-2.4537100000e-02,-5.8213800000e-02,-9.8375600000e-02,-2.6006400000e-02,4.6402000000e-01,6.4887000000e-01},
        std::vector<double>{2.5880000000e+02,6.0890000000e+01,1.9140000000e+01,6.8810000000e+00,2.5740000000e+00,9.5720000000e-01,2.0990000000e-01,5.9860000000e-02});
    basis_map.at(13).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{5.9860000000e-02});
    basis_map.at(13).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{9.2800000000e-01});
    basis_map.at(13).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.5300000000e-02});
    basis_map.at(13).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.7920000000e+00});
    basis_map.at(13).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.8900000000e-01});
    basis_map.at(13).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{5.3500000000e-02});

    // Z = 14
    basis_map.emplace(14, chemist::AtomicBasisSet<double>("aug-cc-pwcvdz", 14, 0.0, 0.0, 0.0));

    basis_map.at(14).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{2.7044300000e-04,2.0971700000e-03,1.0850600000e-02,4.3675400000e-02,1.3765300000e-01,3.1664400000e-01,4.1858100000e-01,2.1021200000e-01,1.4495200000e-02,-2.0359000000e-03,6.2418600000e-04,-2.8287200000e-04},
        std::vector<double>{7.8860000000e+04,1.1820000000e+04,2.6920000000e+03,7.6340000000e+02,2.4960000000e+02,9.0280000000e+01,3.5290000000e+01,1.4510000000e+01,4.0530000000e+00,1.4820000000e+00,2.5170000000e-01,9.2430000000e-02});
    basis_map.at(14).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{-7.2317700000e-05,-5.5511600000e-04,-2.9380500000e-03,-1.1768700000e-02,-4.0290700000e-02,-1.0060900000e-01,-1.9652800000e-01,-1.0238200000e-01,5.2719000000e-01,5.9325100000e-01,3.3265200000e-02,-9.7366200000e-03},
        std::vector<double>{7.8860000000e+04,1.1820000000e+04,2.6920000000e+03,7.6340000000e+02,2.4960000000e+02,9.0280000000e+01,3.5290000000e+01,1.4510000000e+01,4.0530000000e+00,1.4820000000e+00,2.5170000000e-01,9.2430000000e-02});
    basis_map.at(14).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.8511300000e-05,1.4223600000e-04,7.5218500000e-04,3.0227900000e-03,1.0367700000e-02,2.6256300000e-02,5.2398900000e-02,2.9095900000e-02,-1.7800300000e-01,-3.4687400000e-01,6.2302000000e-01,5.3771200000e-01},
        std::vector<double>{7.8860000000e+04,1.1820000000e+04,2.6920000000e+03,7.6340000000e+02,2.4960000000e+02,9.0280000000e+01,3.5290000000e+01,1.4510000000e+01,4.0530000000e+00,1.4820000000e+00,2.5170000000e-01,9.2430000000e-02});
    basis_map.at(14).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{9.2430000000e-02});
    basis_map.at(14).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.4990000000e+00});
    basis_map.at(14).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.3200000000e-02});
    basis_map.at(14).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{3.9265600000e-03,2.9881100000e-02,1.2721200000e-01,3.2094300000e-01,4.5542900000e-01,2.6856300000e-01,1.8833600000e-02,-2.6243100000e-03},
        std::vector<double>{3.1590000000e+02,7.4420000000e+01,2.3480000000e+01,8.4880000000e+00,3.2170000000e+00,1.2290000000e+00,2.9640000000e-01,8.7680000000e-02});
    basis_map.at(14).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{-8.5830200000e-04,-6.3032800000e-03,-2.8825500000e-02,-6.9456000000e-02,-1.1949300000e-01,-1.9958100000e-02,5.1026800000e-01,6.0038200000e-01},
        std::vector<double>{3.1590000000e+02,7.4420000000e+01,2.3480000000e+01,8.4880000000e+00,3.2170000000e+00,1.2290000000e+00,2.9640000000e-01,8.7680000000e-02});
    basis_map.at(14).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{8.7680000000e-02});
    basis_map.at(14).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.3880000000e+00});
    basis_map.at(14).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.5000000000e-02});
    basis_map.at(14).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.2770000000e+00});
    basis_map.at(14).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.7500000000e-01});
    basis_map.at(14).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{8.2300000000e-02});

    // Z = 15
    basis_map.emplace(15, chemist::AtomicBasisSet<double>("aug-cc-pwcvdz", 15, 0.0, 0.0, 0.0));

    basis_map.at(15).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{2.5550900000e-04,1.9819300000e-03,1.0276000000e-02,4.1482300000e-02,1.3198400000e-01,3.0866200000e-01,4.2064700000e-01,2.2287800000e-01,1.6403500000e-02,-2.5425500000e-03,7.4805000000e-04,-3.3096300000e-04},
        std::vector<double>{9.4840000000e+04,1.4220000000e+04,3.2360000000e+03,9.1710000000e+02,2.9950000000e+02,1.0810000000e+02,4.2180000000e+01,1.7280000000e+01,4.8580000000e+00,1.8180000000e+00,3.3720000000e-01,1.2320000000e-01});
    basis_map.at(15).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{-6.9693900000e-05,-5.3526600000e-04,-2.8370900000e-03,-1.1398300000e-02,-3.9292900000e-02,-9.9636400000e-02,-1.9798300000e-01,-1.1486000000e-01,5.1859500000e-01,6.0184700000e-01,3.6861200000e-02,-9.7075900000e-03},
        std::vector<double>{9.4840000000e+04,1.4220000000e+04,3.2360000000e+03,9.1710000000e+02,2.9950000000e+02,1.0810000000e+02,4.2180000000e+01,1.7280000000e+01,4.8580000000e+00,1.8180000000e+00,3.3720000000e-01,1.2320000000e-01});
    basis_map.at(15).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.9119900000e-05,1.4722300000e-04,7.7791200000e-04,3.1454600000e-03,1.0820000000e-02,2.7995700000e-02,5.6397800000e-02,3.5819000000e-02,-1.9338700000e-01,-3.7209700000e-01,6.2424600000e-01,5.5172100000e-01},
        std::vector<double>{9.4840000000e+04,1.4220000000e+04,3.2360000000e+03,9.1710000000e+02,2.9950000000e+02,1.0810000000e+02,4.2180000000e+01,1.7280000000e+01,4.8580000000e+00,1.8180000000e+00,3.3720000000e-01,1.2320000000e-01});
    basis_map.at(15).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.2320000000e-01});
    basis_map.at(15).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{5.4740000000e+00});
    basis_map.at(15).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.1700000000e-02});
    basis_map.at(15).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{3.9500500000e-03,3.0249200000e-02,1.2955400000e-01,3.2759400000e-01,4.5699200000e-01,2.5308600000e-01,1.6879800000e-02,-2.0709300000e-03},
        std::vector<double>{3.7050000000e+02,8.7330000000e+01,2.7590000000e+01,1.0000000000e+01,3.8250000000e+00,1.4940000000e+00,3.9210000000e-01,1.1860000000e-01});
    basis_map.at(15).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{-9.5983200000e-04,-7.1117700000e-03,-3.2712200000e-02,-7.9578400000e-02,-1.3501600000e-01,-9.1058500000e-03,5.3780200000e-01,5.6906600000e-01},
        std::vector<double>{3.7050000000e+02,8.7330000000e+01,2.7590000000e+01,1.0000000000e+01,3.8250000000e+00,1.4940000000e+00,3.9210000000e-01,1.1860000000e-01});
    basis_map.at(15).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.1860000000e-01});
    basis_map.at(15).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.9390000000e+00});
    basis_map.at(15).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.4300000000e-02});
    basis_map.at(15).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.7600000000e+00});
    basis_map.at(15).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.7300000000e-01});
    basis_map.at(15).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.1300000000e-01});

    // Z = 16
    basis_map.emplace(16, chemist::AtomicBasisSet<double>("aug-cc-pwcvdz", 16, 0.0, 0.0, 0.0));

    basis_map.at(16).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{2.4763500000e-04,1.9202600000e-03,9.9619200000e-03,4.0297500000e-02,1.2860400000e-01,3.0348000000e-01,4.2143200000e-01,2.3078100000e-01,1.7897100000e-02,-2.9751600000e-03,8.4952200000e-04,-3.6793600000e-04},
        std::vector<double>{1.1080000000e+05,1.6610000000e+04,3.7810000000e+03,1.0710000000e+03,3.4980000000e+02,1.2630000000e+02,4.9260000000e+01,2.0160000000e+01,5.7200000000e+00,2.1820000000e+00,4.3270000000e-01,1.5700000000e-01});
    basis_map.at(16).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{-6.8703900000e-05,-5.2768100000e-04,-2.7967100000e-03,-1.1265100000e-02,-3.8883400000e-02,-9.9502500000e-02,-1.9974000000e-01,-1.2336000000e-01,5.1319400000e-01,6.0712000000e-01,3.9675300000e-02,-9.4686400000e-03},
        std::vector<double>{1.1080000000e+05,1.6610000000e+04,3.7810000000e+03,1.0710000000e+03,3.4980000000e+02,1.2630000000e+02,4.9260000000e+01,2.0160000000e+01,5.7200000000e+00,2.1820000000e+00,4.3270000000e-01,1.5700000000e-01});
    basis_map.at(16).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.9907700000e-05,1.5348300000e-04,8.0950300000e-04,3.2897400000e-03,1.1296700000e-02,2.9638500000e-02,5.9985100000e-02,4.1324800000e-02,-2.0747400000e-01,-3.9288900000e-01,6.3284000000e-01,5.5692400000e-01},
        std::vector<double>{1.1080000000e+05,1.6610000000e+04,3.7810000000e+03,1.0710000000e+03,3.4980000000e+02,1.2630000000e+02,4.9260000000e+01,2.0160000000e+01,5.7200000000e+00,2.1820000000e+00,4.3270000000e-01,1.5700000000e-01});
    basis_map.at(16).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.5700000000e-01});
    basis_map.at(16).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{6.5010000000e+00});
    basis_map.at(16).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{5.0700000000e-02});
    basis_map.at(16).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{4.4754100000e-03,3.4170800000e-02,1.4425000000e-01,3.5392800000e-01,4.5908500000e-01,2.0638300000e-01,1.0214100000e-02,-6.0312200000e-05},
        std::vector<double>{3.9970000000e+02,9.4190000000e+01,2.9750000000e+01,1.0770000000e+01,4.1190000000e+00,1.6250000000e+00,4.7260000000e-01,1.4070000000e-01});
    basis_map.at(16).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{-1.1625100000e-03,-8.6566400000e-03,-3.9088600000e-02,-9.3462500000e-02,-1.4799400000e-01,3.0190400000e-02,5.6157300000e-01,5.3477600000e-01},
        std::vector<double>{3.9970000000e+02,9.4190000000e+01,2.9750000000e+01,1.0770000000e+01,4.1190000000e+00,1.6250000000e+00,4.7260000000e-01,1.4070000000e-01});
    basis_map.at(16).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.4070000000e-01});
    basis_map.at(16).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.5030000000e+00});
    basis_map.at(16).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.9900000000e-02});
    basis_map.at(16).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.2110000000e+00});
    basis_map.at(16).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.7900000000e-01});
    basis_map.at(16).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.5200000000e-01});

    // Z = 17
    basis_map.emplace(17, chemist::AtomicBasisSet<double>("aug-cc-pwcvdz", 17, 0.0, 0.0, 0.0));

    basis_map.at(17).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{2.4115300000e-04,1.8709500000e-03,9.7082700000e-03,3.9315300000e-02,1.2593200000e-01,2.9934100000e-01,4.2188600000e-01,2.3720100000e-01,1.9153100000e-02,-3.3479200000e-03,9.2988300000e-04,-3.9637900000e-04},
        std::vector<double>{1.2790000000e+05,1.9170000000e+04,4.3630000000e+03,1.2360000000e+03,4.0360000000e+02,1.4570000000e+02,5.6810000000e+01,2.3230000000e+01,6.6440000000e+00,2.5750000000e+00,5.3710000000e-01,1.9380000000e-01});
    basis_map.at(17).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{-6.7892200000e-05,-5.2183600000e-04,-2.7651300000e-03,-1.1153700000e-02,-3.8591900000e-02,-9.9484800000e-02,-2.0139200000e-01,-1.3031300000e-01,5.0944300000e-01,6.1072500000e-01,4.2154900000e-02,-9.2342700000e-03},
        std::vector<double>{1.2790000000e+05,1.9170000000e+04,4.3630000000e+03,1.2360000000e+03,4.0360000000e+02,1.4570000000e+02,5.6810000000e+01,2.3230000000e+01,6.6440000000e+00,2.5750000000e+00,5.3710000000e-01,1.9380000000e-01});
    basis_map.at(17).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{2.0498600000e-05,1.5829800000e-04,8.3363900000e-04,3.3988000000e-03,1.1673800000e-02,3.0962200000e-02,6.2953300000e-02,4.6025700000e-02,-2.1931200000e-01,-4.0877300000e-01,6.3846500000e-01,5.6236200000e-01},
        std::vector<double>{1.2790000000e+05,1.9170000000e+04,4.3630000000e+03,1.2360000000e+03,4.0360000000e+02,1.4570000000e+02,5.6810000000e+01,2.3230000000e+01,6.6440000000e+00,2.5750000000e+00,5.3710000000e-01,1.9380000000e-01});
    basis_map.at(17).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.9380000000e-01});
    basis_map.at(17).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{7.5910000000e+00});
    basis_map.at(17).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{6.0800000000e-02});
    basis_map.at(17).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{5.2598200000e-03,3.9833200000e-02,1.6465500000e-01,3.8732200000e-01,4.5707200000e-01,1.5163600000e-01,1.8161500000e-03,1.8829600000e-03},
        std::vector<double>{4.1760000000e+02,9.8330000000e+01,3.1040000000e+01,1.1190000000e+01,4.2490000000e+00,1.6240000000e+00,5.3220000000e-01,1.6200000000e-01});
    basis_map.at(17).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{-1.4357000000e-03,-1.0779600000e-02,-4.7007500000e-02,-1.1103000000e-01,-1.5327500000e-01,8.9460900000e-02,5.7944400000e-01,4.8327200000e-01},
        std::vector<double>{4.1760000000e+02,9.8330000000e+01,3.1040000000e+01,1.1190000000e+01,4.2490000000e+00,1.6240000000e+00,5.3220000000e-01,1.6200000000e-01});
    basis_map.at(17).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.6200000000e-01});
    basis_map.at(17).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.0920000000e+00});
    basis_map.at(17).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.6600000000e-02});
    basis_map.at(17).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.7010000000e+00});
    basis_map.at(17).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{6.0000000000e-01});
    basis_map.at(17).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.9600000000e-01});

    // Z = 18
    basis_map.emplace(18, chemist::AtomicBasisSet<double>("aug-cc-pwcvdz", 18, 0.0, 0.0, 0.0));

    basis_map.at(18).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{2.3670000000e-04,1.8352300000e-03,9.5286000000e-03,3.8628300000e-02,1.2408100000e-01,2.9647100000e-01,4.2206800000e-01,2.4171100000e-01,2.0050900000e-02,-3.6100000000e-03,9.7560700000e-04,-4.1131600000e-04},
        std::vector<double>{1.4570000000e+05,2.1840000000e+04,4.9720000000e+03,1.4080000000e+03,4.5970000000e+02,1.6590000000e+02,6.4690000000e+01,2.6440000000e+01,7.6280000000e+00,2.9960000000e+00,6.5040000000e-01,2.3370000000e-01});
    basis_map.at(18).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{-6.7491000000e-05,-5.1852200000e-04,-2.7482500000e-03,-1.1100700000e-02,-3.8482000000e-02,-9.9759900000e-02,-2.0308800000e-01,-1.3560800000e-01,5.0719500000e-01,6.1289800000e-01,4.4296800000e-02,-8.9927800000e-03},
        std::vector<double>{1.4570000000e+05,2.1840000000e+04,4.9720000000e+03,1.4080000000e+03,4.5970000000e+02,1.6590000000e+02,6.4690000000e+01,2.6440000000e+01,7.6280000000e+00,2.9960000000e+00,6.5040000000e-01,2.3370000000e-01});
    basis_map.at(18).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{2.1045700000e-05,1.6256500000e-04,8.5546300000e-04,3.4974500000e-03,1.2015600000e-02,3.2136800000e-02,6.5527900000e-02,4.9937000000e-02,-2.2976900000e-01,-4.2100600000e-01,6.4233100000e-01,5.6754000000e-01},
        std::vector<double>{1.4570000000e+05,2.1840000000e+04,4.9720000000e+03,1.4080000000e+03,4.5970000000e+02,1.6590000000e+02,6.4690000000e+01,2.6440000000e+01,7.6280000000e+00,2.9960000000e+00,6.5040000000e-01,2.3370000000e-01});
    basis_map.at(18).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.3370000000e-01});
    basis_map.at(18).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{8.7420000000e+00});
    basis_map.at(18).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{7.0900000000e-02});
    basis_map.at(18).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{5.7055500000e-03,4.3046000000e-02,1.7659100000e-01,4.0686300000e-01,4.5254900000e-01,1.2280100000e-01,-4.4599600000e-03,2.0522500000e-03},
        std::vector<double>{4.5370000000e+02,1.0680000000e+02,3.3730000000e+01,1.2130000000e+01,4.5940000000e+00,1.6780000000e+00,5.9090000000e-01,1.8520000000e-01});
    basis_map.at(18).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{-1.6065500000e-03,-1.2171400000e-02,-5.2078900000e-02,-1.2373700000e-01,-1.5161900000e-01,1.4242500000e-01,5.8450100000e-01,4.3754000000e-01},
        std::vector<double>{4.5370000000e+02,1.0680000000e+02,3.3730000000e+01,1.2130000000e+01,4.5940000000e+00,1.6780000000e+00,5.9090000000e-01,1.8520000000e-01});
    basis_map.at(18).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.8520000000e-01});
    basis_map.at(18).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.7110000000e+00});
    basis_map.at(18).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{5.3300000000e-02});
    basis_map.at(18).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.2510000000e+00});
    basis_map.at(18).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{7.3800000000e-01});
    basis_map.at(18).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.4000000000e-01});


    bsm.insert("aug-cc-pwcvdz", basis_map);
} // function load_aug_dash_cc_dash_pwcvdz

} // namespace chemcache::basis_sets
