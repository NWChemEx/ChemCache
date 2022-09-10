/*
 * This file is autogenerated by: generate_basis.py
 *
 * NOTE: Any modifications made in this file will be lost next time
 *       generate_basis.py is run.
 */

#include "../basis_set_list.hpp"
#include <chemist/chemist.hpp>

namespace chemcache::basis_sets {

void load_six_dash_311g_oparen_2df_comma_2pd_cparen_(
  chemist::BasisSetManager& bsm) {
    chemist::BasisSetManager::ao_basis_map basis_map;

    // Z = 1
    basis_map.emplace(
      1, chemist::AtomicBasisSet<double>("6-311g(2df,2pd)", 1, 0.0, 0.0, 0.0));

    basis_map.at(1).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{2.5493800000e-02, 1.9037300000e-01, 8.5216100000e-01},
      std::vector<double>{3.3865000000e+01, 5.0947900000e+00,
                          1.1587900000e+00});
    basis_map.at(1).add_shell(chemist::ShellType::pure, 0,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{3.2584000000e-01});
    basis_map.at(1).add_shell(chemist::ShellType::pure, 0,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{1.0274100000e-01});
    basis_map.at(1).add_shell(chemist::ShellType::pure, 1,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{1.5000000000e+00});
    basis_map.at(1).add_shell(chemist::ShellType::pure, 1,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{3.7500000000e-01});
    basis_map.at(1).add_shell(chemist::ShellType::pure, 2,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{1.0000000000e+00});

    // Z = 2
    basis_map.emplace(
      2, chemist::AtomicBasisSet<double>("6-311g(2df,2pd)", 2, 0.0, 0.0, 0.0));

    basis_map.at(2).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{2.8745200000e-02, 2.0806100000e-01, 8.3763500000e-01},
      std::vector<double>{9.8124300000e+01, 1.4768900000e+01,
                          3.3188300000e+00});
    basis_map.at(2).add_shell(chemist::ShellType::pure, 0,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{8.7404700000e-01});
    basis_map.at(2).add_shell(chemist::ShellType::pure, 0,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{2.4456400000e-01});
    basis_map.at(2).add_shell(chemist::ShellType::pure, 1,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{1.5000000000e+00});
    basis_map.at(2).add_shell(chemist::ShellType::pure, 1,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{3.7500000000e-01});
    basis_map.at(2).add_shell(chemist::ShellType::pure, 2,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{2.0000000000e+00});

    // Z = 3
    basis_map.emplace(
      3, chemist::AtomicBasisSet<double>("6-311g(2df,2pd)", 3, 0.0, 0.0, 0.0));

    basis_map.at(3).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{2.2870400000e-03, 1.7635000000e-02, 8.7343400000e-02,
                          2.8097700000e-01, 6.5874100000e-01, 1.1871200000e-01},
      std::vector<double>{9.0046000000e+02, 1.3443300000e+02, 3.0436500000e+01,
                          8.6263900000e+00, 2.4833200000e+00,
                          3.0317900000e-01});
    basis_map.at(3).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{9.3329300000e-02, 9.4304500000e-01,
                          -2.7982700000e-03},
      std::vector<double>{4.8689000000e+00, 8.5692400000e-01,
                          2.4322700000e-01});
    basis_map.at(3).add_shell(
      chemist::ShellType::pure, 1,
      std::vector<double>{3.2766100000e-02, 1.5979200000e-01, 8.8566700000e-01},
      std::vector<double>{4.8689000000e+00, 8.5692400000e-01,
                          2.4322700000e-01});
    basis_map.at(3).add_shell(chemist::ShellType::pure, 0,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{6.3507000000e-02});
    basis_map.at(3).add_shell(chemist::ShellType::pure, 1,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{6.3507000000e-02});
    basis_map.at(3).add_shell(chemist::ShellType::pure, 0,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{2.4368300000e-02});
    basis_map.at(3).add_shell(chemist::ShellType::pure, 1,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{2.4368300000e-02});
    basis_map.at(3).add_shell(chemist::ShellType::pure, 2,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{4.0000000000e-01});
    basis_map.at(3).add_shell(chemist::ShellType::pure, 2,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{1.0000000000e-01});
    basis_map.at(3).add_shell(chemist::ShellType::pure, 3,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{1.5000000000e-01});

    // Z = 4
    basis_map.emplace(
      4, chemist::AtomicBasisSet<double>("6-311g(2df,2pd)", 4, 0.0, 0.0, 0.0));

    basis_map.at(4).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{2.2857400000e-03, 1.7593800000e-02, 8.6331500000e-02,
                          2.8183500000e-01, 6.4059400000e-01, 1.4446700000e-01},
      std::vector<double>{1.6828000000e+03, 2.5171500000e+02, 5.7411600000e+01,
                          1.6517100000e+01, 4.8536400000e+00,
                          6.2686300000e-01});
    basis_map.at(4).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{1.0862100000e-01, 9.2730100000e-01,
                          -2.9716900000e-03},
      std::vector<double>{8.3093800000e+00, 1.7407500000e+00,
                          4.8581600000e-01});
    basis_map.at(4).add_shell(
      chemist::ShellType::pure, 1,
      std::vector<double>{3.6134400000e-02, 2.1695800000e-01, 8.4183900000e-01},
      std::vector<double>{8.3093800000e+00, 1.7407500000e+00,
                          4.8581600000e-01});
    basis_map.at(4).add_shell(chemist::ShellType::pure, 0,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{1.6361300000e-01});
    basis_map.at(4).add_shell(chemist::ShellType::pure, 1,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{1.6361300000e-01});
    basis_map.at(4).add_shell(chemist::ShellType::pure, 0,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{5.6728500000e-02});
    basis_map.at(4).add_shell(chemist::ShellType::pure, 1,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{5.6728500000e-02});
    basis_map.at(4).add_shell(chemist::ShellType::pure, 2,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{5.1000000000e-01});
    basis_map.at(4).add_shell(chemist::ShellType::pure, 2,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{1.2750000000e-01});
    basis_map.at(4).add_shell(chemist::ShellType::pure, 3,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{2.6000000000e-01});

    // Z = 5
    basis_map.emplace(
      5, chemist::AtomicBasisSet<double>("6-311g(2df,2pd)", 5, 0.0, 0.0, 0.0));

    basis_map.at(5).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{2.1537500000e-03, 1.6582300000e-02, 8.2187000000e-02,
                          2.7661800000e-01, 6.2931600000e-01, 1.7377000000e-01},
      std::vector<double>{2.8588900000e+03, 4.2814000000e+02, 9.7528200000e+01,
                          2.7969300000e+01, 8.2157700000e+00,
                          1.1127800000e+00});
    basis_map.at(5).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{1.1744300000e-01, 9.1800200000e-01,
                          -2.6510500000e-03},
      std::vector<double>{1.3241500000e+01, 3.0016600000e+00,
                          9.1285600000e-01});
    basis_map.at(5).add_shell(
      chemist::ShellType::pure, 1,
      std::vector<double>{4.1810000000e-02, 2.3657500000e-01, 8.1621400000e-01},
      std::vector<double>{1.3241500000e+01, 3.0016600000e+00,
                          9.1285600000e-01});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 0,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{3.1545400000e-01});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 1,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{3.1545400000e-01});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 0,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{9.8856300000e-02});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 1,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{9.8856300000e-02});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 2,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{8.0200000000e-01});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 2,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{2.0050000000e-01});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 3,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{5.0000000000e-01});

    // Z = 6
    basis_map.emplace(
      6, chemist::AtomicBasisSet<double>("6-311g(2df,2pd)", 6, 0.0, 0.0, 0.0));

    basis_map.at(6).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{1.9666500000e-03, 1.5230600000e-02, 7.6126900000e-02,
                          2.6080100000e-01, 6.1646200000e-01, 2.2100600000e-01},
      std::vector<double>{4.5632400000e+03, 6.8202400000e+02, 1.5497300000e+02,
                          4.4455300000e+01, 1.3029000000e+01,
                          1.8277300000e+00});
    basis_map.at(6).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{1.1466000000e-01, 9.1999900000e-01,
                          -3.0306800000e-03},
      std::vector<double>{2.0964200000e+01, 4.8033100000e+00,
                          1.4593300000e+00});
    basis_map.at(6).add_shell(
      chemist::ShellType::pure, 1,
      std::vector<double>{4.0248700000e-02, 2.3759400000e-01, 8.1585400000e-01},
      std::vector<double>{2.0964200000e+01, 4.8033100000e+00,
                          1.4593300000e+00});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 0,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{4.8345600000e-01});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 1,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{4.8345600000e-01});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 0,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{1.4558500000e-01});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 1,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{1.4558500000e-01});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 2,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{1.2520000000e+00});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 2,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{3.1300000000e-01});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 3,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{8.0000000000e-01});

    // Z = 7
    basis_map.emplace(
      7, chemist::AtomicBasisSet<double>("6-311g(2df,2pd)", 7, 0.0, 0.0, 0.0));

    basis_map.at(7).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{1.9697900000e-03, 1.4961300000e-02, 7.3500600000e-02,
                          2.4893700000e-01, 6.0246000000e-01, 2.5620200000e-01},
      std::vector<double>{6.2934800000e+03, 9.4904400000e+02, 2.1877600000e+02,
                          6.3691600000e+01, 1.8828200000e+01,
                          2.7202300000e+00});
    basis_map.at(7).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{1.1190600000e-01, 9.2166600000e-01,
                          -2.5691900000e-03},
      std::vector<double>{3.0633100000e+01, 7.0261400000e+00,
                          2.1120500000e+00});
    basis_map.at(7).add_shell(
      chemist::ShellType::pure, 1,
      std::vector<double>{3.8311900000e-02, 2.3740300000e-01, 8.1759200000e-01},
      std::vector<double>{3.0633100000e+01, 7.0261400000e+00,
                          2.1120500000e+00});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 0,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{6.8400900000e-01});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 1,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{6.8400900000e-01});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 0,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{2.0087800000e-01});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 1,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{2.0087800000e-01});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 2,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{1.8260000000e+00});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 2,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{4.5650000000e-01});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 3,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{1.0000000000e+00});

    // Z = 8
    basis_map.emplace(
      8, chemist::AtomicBasisSet<double>("6-311g(2df,2pd)", 8, 0.0, 0.0, 0.0));

    basis_map.at(8).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{1.8951500000e-03, 1.4385900000e-02, 7.0732000000e-02,
                          2.4000100000e-01, 5.9479700000e-01, 2.8080200000e-01},
      std::vector<double>{8.5885000000e+03, 1.2972300000e+03, 2.9929600000e+02,
                          8.7377100000e+01, 2.5678900000e+01,
                          3.7400400000e+00});
    basis_map.at(8).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{1.1388900000e-01, 9.2081100000e-01,
                          -3.2744700000e-03},
      std::vector<double>{4.2117500000e+01, 9.6283700000e+00,
                          2.8533200000e+00});
    basis_map.at(8).add_shell(
      chemist::ShellType::pure, 1,
      std::vector<double>{3.6511400000e-02, 2.3715300000e-01, 8.1970200000e-01},
      std::vector<double>{4.2117500000e+01, 9.6283700000e+00,
                          2.8533200000e+00});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 0,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{9.0566100000e-01});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 1,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{9.0566100000e-01});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 0,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{2.5561100000e-01});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 1,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{2.5561100000e-01});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 2,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{2.5840000000e+00});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 2,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{6.4600000000e-01});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 3,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{1.4000000000e+00});

    // Z = 9
    basis_map.emplace(
      9, chemist::AtomicBasisSet<double>("6-311g(2df,2pd)", 9, 0.0, 0.0, 0.0));

    basis_map.at(9).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{1.8009300000e-03, 1.3741900000e-02, 6.8133400000e-02,
                          2.3332500000e-01, 5.8908600000e-01, 2.9950500000e-01},
      std::vector<double>{1.1427100000e+04, 1.7223500000e+03, 3.9574600000e+02,
                          1.1513900000e+02, 3.3602600000e+01,
                          4.9190100000e+00});
    basis_map.at(9).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{1.1453600000e-01, 9.2051200000e-01,
                          -3.3780400000e-03},
      std::vector<double>{5.5444100000e+01, 1.2632300000e+01,
                          3.7175600000e+00});
    basis_map.at(9).add_shell(
      chemist::ShellType::pure, 1,
      std::vector<double>{3.5460900000e-02, 2.3745100000e-01, 8.2045800000e-01},
      std::vector<double>{5.5444100000e+01, 1.2632300000e+01,
                          3.7175600000e+00});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 0,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{1.1654500000e+00});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 1,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{1.1654500000e+00});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 0,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{3.2189200000e-01});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 1,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{3.2189200000e-01});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 2,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{3.5000000000e+00});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 2,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{8.7500000000e-01});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 3,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{1.8500000000e+00});

    // Z = 10
    basis_map.emplace(10, chemist::AtomicBasisSet<double>("6-311g(2df,2pd)", 10,
                                                          0.0, 0.0, 0.0));

    basis_map.at(10).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{1.8327600000e-03, 1.3882700000e-02, 6.8068700000e-02,
                          2.3132800000e-01, 5.8589000000e-01, 3.0588300000e-01},
      std::vector<double>{1.3995700000e+04, 2.1171000000e+03, 4.9042500000e+02,
                          1.4383300000e+02, 4.1926500000e+01,
                          6.1568400000e+00});
    basis_map.at(10).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{1.1914900000e-01, 9.1737500000e-01,
                          -4.0583900000e-03},
      std::vector<double>{6.9121100000e+01, 1.5835000000e+01,
                          4.6732600000e+00});
    basis_map.at(10).add_shell(
      chemist::ShellType::pure, 1,
      std::vector<double>{3.5657400000e-02, 2.3947700000e-01, 8.1846100000e-01},
      std::vector<double>{6.9121100000e+01, 1.5835000000e+01,
                          4.6732600000e+00});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.4575600000e+00});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.4575600000e+00});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{3.9705700000e-01});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{3.9705700000e-01});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 2,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{4.6080000000e+00});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 2,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.1520000000e+00});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 3,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{2.5000000000e+00});

    // Z = 19
    basis_map.emplace(19, chemist::AtomicBasisSet<double>("6-311g(2df,2pd)", 19,
                                                          0.0, 0.0, 0.0));

    basis_map.at(19).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{1.7664000000e-03, 9.1949700000e-03, 3.7455100000e-02,
                          1.2204500000e-01, 2.9899000000e-01},
      std::vector<double>{2.7369000000e+04, 6.2291700000e+03, 1.7645800000e+03,
                          5.7705100000e+02, 2.1024900000e+02});
    basis_map.at(19).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{4.0514700000e-01, 2.9253200000e-01},
      std::vector<double>{8.2617800000e+01, 3.3233200000e+01});
    basis_map.at(19).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{8.1064900000e+00});
    basis_map.at(19).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{3.3340300000e+00});
    basis_map.at(19).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{8.4554400000e-01});
    basis_map.at(19).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{3.2821600000e-01});
    basis_map.at(19).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{3.6403500000e-02});
    basis_map.at(19).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.7646300000e-02});
    basis_map.at(19).add_shell(
      chemist::ShellType::pure, 1,
      std::vector<double>{2.1842900000e-03, 1.7589100000e-02, 8.1777500000e-02},
      std::vector<double>{8.9105400000e+02, 2.1101600000e+02,
                          6.7671400000e+01});
    basis_map.at(19).add_shell(
      chemist::ShellType::pure, 1,
      std::vector<double>{2.4565600000e-01, 4.3398400000e-01, 3.6237700000e-01},
      std::vector<double>{2.5271500000e+01, 1.0139000000e+01,
                          4.2018600000e+00});
    basis_map.at(19).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.6250700000e+00});
    basis_map.at(19).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{6.4377000000e-01});
    basis_map.at(19).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{2.4613000000e-01});
    basis_map.at(19).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{4.5440000000e-02});
    basis_map.at(19).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.6160000000e-02});
    basis_map.at(19).add_shell(
      chemist::ShellType::pure, 2,
      std::vector<double>{3.1601600000e-02, 1.5687900000e-01, 3.9058200000e-01},
      std::vector<double>{1.3370000000e+01, 3.4210000000e+00,
                          1.0630000000e+00});
    basis_map.at(19).add_shell(chemist::ShellType::pure, 2,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{4.5800000000e-01});
    basis_map.at(19).add_shell(chemist::ShellType::pure, 2,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.1450000000e-01});
    basis_map.at(19).add_shell(chemist::ShellType::pure, 3,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.1100000000e+00});

    // Z = 20
    basis_map.emplace(20, chemist::AtomicBasisSet<double>("6-311g(2df,2pd)", 20,
                                                          0.0, 0.0, 0.0));

    basis_map.at(20).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{1.7293200000e-03, 9.0022600000e-03, 3.6669900000e-02,
                          1.1941000000e-01, 2.9182500000e-01},
      std::vector<double>{3.0382500000e+04, 6.9150800000e+03, 1.9590200000e+03,
                          6.4093600000e+02, 2.3397700000e+02});
    basis_map.at(20).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{4.0441500000e-01, 2.9631300000e-01},
      std::vector<double>{9.2289200000e+01, 3.7254500000e+01});
    basis_map.at(20).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{9.1319800000e+00});
    basis_map.at(20).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{3.8177900000e+00});
    basis_map.at(20).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.0493500000e+00});
    basis_map.at(20).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{4.2866000000e-01});
    basis_map.at(20).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{6.2822600000e-02});
    basis_map.at(20).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{2.6016200000e-02});
    basis_map.at(20).add_shell(
      chemist::ShellType::pure, 1,
      std::vector<double>{2.0598600000e-03, 1.6650100000e-02, 7.7764600000e-02},
      std::vector<double>{1.0197600000e+03, 2.4159600000e+02,
                          7.7637000000e+01});
    basis_map.at(20).add_shell(
      chemist::ShellType::pure, 1,
      std::vector<double>{2.4180600000e-01, 4.3257800000e-01, 3.6732500000e-01},
      std::vector<double>{2.9115400000e+01, 1.1762600000e+01,
                          4.9228900000e+00});
    basis_map.at(20).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.9064500000e+00});
    basis_map.at(20).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{7.3690000000e-01});
    basis_map.at(20).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{2.7642000000e-01});
    basis_map.at(20).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{6.0270000000e-02});
    basis_map.at(20).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.7910000000e-02});
    basis_map.at(20).add_shell(
      chemist::ShellType::pure, 2,
      std::vector<double>{3.6894700000e-02, 1.7782000000e-01, 4.2551300000e-01},
      std::vector<double>{1.5080000000e+01, 3.9260000000e+00,
                          1.2330000000e+00});
    basis_map.at(20).add_shell(chemist::ShellType::pure, 2,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{5.2000000000e-01});
    basis_map.at(20).add_shell(chemist::ShellType::pure, 2,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.3000000000e-01});
    basis_map.at(20).add_shell(chemist::ShellType::pure, 3,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.3300000000e+00});

    bsm.insert("6-311g(2df,2pd)", basis_map);
} // function load_six_dash_311g_oparen_2df_comma_2pd_cparen_

} // namespace chemcache::basis_sets
