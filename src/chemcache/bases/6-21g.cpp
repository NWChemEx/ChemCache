/*
 * This file is autogenerated by: generate_basis.py
 *
 * NOTE: Any modifications made in this file will be lost next time
 *       generate_basis.py is run.
 */

#include "../basis_set_list.hpp"
#include <chemist/chemist.hpp>

namespace chemcache::basis_sets {

void load_six_dash_21g(chemist::BasisSetManager& bsm) {
    chemist::BasisSetManager::ao_basis_map basis_map;

    // Z = 1
    basis_map.emplace(
      1, chemist::AtomicBasisSet<double>("6-21g", 1, 0.0, 0.0, 0.0));

    basis_map.at(1).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{1.5628497870e-01, 9.0469087670e-01},
      std::vector<double>{5.4471780000e+00, 8.2454724000e-01});
    basis_map.at(1).add_shell(chemist::ShellType::pure, 0,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{1.8319158000e-01});

    // Z = 2
    basis_map.emplace(
      2, chemist::AtomicBasisSet<double>("6-21g", 2, 0.0, 0.0, 0.0));

    basis_map.at(2).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{1.7522987180e-01, 8.9348234650e-01},
      std::vector<double>{1.3626700000e+01, 1.9993500000e+00});
    basis_map.at(2).add_shell(chemist::ShellType::pure, 0,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{3.8299300000e-01});

    // Z = 3
    basis_map.emplace(
      3, chemist::AtomicBasisSet<double>("6-21g", 3, 0.0, 0.0, 0.0));

    basis_map.at(3).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{2.1509588790e-03, 1.6267691520e-02, 7.7638259540e-02,
                          2.4649487150e-01, 4.6750575630e-01, 3.4691481920e-01},
      std::vector<double>{6.4241800000e+02, 9.6516400000e+01, 2.2017400000e+01,
                          6.1764500000e+00, 1.9351100000e+00,
                          6.3957700000e-01});
    basis_map.at(3).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{-2.6312640580e-01, 1.1433874180e+00},
      std::vector<double>{5.4020500000e-01, 1.0225500000e-01});
    basis_map.at(3).add_shell(
      chemist::ShellType::pure, 1,
      std::vector<double>{1.6154597080e-01, 9.1566283470e-01},
      std::vector<double>{5.4020500000e-01, 1.0225500000e-01});
    basis_map.at(3).add_shell(chemist::ShellType::pure, 0,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{2.8564500000e-02});
    basis_map.at(3).add_shell(chemist::ShellType::pure, 1,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{2.8564500000e-02});

    // Z = 4
    basis_map.emplace(
      4, chemist::AtomicBasisSet<double>("6-21g", 4, 0.0, 0.0, 0.0));

    basis_map.at(4).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{1.9433606270e-03, 1.4825104780e-02, 7.2066223250e-02,
                          2.3702207650e-01, 4.6878915120e-01, 3.5638211500e-01},
      std::vector<double>{1.2645000000e+03, 1.8993000000e+02, 4.3127500000e+01,
                          1.2088900000e+01, 3.8079000000e+00,
                          1.2826600000e+00});
    basis_map.at(4).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{-4.2106406590e-01, 1.2240701920e+00},
      std::vector<double>{1.2954800000e+00, 2.6888100000e-01});
    basis_map.at(4).add_shell(
      chemist::ShellType::pure, 1,
      std::vector<double>{2.0513192370e-01, 8.8252767190e-01},
      std::vector<double>{1.2954800000e+00, 2.6888100000e-01});
    basis_map.at(4).add_shell(chemist::ShellType::pure, 0,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{7.7350100000e-02});
    basis_map.at(4).add_shell(chemist::ShellType::pure, 1,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{7.7350100000e-02});

    // Z = 5
    basis_map.emplace(
      5, chemist::AtomicBasisSet<double>("6-21g", 5, 0.0, 0.0, 0.0));

    basis_map.at(5).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{1.8498587790e-03, 1.4127690670e-02, 6.9269654270e-02,
                          2.3239284660e-01, 4.7015368960e-01, 3.6028776210e-01},
      std::vector<double>{2.0821200000e+03, 3.1231000000e+02, 7.0887400000e+01,
                          1.9852500000e+01, 6.2916100000e+00,
                          2.1286200000e+00});
    basis_map.at(5).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{-3.6866347730e-01, 1.1994448060e+00},
      std::vector<double>{2.2818700000e+00, 4.6524800000e-01});
    basis_map.at(5).add_shell(
      chemist::ShellType::pure, 1,
      std::vector<double>{2.3115190230e-01, 8.6676363370e-01},
      std::vector<double>{2.2818700000e+00, 4.6524800000e-01});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 0,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{1.2432800000e-01});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 1,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{1.2432800000e-01});

    // Z = 6
    basis_map.emplace(
      6, chemist::AtomicBasisSet<double>("6-21g", 6, 0.0, 0.0, 0.0));

    basis_map.at(6).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{1.8258801230e-03, 1.4056600940e-02, 6.8757004620e-02,
                          2.3042201550e-01, 4.6846303150e-01, 3.6278002440e-01},
      std::vector<double>{3.0475200000e+03, 4.5642400000e+02, 1.0365300000e+02,
                          2.9225800000e+01, 9.3486300000e+00,
                          3.1890400000e+00});
    basis_map.at(6).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{-3.9589516210e-01, 1.2158343560e+00},
      std::vector<double>{3.6649800000e+00, 7.7054500000e-01});
    basis_map.at(6).add_shell(
      chemist::ShellType::pure, 1,
      std::vector<double>{2.3645994660e-01, 8.6061880570e-01},
      std::vector<double>{3.6649800000e+00, 7.7054500000e-01});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 0,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{1.9585700000e-01});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 1,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{1.9585700000e-01});

    // Z = 7
    basis_map.emplace(
      7, chemist::AtomicBasisSet<double>("6-21g", 7, 0.0, 0.0, 0.0));

    basis_map.at(7).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{1.8454099110e-03, 1.4164499320e-02, 6.8632496710e-02,
                          2.2857398900e-01, 4.6616197760e-01, 3.6567198250e-01},
      std::vector<double>{4.1501100000e+03, 6.2008400000e+02, 1.4168800000e+02,
                          4.0336700000e+01, 1.3026700000e+01,
                          4.4700300000e+00});
    basis_map.at(7).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{-4.1330007740e-01, 1.2244172670e+00},
      std::vector<double>{5.4252200000e+00, 1.1491500000e+00});
    basis_map.at(7).add_shell(
      chemist::ShellType::pure, 1,
      std::vector<double>{2.3797201620e-01, 8.5895305860e-01},
      std::vector<double>{5.4252200000e+00, 1.1491500000e+00});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 0,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{2.8320500000e-01});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 1,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{2.8320500000e-01});

    // Z = 8
    basis_map.emplace(
      8, chemist::AtomicBasisSet<double>("6-21g", 8, 0.0, 0.0, 0.0));

    basis_map.at(8).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{1.8321688100e-03, 1.4104690840e-02, 6.8626155420e-02,
                          2.2937585100e-01, 4.6639869700e-01, 3.6417276340e-01},
      std::vector<double>{5.4722700000e+03, 8.1780600000e+02, 1.8644600000e+02,
                          5.3023000000e+01, 1.7180000000e+01,
                          5.9119600000e+00});
    basis_map.at(8).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{-4.0445358320e-01, 1.2215617610e+00},
      std::vector<double>{7.4029400000e+00, 1.5762000000e+00});
    basis_map.at(8).add_shell(
      chemist::ShellType::pure, 1,
      std::vector<double>{2.4458610700e-01, 8.5395537350e-01},
      std::vector<double>{7.4029400000e+00, 1.5762000000e+00});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 0,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{3.7368400000e-01});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 1,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{3.7368400000e-01});

    // Z = 9
    basis_map.emplace(
      9, chemist::AtomicBasisSet<double>("6-21g", 9, 0.0, 0.0, 0.0));

    basis_map.at(9).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{1.8846298750e-03, 1.3812099090e-02, 6.6249295620e-02,
                          2.2187498530e-01, 4.6084196950e-01, 3.7845297500e-01},
      std::vector<double>{6.7831900000e+03, 1.0424400000e+03, 2.4239800000e+02,
                          6.9632000000e+01, 2.2689400000e+01,
                          7.7963600000e+00});
    basis_map.at(9).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{-4.0732627770e-01, 1.2231378310e+00},
      std::vector<double>{9.7775900000e+00, 2.0861700000e+00});
    basis_map.at(9).add_shell(
      chemist::ShellType::pure, 1,
      std::vector<double>{2.4668000320e-01, 8.5232101100e-01},
      std::vector<double>{9.7775900000e+00, 2.0861700000e+00});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 0,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{4.8238300000e-01});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 1,
                              std::vector<double>{1.0000000000e+00},
                              std::vector<double>{4.8238300000e-01});

    // Z = 10
    basis_map.emplace(
      10, chemist::AtomicBasisSet<double>("6-21g", 10, 0.0, 0.0, 0.0));

    basis_map.at(10).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{1.9376592770e-03, 1.4806994480e-02, 7.2705472880e-02,
                          2.5262890580e-01, 4.9324181600e-01, 3.1316888320e-01},
      std::vector<double>{9.9932000000e+03, 1.4998900000e+03, 3.4195100000e+02,
                          9.4679600000e+01, 2.9734500000e+01,
                          1.0006300000e+01});
    basis_map.at(10).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{-4.0992232080e-01, 1.2243109580e+00},
      std::vector<double>{1.2483000000e+01, 2.6645100000e+00});
    basis_map.at(10).add_shell(
      chemist::ShellType::pure, 1,
      std::vector<double>{2.4745998360e-01, 8.5174294350e-01},
      std::vector<double>{1.2483000000e+01, 2.6645100000e+00});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{6.0625000000e-01});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{6.0625000000e-01});

    // Z = 11
    basis_map.emplace(
      11, chemist::AtomicBasisSet<double>("6-21g", 11, 0.0, 0.0, 0.0));

    basis_map.at(11).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{1.9376592770e-03, 1.4806994480e-02, 7.2705472880e-02,
                          2.5262890580e-01, 4.9324181600e-01, 3.1316888320e-01},
      std::vector<double>{9.9932000000e+03, 1.4998900000e+03, 3.4195100000e+02,
                          9.4679600000e+01, 2.9734500000e+01,
                          1.0006300000e+01});
    basis_map.at(11).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{-3.5420835040e-03, -4.3958843480e-02,
                          -1.0975210860e-01, 1.8739818540e-01, 6.4669963970e-01,
                          3.0605830270e-01},
      std::vector<double>{1.5096300000e+02, 3.5587800000e+01, 1.1168300000e+01,
                          3.9020100000e+00, 1.3817700000e+00,
                          4.6638200000e-01});
    basis_map.at(11).add_shell(
      chemist::ShellType::pure, 1,
      std::vector<double>{5.0016597100e-03, 3.5510897940e-02, 1.4282499170e-01,
                          3.3861998030e-01, 4.5157897380e-01, 2.7327098410e-01},
      std::vector<double>{1.5096300000e+02, 3.5587800000e+01, 1.1168300000e+01,
                          3.9020100000e+00, 1.3817700000e+00,
                          4.6638200000e-01});
    basis_map.at(11).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{-2.1966049750e-01, 1.0891224670e+00},
      std::vector<double>{5.0182400000e-01, 6.0945800000e-02});
    basis_map.at(11).add_shell(
      chemist::ShellType::pure, 1,
      std::vector<double>{9.0664879580e-03, 9.9720177540e-01},
      std::vector<double>{5.0182400000e-01, 6.0945800000e-02});
    basis_map.at(11).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{2.4434900000e-02});
    basis_map.at(11).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{2.4434900000e-02});

    // Z = 12
    basis_map.emplace(
      12, chemist::AtomicBasisSet<double>("6-21g", 12, 0.0, 0.0, 0.0));

    basis_map.at(12).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{1.9778293170e-03, 1.5113994780e-02, 7.3910774480e-02,
                          2.4919091400e-01, 4.8792783160e-01, 3.1966188960e-01},
      std::vector<double>{1.1722800000e+04, 1.7599300000e+03, 4.0084600000e+02,
                          1.1280700000e+02, 3.5999700000e+01,
                          1.2182800000e+01});
    basis_map.at(12).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{-3.2371704710e-03, -4.1007905970e-02,
                          -1.1260001640e-01, 1.4863302160e-01, 6.1649708980e-01,
                          3.6482905310e-01},
      std::vector<double>{1.8918000000e+02, 4.5211900000e+01, 1.4356300000e+01,
                          5.1388600000e+00, 1.9065200000e+00,
                          7.0588700000e-01});
    basis_map.at(12).add_shell(
      chemist::ShellType::pure, 1,
      std::vector<double>{4.9281299210e-03, 3.4988799440e-02, 1.4072499770e-01,
                          3.3364199470e-01, 4.4493999290e-01, 2.6925399570e-01},
      std::vector<double>{1.8918000000e+02, 4.5211900000e+01, 1.4356300000e+01,
                          5.1388600000e+00, 1.9065200000e+00,
                          7.0588700000e-01});
    basis_map.at(12).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{-3.6110259330e-01, 1.2150553610e+00},
      std::vector<double>{6.1134900000e-01, 1.4184100000e-01});
    basis_map.at(12).add_shell(
      chemist::ShellType::pure, 1,
      std::vector<double>{2.4263309200e-02, 9.8667337390e-01},
      std::vector<double>{6.1134900000e-01, 1.4184100000e-01});
    basis_map.at(12).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{4.6401100000e-02});
    basis_map.at(12).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{4.6401100000e-02});

    // Z = 13
    basis_map.emplace(
      13, chemist::AtomicBasisSet<double>("6-21g", 13, 0.0, 0.0, 0.0));

    basis_map.at(13).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{1.9426699470e-03, 1.4859899590e-02, 7.2849398000e-02,
                          2.4682999320e-01, 4.8725798660e-01, 3.2349599110e-01},
      std::vector<double>{1.3983100000e+04, 2.0987500000e+03, 4.7770500000e+02,
                          1.3436000000e+02, 4.2870900000e+01,
                          1.4518900000e+01});
    basis_map.at(13).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{-2.9261900280e-03, -3.7408300360e-02,
                          -1.1448700110e-01, 1.1563500110e-01, 6.1259500580e-01,
                          3.9379900370e-01},
      std::vector<double>{2.3966800000e+02, 5.7441900000e+01, 1.8285900000e+01,
                          6.5991400000e+00, 2.4904900000e+00,
                          9.4454500000e-01});
    basis_map.at(13).add_shell(
      chemist::ShellType::pure, 1,
      std::vector<double>{4.6028455820e-03, 3.3198968130e-02, 1.3628186920e-01,
                          3.3047568280e-01, 4.4914556890e-01, 2.6570374500e-01},
      std::vector<double>{2.3966800000e+02, 5.7441900000e+01, 1.8285900000e+01,
                          6.5991400000e+00, 2.4904900000e+00,
                          9.4454500000e-01});
    basis_map.at(13).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{-3.2032690910e-01, 1.1841196640e+00},
      std::vector<double>{9.4616000000e-01, 2.0250600000e-01});
    basis_map.at(13).add_shell(
      chemist::ShellType::pure, 1,
      std::vector<double>{5.1938280870e-02, 9.7265964170e-01},
      std::vector<double>{9.4616000000e-01, 2.0250600000e-01});
    basis_map.at(13).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{6.3908800000e-02});
    basis_map.at(13).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{6.3908800000e-02});

    // Z = 14
    basis_map.emplace(
      14, chemist::AtomicBasisSet<double>("6-21g", 14, 0.0, 0.0, 0.0));

    basis_map.at(14).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{1.9594802160e-03, 1.4928801640e-02, 7.2847808010e-02,
                          2.4613002710e-01, 4.8591405350e-01, 3.2500203580e-01},
      std::vector<double>{1.6115900000e+04, 2.4255800000e+03, 5.5386700000e+02,
                          1.5634000000e+02, 5.0068300000e+01,
                          1.7017800000e+01});
    basis_map.at(14).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{-2.7809414150e-03, -3.5714618170e-02,
                          -1.1498505850e-01, 9.3563447600e-02, 6.0301730680e-01,
                          4.1895921310e-01},
      std::vector<double>{2.9271800000e+02, 6.9873100000e+01, 2.2336300000e+01,
                          8.1503900000e+00, 3.1345800000e+00,
                          1.2254300000e+00});
    basis_map.at(14).add_shell(
      chemist::ShellType::pure, 1,
      std::vector<double>{4.4382645210e-03, 3.2667933280e-02, 1.3472113720e-01,
                          3.2867833480e-01, 4.4964045800e-01, 2.6137226620e-01},
      std::vector<double>{2.9271800000e+02, 6.9873100000e+01, 2.2336300000e+01,
                          8.1503900000e+00, 3.1345800000e+00,
                          1.2254300000e+00});
    basis_map.at(14).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{-3.7610787950e-01, 1.2516495990e+00},
      std::vector<double>{1.0791300000e+00, 3.0242200000e-01});
    basis_map.at(14).add_shell(
      chemist::ShellType::pure, 1,
      std::vector<double>{6.7102991120e-02, 9.5688287340e-01},
      std::vector<double>{1.0791300000e+00, 3.0242200000e-01});
    basis_map.at(14).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{9.3339200000e-02});
    basis_map.at(14).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{9.3339200000e-02});

    // Z = 15
    basis_map.emplace(
      15, chemist::AtomicBasisSet<double>("6-21g", 15, 0.0, 0.0, 0.0));

    basis_map.at(15).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{1.8515989230e-03, 1.4206191740e-02, 6.9999459280e-02,
                          2.4007886030e-01, 4.8476171800e-01, 3.3519980500e-01},
      std::vector<double>{1.9413300000e+04, 2.9094200000e+03, 6.6136400000e+02,
                          1.8575900000e+02, 5.9194300000e+01,
                          2.0031000000e+01});
    basis_map.at(15).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{-2.7821701050e-03, -3.6049901350e-02,
                          -1.1663100440e-01, 9.6832803640e-02, 6.1441802310e-01,
                          4.0379801520e-01},
      std::vector<double>{3.3947800000e+02, 8.1010100000e+01, 2.5878000000e+01,
                          9.4522100000e+00, 3.6656600000e+00,
                          1.4674600000e+00});
    basis_map.at(15).add_shell(
      chemist::ShellType::pure, 1,
      std::vector<double>{4.5646161910e-03, 3.3693571880e-02, 1.3975488340e-01,
                          3.3936171680e-01, 4.5092062370e-01, 2.3858580090e-01},
      std::vector<double>{3.3947800000e+02, 8.1010100000e+01, 2.5878000000e+01,
                          9.4522100000e+00, 3.6656600000e+00,
                          1.4674600000e+00});
    basis_map.at(15).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{-3.7149602190e-01, 1.2709934960e+00},
      std::vector<double>{1.2186500000e+00, 3.9554600000e-01});
    basis_map.at(15).add_shell(
      chemist::ShellType::pure, 1,
      std::vector<double>{9.1582310220e-02, 9.3492410430e-01},
      std::vector<double>{1.2186500000e+00, 3.9554600000e-01});
    basis_map.at(15).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.2281100000e-01});
    basis_map.at(15).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.2281100000e-01});

    // Z = 16
    basis_map.emplace(
      16, chemist::AtomicBasisSet<double>("6-21g", 16, 0.0, 0.0, 0.0));

    basis_map.at(16).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{1.8692408490e-03, 1.4230306460e-02, 6.9696231660e-02,
                          2.3848710830e-01, 4.8330721950e-01, 3.3807415360e-01},
      std::vector<double>{2.1917100000e+04, 3.3014900000e+03, 7.5414600000e+02,
                          2.1271100000e+02, 6.7989600000e+01,
                          2.3051500000e+01});
    basis_map.at(16).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{-2.3767704990e-03, -3.1693006650e-02,
                          -1.1331702380e-01, 5.6090011770e-02, 5.9225512430e-01,
                          4.5500609550e-01},
      std::vector<double>{4.2373500000e+02, 1.0071000000e+02, 3.2159900000e+01,
                          1.1807900000e+01, 4.6311000000e+00,
                          1.8702500000e+00});
    basis_map.at(16).add_shell(
      chemist::ShellType::pure, 1,
      std::vector<double>{4.0610099820e-03, 3.0681299860e-02, 1.3045199940e-01,
                          3.2720499850e-01, 4.5285099800e-01, 2.5604199890e-01},
      std::vector<double>{4.2373500000e+02, 1.0071000000e+02, 3.2159900000e+01,
                          1.1807900000e+01, 4.6311000000e+00,
                          1.8702500000e+00});
    basis_map.at(16).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{-2.8608885370e-01, 1.2280593720e+00},
      std::vector<double>{1.2238400000e+00, 4.5730300000e-01});
    basis_map.at(16).add_shell(
      chemist::ShellType::pure, 1,
      std::vector<double>{1.6477699470e-01, 8.7085497220e-01},
      std::vector<double>{1.2238400000e+00, 4.5730300000e-01});
    basis_map.at(16).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.4226900000e-01});
    basis_map.at(16).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.4226900000e-01});

    // Z = 17
    basis_map.emplace(
      17, chemist::AtomicBasisSet<double>("6-21g", 17, 0.0, 0.0, 0.0));

    basis_map.at(17).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{1.8329598480e-03, 1.4034198830e-02, 6.9097394260e-02,
                          2.3745198030e-01, 4.8303395990e-01, 3.3985597180e-01},
      std::vector<double>{2.5180100000e+04, 3.7803500000e+03, 8.6047400000e+02,
                          2.4214500000e+02, 7.7334900000e+01,
                          2.6247000000e+01});
    basis_map.at(17).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{-2.2973914170e-03, -3.0713718940e-02,
                          -1.1252806940e-01, 4.5016327760e-02, 5.8935336340e-01,
                          4.6520628680e-01},
      std::vector<double>{4.9176500000e+02, 1.1698400000e+02, 3.7415300000e+01,
                          1.3783400000e+01, 5.4521500000e+00,
                          2.2258800000e+00});
    basis_map.at(17).add_shell(
      chemist::ShellType::pure, 1,
      std::vector<double>{3.9894008790e-03, 3.0317706680e-02, 1.2988002860e-01,
                          3.2795107230e-01, 4.5352710000e-01, 2.5215405560e-01},
      std::vector<double>{4.9176500000e+02, 1.1698400000e+02, 3.7415300000e+01,
                          1.3783400000e+01, 5.4521500000e+00,
                          2.2258800000e+00});
    basis_map.at(17).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{-2.2240148410e-01, 1.1825225740e+00},
      std::vector<double>{1.3529900000e+00, 5.2695500000e-01});
    basis_map.at(17).add_shell(
      chemist::ShellType::pure, 1,
      std::vector<double>{2.1921579720e-01, 8.2232023930e-01},
      std::vector<double>{1.3529900000e+00, 5.2695500000e-01});
    basis_map.at(17).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.6671400000e-01});
    basis_map.at(17).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.6671400000e-01});

    // Z = 18
    basis_map.emplace(
      18, chemist::AtomicBasisSet<double>("6-21g", 18, 0.0, 0.0, 0.0));

    basis_map.at(18).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{1.8252601920e-03, 1.3968601470e-02, 6.8707307230e-02,
                          2.3620402490e-01, 4.8221405080e-01, 3.4204303600e-01},
      std::vector<double>{2.8348300000e+04, 4.2576200000e+03, 9.6985700000e+02,
                          2.7326300000e+02, 8.7369500000e+01,
                          2.9686700000e+01});
    basis_map.at(18).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{-2.1597208950e-03, -2.9077512060e-02,
                          -1.1082704600e-01, 2.7699911480e-02, 5.7761323950e-01,
                          4.8868820260e-01},
      std::vector<double>{5.7589100000e+02, 1.3681600000e+02, 4.3809800000e+01,
                          1.6209400000e+01, 6.4608400000e+00,
                          2.6511400000e+00});
    basis_map.at(18).add_shell(
      chemist::ShellType::pure, 1,
      std::vector<double>{3.8066498420e-03, 2.9230498790e-02, 1.2646699480e-01,
                          3.2350998660e-01, 4.5489598110e-01, 2.5662998940e-01},
      std::vector<double>{5.7589100000e+02, 1.3681600000e+02, 4.3809800000e+01,
                          1.6209400000e+01, 6.4608400000e+00,
                          2.6511400000e+00});
    basis_map.at(18).add_shell(
      chemist::ShellType::pure, 0,
      std::vector<double>{-1.7686556850e-01, 1.1468972020e+00},
      std::vector<double>{1.5420900000e+00, 6.0726700000e-01});
    basis_map.at(18).add_shell(
      chemist::ShellType::pure, 1,
      std::vector<double>{2.5568701300e-01, 7.8984204010e-01},
      std::vector<double>{1.5420900000e+00, 6.0726700000e-01});
    basis_map.at(18).add_shell(chemist::ShellType::pure, 0,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.9537300000e-01});
    basis_map.at(18).add_shell(chemist::ShellType::pure, 1,
                               std::vector<double>{1.0000000000e+00},
                               std::vector<double>{1.9537300000e-01});

    bsm.insert("6-21g", basis_map);
} // function load_six_dash_21g

} // namespace chemcache::basis_sets
