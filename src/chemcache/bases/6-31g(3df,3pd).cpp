/*
 * This file is autogenerated by: generate_basis.py 
 * 
 * NOTE: Any modifications made in this file will be lost next time 
 *       generate_basis.py is run.
 */

#include "../basis_set_list.hpp"
# include <chemist/chemist.hpp>

namespace chemcache::basis_sets {

void load_six_dash_31g_oparen_3df_comma_3pd_cparen_(chemist::BasisSetManager& bsm) {
    chemist::BasisSetManager::ao_basis_map basis_map;

    // Z = 1
    basis_map.emplace(1, chemist::Center<double>(0.0, 0.0, 0.0));

    basis_map.at(1).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{3.3494600000e-02,2.3472695000e-01,8.1375733000e-01},
        std::vector<double>{1.8731137000e+01,2.8253937000e+00,6.4012170000e-01});
    basis_map.at(1).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.6127780000e-01});
    basis_map.at(1).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.0000000000e+00});
    basis_map.at(1).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{7.5000000000e-01});
    basis_map.at(1).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.8750000000e-01});
    basis_map.at(1).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.0000000000e+00});

    // Z = 2
    basis_map.emplace(2, chemist::Center<double>(0.0, 0.0, 0.0));

    basis_map.at(2).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{2.3766000000e-02,1.5467900000e-01,4.6963000000e-01},
        std::vector<double>{3.8421634000e+01,5.7780300000e+00,1.2417740000e+00});
    basis_map.at(2).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.9796400000e-01});
    basis_map.at(2).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.0000000000e+00});
    basis_map.at(2).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{7.5000000000e-01});
    basis_map.at(2).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.8750000000e-01});
    basis_map.at(2).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.0000000000e+00});

    // Z = 3
    basis_map.emplace(3, chemist::Center<double>(0.0, 0.0, 0.0));

    basis_map.at(3).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{2.1426000000e-03,1.6208900000e-02,7.7315600000e-02,2.4578600000e-01,4.7018900000e-01,3.4547080000e-01},
        std::vector<double>{6.4241892000e+02,9.6798515000e+01,2.2091121000e+01,6.2010703000e+00,1.9351177000e+00,6.3673580000e-01});
    basis_map.at(3).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{-3.5091700000e-02,-1.9123280000e-01,1.0839878000e+00},
        std::vector<double>{2.3249184000e+00,6.3243060000e-01,7.9053400000e-02});
    basis_map.at(3).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{8.9415000000e-03,1.4100950000e-01,9.4536370000e-01},
        std::vector<double>{2.3249184000e+00,6.3243060000e-01,7.9053400000e-02});
    basis_map.at(3).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.5962000000e-02});
    basis_map.at(3).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.5962000000e-02});
    basis_map.at(3).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{8.0000000000e-01});
    basis_map.at(3).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.0000000000e-01});
    basis_map.at(3).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{5.0000000000e-02});
    basis_map.at(3).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.5000000000e-01});

    // Z = 4
    basis_map.emplace(4, chemist::Center<double>(0.0, 0.0, 0.0));

    basis_map.at(4).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.9448000000e-03,1.4835100000e-02,7.2090600000e-02,2.3715420000e-01,4.6919870000e-01,3.5652020000e-01},
        std::vector<double>{1.2645857000e+03,1.8993681000e+02,4.3159089000e+01,1.2098663000e+01,3.8063232000e+00,1.2728903000e+00});
    basis_map.at(4).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{-1.1264870000e-01,-2.2950640000e-01,1.1869167000e+00},
        std::vector<double>{3.1964631000e+00,7.4781330000e-01,2.1996630000e-01});
    basis_map.at(4).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{5.5980200000e-02,2.6155060000e-01,7.9397230000e-01},
        std::vector<double>{3.1964631000e+00,7.4781330000e-01,2.1996630000e-01});
    basis_map.at(4).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{8.2309900000e-02});
    basis_map.at(4).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{8.2309900000e-02});
    basis_map.at(4).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.0200000000e+00});
    basis_map.at(4).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.5500000000e-01});
    basis_map.at(4).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{6.3750000000e-02});
    basis_map.at(4).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.6000000000e-01});

    // Z = 5
    basis_map.emplace(5, chemist::Center<double>(0.0, 0.0, 0.0));

    basis_map.at(5).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.8663000000e-03,1.4251500000e-02,6.9551600000e-02,2.3257290000e-01,4.6707870000e-01,3.6343140000e-01},
        std::vector<double>{2.0688823000e+03,3.1064957000e+02,7.0683033000e+01,1.9861080000e+01,6.2993048000e+00,2.1270270000e+00});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{-1.3039380000e-01,-1.3078890000e-01,1.1309444000e+00},
        std::vector<double>{4.7279710000e+00,1.1903377000e+00,3.5941170000e-01});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{7.4597600000e-02,3.0784670000e-01,7.4345680000e-01},
        std::vector<double>{4.7279710000e+00,1.1903377000e+00,3.5941170000e-01});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.2675120000e-01});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.2675120000e-01});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.6040000000e+00});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.0100000000e-01});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.0025000000e-01});
    basis_map.at(5).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{5.0000000000e-01});

    // Z = 6
    basis_map.emplace(6, chemist::Center<double>(0.0, 0.0, 0.0));

    basis_map.at(6).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.8347000000e-03,1.4037300000e-02,6.8842600000e-02,2.3218440000e-01,4.6794130000e-01,3.6231200000e-01},
        std::vector<double>{3.0475249000e+03,4.5736951000e+02,1.0394869000e+02,2.9210155000e+01,9.2866630000e+00,3.1639270000e+00});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{-1.1933240000e-01,-1.6085420000e-01,1.1434564000e+00},
        std::vector<double>{7.8682724000e+00,1.8812885000e+00,5.4424930000e-01});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{6.8999100000e-02,3.1642400000e-01,7.4430830000e-01},
        std::vector<double>{7.8682724000e+00,1.8812885000e+00,5.4424930000e-01});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.6871440000e-01});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.6871440000e-01});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.5040000000e+00});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{6.2600000000e-01});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.5650000000e-01});
    basis_map.at(6).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{8.0000000000e-01});

    // Z = 7
    basis_map.emplace(7, chemist::Center<double>(0.0, 0.0, 0.0));

    basis_map.at(7).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.8348000000e-03,1.3995000000e-02,6.8587000000e-02,2.3224100000e-01,4.6907000000e-01,3.6045500000e-01},
        std::vector<double>{4.1735110000e+03,6.2745790000e+02,1.4290210000e+02,4.0234330000e+01,1.2820210000e+01,4.3904370000e+00});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{-1.1496100000e-01,-1.6911800000e-01,1.1458520000e+00},
        std::vector<double>{1.1626358000e+01,2.7162800000e+00,7.7221800000e-01});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{6.7580000000e-02,3.2390700000e-01,7.4089500000e-01},
        std::vector<double>{1.1626358000e+01,2.7162800000e+00,7.7221800000e-01});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.1203130000e-01});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.1203130000e-01});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.6520000000e+00});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{9.1300000000e-01});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.2825000000e-01});
    basis_map.at(7).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.0000000000e+00});

    // Z = 8
    basis_map.emplace(8, chemist::Center<double>(0.0, 0.0, 0.0));

    basis_map.at(8).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.8311000000e-03,1.3950100000e-02,6.8445100000e-02,2.3271430000e-01,4.7019300000e-01,3.5852090000e-01},
        std::vector<double>{5.4846717000e+03,8.2523495000e+02,1.8804696000e+02,5.2964500000e+01,1.6897570000e+01,5.7996353000e+00});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{-1.1077750000e-01,-1.4802630000e-01,1.1307670000e+00},
        std::vector<double>{1.5539616000e+01,3.5999336000e+00,1.0137618000e+00});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{7.0874300000e-02,3.3975280000e-01,7.2715860000e-01},
        std::vector<double>{1.5539616000e+01,3.5999336000e+00,1.0137618000e+00});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.7000580000e-01});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.7000580000e-01});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{5.1600000000e+00});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.2920000000e+00});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.2250000000e-01});
    basis_map.at(8).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.4000000000e+00});

    // Z = 9
    basis_map.emplace(9, chemist::Center<double>(0.0, 0.0, 0.0));

    basis_map.at(9).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.8196169000e-03,1.3916079600e-02,6.8405324500e-02,2.3318576000e-01,4.7126743900e-01,3.5661854600e-01},
        std::vector<double>{7.0017130900e+03,1.0513660900e+03,2.3928569000e+02,6.7397445300e+01,2.1519957300e+01,7.4031013000e+00});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{-1.0850697500e-01,-1.4645165800e-01,1.1286885800e+00},
        std::vector<double>{2.0847952800e+01,4.8083083400e+00,1.3440698600e+00});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{7.1628724300e-02,3.4591210300e-01,7.2246995700e-01},
        std::vector<double>{2.0847952800e+01,4.8083083400e+00,1.3440698600e+00});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.5815139300e-01});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.5815139300e-01});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{7.0000000000e+00});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.7500000000e+00});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.3750000000e-01});
    basis_map.at(9).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.8500000000e+00});

    // Z = 10
    basis_map.emplace(10, chemist::Center<double>(0.0, 0.0, 0.0));

    basis_map.at(10).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.8843481000e-03,1.4336899400e-02,7.0109623300e-02,2.3737326600e-01,4.7300712600e-01,3.4840124100e-01},
        std::vector<double>{8.4258515300e+03,1.2685194000e+03,2.8962141400e+02,8.1859004000e+01,2.6251507900e+01,9.0947205100e+00});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{-1.0711828700e-01,-1.4616382100e-01,1.1277735000e+00},
        std::vector<double>{2.6532131000e+01,6.1017550100e+00,1.6962715300e+00});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{7.1909588500e-02,3.4951337200e-01,7.1994051200e-01},
        std::vector<double>{2.6532131000e+01,6.1017550100e+00,1.6962715300e+00});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.4581870000e-01});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.4581870000e-01});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{9.2160000000e+00});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.3040000000e+00});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{5.7600000000e-01});
    basis_map.at(10).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.5000000000e+00});

    // Z = 11
    basis_map.emplace(11, chemist::Center<double>(0.0, 0.0, 0.0));

    basis_map.at(11).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.9377000000e-03,1.4807000000e-02,7.2706000000e-02,2.5262900000e-01,4.9324200000e-01,3.1316900000e-01},
        std::vector<double>{9.9932000000e+03,1.4998900000e+03,3.4195100000e+02,9.4679700000e+01,2.9734500000e+01,1.0006300000e+01});
    basis_map.at(11).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{-3.5421000000e-03,-4.3959000000e-02,-1.0975210000e-01,1.8739800000e-01,6.4669900000e-01,3.0605800000e-01},
        std::vector<double>{1.5096300000e+02,3.5587800000e+01,1.1168300000e+01,3.9020100000e+00,1.3817700000e+00,4.6638200000e-01});
    basis_map.at(11).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{5.0017000000e-03,3.5511000000e-02,1.4282500000e-01,3.3862000000e-01,4.5157900000e-01,2.7327100000e-01},
        std::vector<double>{1.5096300000e+02,3.5587800000e+01,1.1168300000e+01,3.9020100000e+00,1.3817700000e+00,4.6638200000e-01});
    basis_map.at(11).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{-2.4850300000e-01,-1.3170400000e-01,1.2335200000e+00},
        std::vector<double>{4.9796600000e-01,8.4353000000e-02,6.6635000000e-02});
    basis_map.at(11).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{-2.3023000000e-02,9.5035900000e-01,5.9858000000e-02},
        std::vector<double>{4.9796600000e-01,8.4353000000e-02,6.6635000000e-02});
    basis_map.at(11).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.5954400000e-02});
    basis_map.at(11).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.5954400000e-02});
    basis_map.at(11).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{7.0000000000e-01});
    basis_map.at(11).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.7500000000e-01});
    basis_map.at(11).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.3750000000e-02});
    basis_map.at(11).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.5000000000e-01});

    // Z = 12
    basis_map.emplace(12, chemist::Center<double>(0.0, 0.0, 0.0));

    basis_map.at(12).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.9778000000e-03,1.5114000000e-02,7.3911000000e-02,2.4919100000e-01,4.8792800000e-01,3.1966200000e-01},
        std::vector<double>{1.1722800000e+04,1.7599300000e+03,4.0084600000e+02,1.1280700000e+02,3.5999700000e+01,1.2182800000e+01});
    basis_map.at(12).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{-3.2372000000e-03,-4.1008000000e-02,-1.1260000000e-01,1.4863300000e-01,6.1649700000e-01,3.6482900000e-01},
        std::vector<double>{1.8918000000e+02,4.5211900000e+01,1.4356300000e+01,5.1388600000e+00,1.9065200000e+00,7.0588700000e-01});
    basis_map.at(12).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{4.9281000000e-03,3.4989000000e-02,1.4072500000e-01,3.3364200000e-01,4.4494000000e-01,2.6925400000e-01},
        std::vector<double>{1.8918000000e+02,4.5211900000e+01,1.4356300000e+01,5.1388600000e+00,1.9065200000e+00,7.0588700000e-01});
    basis_map.at(12).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{-2.1229000000e-01,-1.0798500000e-01,1.1758400000e+00},
        std::vector<double>{9.2934000000e-01,2.6903500000e-01,1.1737900000e-01});
    basis_map.at(12).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{-2.2419000000e-02,1.9227000000e-01,8.4618100000e-01},
        std::vector<double>{9.2934000000e-01,2.6903500000e-01,1.1737900000e-01});
    basis_map.at(12).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.2106100000e-02});
    basis_map.at(12).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.2106100000e-02});
    basis_map.at(12).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{7.0000000000e-01});
    basis_map.at(12).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.7500000000e-01});
    basis_map.at(12).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.3750000000e-02});
    basis_map.at(12).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.0000000000e-01});

    // Z = 13
    basis_map.emplace(13, chemist::Center<double>(0.0, 0.0, 0.0));

    basis_map.at(13).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.9426700000e-03,1.4859900000e-02,7.2849400000e-02,2.4683000000e-01,4.8725800000e-01,3.2349600000e-01},
        std::vector<double>{1.3983100000e+04,2.0987500000e+03,4.7770500000e+02,1.3436000000e+02,4.2870900000e+01,1.4518900000e+01});
    basis_map.at(13).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{-2.9261900000e-03,-3.7408000000e-02,-1.1448700000e-01,1.1563500000e-01,6.1259500000e-01,3.9379900000e-01},
        std::vector<double>{2.3966800000e+02,5.7441900000e+01,1.8285900000e+01,6.5991400000e+00,2.4904900000e+00,9.4454000000e-01});
    basis_map.at(13).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{4.6028500000e-03,3.3199000000e-02,1.3628200000e-01,3.3047600000e-01,4.4914600000e-01,2.6570400000e-01},
        std::vector<double>{2.3966800000e+02,5.7441900000e+01,1.8285900000e+01,6.5991400000e+00,2.4904900000e+00,9.4454000000e-01});
    basis_map.at(13).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{-2.2760600000e-01,1.4458300000e-03,1.0927900000e+00},
        std::vector<double>{1.2779000000e+00,3.9759000000e-01,1.6009500000e-01});
    basis_map.at(13).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{-1.7513000000e-02,2.4453300000e-01,8.0493400000e-01},
        std::vector<double>{1.2779000000e+00,3.9759000000e-01,1.6009500000e-01});
    basis_map.at(13).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{5.5657700000e-02});
    basis_map.at(13).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{5.5657700000e-02});
    basis_map.at(13).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.3000000000e+00});
    basis_map.at(13).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.2500000000e-01});
    basis_map.at(13).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{8.1250000000e-02});
    basis_map.at(13).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.5000000000e-01});

    // Z = 14
    basis_map.emplace(14, chemist::Center<double>(0.0, 0.0, 0.0));

    basis_map.at(14).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.9594800000e-03,1.4928800000e-02,7.2847800000e-02,2.4613000000e-01,4.8591400000e-01,3.2500200000e-01},
        std::vector<double>{1.6115900000e+04,2.4255800000e+03,5.5386700000e+02,1.5634000000e+02,5.0068300000e+01,1.7017800000e+01});
    basis_map.at(14).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{-2.7809400000e-03,-3.5714600000e-02,-1.1498500000e-01,9.3563400000e-02,6.0301700000e-01,4.1895900000e-01},
        std::vector<double>{2.9271800000e+02,6.9873100000e+01,2.2336300000e+01,8.1503900000e+00,3.1345800000e+00,1.2254300000e+00});
    basis_map.at(14).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{4.4382600000e-03,3.2667900000e-02,1.3472100000e-01,3.2867800000e-01,4.4964000000e-01,2.6137200000e-01},
        std::vector<double>{2.9271800000e+02,6.9873100000e+01,2.2336300000e+01,8.1503900000e+00,3.1345800000e+00,1.2254300000e+00});
    basis_map.at(14).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{-2.4463000000e-01,4.3157200000e-03,1.0981800000e+00},
        std::vector<double>{1.7273800000e+00,5.7292200000e-01,2.2219200000e-01});
    basis_map.at(14).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{-1.7795100000e-02,2.5353900000e-01,8.0066900000e-01},
        std::vector<double>{1.7273800000e+00,5.7292200000e-01,2.2219200000e-01});
    basis_map.at(14).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{7.7836900000e-02});
    basis_map.at(14).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{7.7836900000e-02});
    basis_map.at(14).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.8000000000e+00});
    basis_map.at(14).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.5000000000e-01});
    basis_map.at(14).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.1250000000e-01});
    basis_map.at(14).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.2000000000e-01});

    // Z = 15
    basis_map.emplace(15, chemist::Center<double>(0.0, 0.0, 0.0));

    basis_map.at(15).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.8516000000e-03,1.4206200000e-02,6.9999500000e-02,2.4007900000e-01,4.8476200000e-01,3.3520000000e-01},
        std::vector<double>{1.9413300000e+04,2.9094200000e+03,6.6136400000e+02,1.8575900000e+02,5.9194300000e+01,2.0031000000e+01});
    basis_map.at(15).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{-2.7821700000e-03,-3.6049900000e-02,-1.1663100000e-01,9.6832800000e-02,6.1441800000e-01,4.0379800000e-01},
        std::vector<double>{3.3947800000e+02,8.1010100000e+01,2.5878000000e+01,9.4522100000e+00,3.6656600000e+00,1.4674600000e+00});
    basis_map.at(15).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{4.5646200000e-03,3.3693600000e-02,1.3975500000e-01,3.3936200000e-01,4.5092100000e-01,2.3858600000e-01},
        std::vector<double>{3.3947800000e+02,8.1010100000e+01,2.5878000000e+01,9.4522100000e+00,3.6656600000e+00,1.4674600000e+00});
    basis_map.at(15).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{-2.5292300000e-01,3.2851700000e-02,1.0812500000e+00},
        std::vector<double>{2.1562300000e+00,7.4899700000e-01,2.8314500000e-01});
    basis_map.at(15).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{-1.7765300000e-02,2.7405800000e-01,7.8542100000e-01},
        std::vector<double>{2.1562300000e+00,7.4899700000e-01,2.8314500000e-01});
    basis_map.at(15).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{9.9831700000e-02});
    basis_map.at(15).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{9.9831700000e-02});
    basis_map.at(15).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.2000000000e+00});
    basis_map.at(15).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{5.5000000000e-01});
    basis_map.at(15).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.3750000000e-01});
    basis_map.at(15).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{4.5000000000e-01});

    // Z = 16
    basis_map.emplace(16, chemist::Center<double>(0.0, 0.0, 0.0));

    basis_map.at(16).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.8690000000e-03,1.4230000000e-02,6.9696000000e-02,2.3848700000e-01,4.8330700000e-01,3.3807400000e-01},
        std::vector<double>{2.1917100000e+04,3.3014900000e+03,7.5414600000e+02,2.1271100000e+02,6.7989600000e+01,2.3051500000e+01});
    basis_map.at(16).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{-2.3767000000e-03,-3.1693000000e-02,-1.1331700000e-01,5.6090000000e-02,5.9225500000e-01,4.5500600000e-01},
        std::vector<double>{4.2373500000e+02,1.0071000000e+02,3.2159900000e+01,1.1807900000e+01,4.6311000000e+00,1.8702500000e+00});
    basis_map.at(16).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{4.0610000000e-03,3.0681000000e-02,1.3045200000e-01,3.2720500000e-01,4.5285100000e-01,2.5604200000e-01},
        std::vector<double>{4.2373500000e+02,1.0071000000e+02,3.2159900000e+01,1.1807900000e+01,4.6311000000e+00,1.8702500000e+00});
    basis_map.at(16).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{-2.5037400000e-01,6.6957000000e-02,1.0545100000e+00},
        std::vector<double>{2.6158400000e+00,9.2216700000e-01,3.4128700000e-01});
    basis_map.at(16).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{-1.4511000000e-02,3.1026300000e-01,7.5448300000e-01},
        std::vector<double>{2.6158400000e+00,9.2216700000e-01,3.4128700000e-01});
    basis_map.at(16).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.1716700000e-01});
    basis_map.at(16).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.1716700000e-01});
    basis_map.at(16).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.6000000000e+00});
    basis_map.at(16).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{6.5000000000e-01});
    basis_map.at(16).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.6250000000e-01});
    basis_map.at(16).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{5.5000000000e-01});

    // Z = 17
    basis_map.emplace(17, chemist::Center<double>(0.0, 0.0, 0.0));

    basis_map.at(17).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.8330000000e-03,1.4034000000e-02,6.9097000000e-02,2.3745200000e-01,4.8303400000e-01,3.3985600000e-01},
        std::vector<double>{2.5180100000e+04,3.7803500000e+03,8.6047400000e+02,2.4214500000e+02,7.7334900000e+01,2.6247000000e+01});
    basis_map.at(17).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{-2.2974000000e-03,-3.0714000000e-02,-1.1252800000e-01,4.5016000000e-02,5.8935300000e-01,4.6520600000e-01},
        std::vector<double>{4.9176500000e+02,1.1698400000e+02,3.7415300000e+01,1.3783400000e+01,5.4521500000e+00,2.2258800000e+00});
    basis_map.at(17).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{3.9894000000e-03,3.0318000000e-02,1.2988000000e-01,3.2795100000e-01,4.5352700000e-01,2.5215400000e-01},
        std::vector<double>{4.9176500000e+02,1.1698400000e+02,3.7415300000e+01,1.3783400000e+01,5.4521500000e+00,2.2258800000e+00});
    basis_map.at(17).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{-2.5183000000e-01,6.1589000000e-02,1.0601800000e+00},
        std::vector<double>{3.1864900000e+00,1.1442700000e+00,4.2037700000e-01});
    basis_map.at(17).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{-1.4299000000e-02,3.2357200000e-01,7.4350700000e-01},
        std::vector<double>{3.1864900000e+00,1.1442700000e+00,4.2037700000e-01});
    basis_map.at(17).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.4265700000e-01});
    basis_map.at(17).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.4265700000e-01});
    basis_map.at(17).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.0000000000e+00});
    basis_map.at(17).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{7.5000000000e-01});
    basis_map.at(17).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.8750000000e-01});
    basis_map.at(17).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{7.0000000000e-01});

    // Z = 18
    basis_map.emplace(18, chemist::Center<double>(0.0, 0.0, 0.0));

    basis_map.at(18).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.8252600000e-03,1.3968600000e-02,6.8707300000e-02,2.3620400000e-01,4.8221400000e-01,3.4204300000e-01},
        std::vector<double>{2.8348300000e+04,4.2576200000e+03,9.6985700000e+02,2.7326300000e+02,8.7369500000e+01,2.9686700000e+01});
    basis_map.at(18).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{-2.1597200000e-03,-2.9077500000e-02,-1.1082700000e-01,2.7699900000e-02,5.7761300000e-01,4.8868800000e-01},
        std::vector<double>{5.7589100000e+02,1.3681600000e+02,4.3809800000e+01,1.6209400000e+01,6.4608400000e+00,2.6511400000e+00});
    basis_map.at(18).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{3.8066500000e-03,2.9230500000e-02,1.2646700000e-01,3.2351000000e-01,4.5489600000e-01,2.5663000000e-01},
        std::vector<double>{5.7589100000e+02,1.3681600000e+02,4.3809800000e+01,1.6209400000e+01,6.4608400000e+00,2.6511400000e+00});
    basis_map.at(18).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{-2.5559200000e-01,3.7806600000e-02,1.0805600000e+00},
        std::vector<double>{3.8602800000e+00,1.4137300000e+00,5.1664600000e-01});
    basis_map.at(18).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{-1.5919700000e-02,3.2464600000e-01,7.4399000000e-01},
        std::vector<double>{3.8602800000e+00,1.4137300000e+00,5.1664600000e-01});
    basis_map.at(18).add_shell(chemist::ShellType::pure, 0,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.7388800000e-01});
    basis_map.at(18).add_shell(chemist::ShellType::pure, 1,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{1.7388800000e-01});
    basis_map.at(18).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{3.4000000000e+00});
    basis_map.at(18).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{8.5000000000e-01});
    basis_map.at(18).add_shell(chemist::ShellType::pure, 2,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{2.1250000000e-01});
    basis_map.at(18).add_shell(chemist::ShellType::pure, 3,
        std::vector<double>{1.0000000000e+00},
        std::vector<double>{8.5000000000e-01});


    bsm.insert("6-31g(3df,3pd)", basis_map);
} // function load_six_dash_31g_oparen_3df_comma_3pd_cparen_

} // namespace chemcache::basis_sets
