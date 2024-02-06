Search.setIndex({"docnames": ["autoapi/data_management/generate_atomicinfo/index", "autoapi/data_management/generate_basis/index", "autoapi/data_management/generate_densities/index", "autoapi/data_management/generate_elec_configs/index", "autoapi/data_management/generate_molecules/index", "autoapi/data_management/helper_fxns/index", "autoapi/data_management/index", "autoapi/data_management/scrape_bse/index", "design", "index", "installation", "module_api/Atom", "module_api/Atom_Isotope", "module_api/NWX_Molecules", "module_api/Symbol_from_Z", "module_api/Z_from_Symbol", "module_api/index", "module_api/sto-3g", "module_api/sto-3g_atomic_basis"], "filenames": ["autoapi/data_management/generate_atomicinfo/index.rst", "autoapi/data_management/generate_basis/index.rst", "autoapi/data_management/generate_densities/index.rst", "autoapi/data_management/generate_elec_configs/index.rst", "autoapi/data_management/generate_molecules/index.rst", "autoapi/data_management/helper_fxns/index.rst", "autoapi/data_management/index.rst", "autoapi/data_management/scrape_bse/index.rst", "design.rst", "index.rst", "installation.rst", "module_api/Atom.rst", "module_api/Atom_Isotope.rst", "module_api/NWX_Molecules.rst", "module_api/Symbol_from_Z.rst", "module_api/Z_from_Symbol.rst", "module_api/index.rst", "module_api/sto-3g.rst", "module_api/sto-3g_atomic_basis.rst"], "titles": ["<code class=\"xref py py-mod docutils literal notranslate\"><span class=\"pre\">data_management.generate_atomicinfo</span></code>", "<code class=\"xref py py-mod docutils literal notranslate\"><span class=\"pre\">data_management.generate_basis</span></code>", "<code class=\"xref py py-mod docutils literal notranslate\"><span class=\"pre\">data_management.generate_densities</span></code>", "<code class=\"xref py py-mod docutils literal notranslate\"><span class=\"pre\">data_management.generate_elec_configs</span></code>", "<code class=\"xref py py-mod docutils literal notranslate\"><span class=\"pre\">data_management.generate_molecules</span></code>", "<code class=\"xref py py-mod docutils literal notranslate\"><span class=\"pre\">data_management.helper_fxns</span></code>", "<code class=\"xref py py-mod docutils literal notranslate\"><span class=\"pre\">data_management</span></code>", "<code class=\"xref py py-mod docutils literal notranslate\"><span class=\"pre\">data_management.scrape_bse</span></code>", "ChemCache Design", "ChemCache", "Installation", "Atom", "Atom Isotope", "NWX Molecules", "Symbol from Z", "Z from Symbol", "Modules API", "sto-3g", "sto-3g atomic basis"], "terms": {"thi": [0, 1, 2, 3, 4, 5, 7, 8, 10, 11, 12, 13, 14, 15, 17, 18], "script": [0, 1, 2, 3, 4, 5, 6, 7, 8], "i": [0, 1, 2, 3, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 17, 18], "us": [0, 1, 3, 5, 7, 8, 10, 11, 12, 13, 14, 15, 17, 18], "creat": [0, 1, 2, 3, 4, 7, 10], "experiment": 0, "data": [0, 1, 2, 3, 4, 6, 7, 8], "look": [0, 3, 5], "up": [0, 5], "tabl": [0, 11, 12, 13, 14, 15, 17, 18], "atom": [0, 1, 2, 3, 4, 7, 9, 13, 15, 16, 17], "origin": 0, "author": 0, "ben": 0, "pritchard": 0, "modifi": 0, "zacheri": 0, "crandal": 0, "In": 0, "order": 0, "run": [0, 1, 2, 4, 5, 10, 11, 12, 13, 14, 15, 17, 18], "need": [0, 6, 7, 10], "know": 0, "locat": [0, 8, 10], "directori": [0, 1, 2, 3, 4, 5, 7, 8, 10], "read": [0, 2, 4], "from": [0, 1, 2, 3, 4, 7, 8, 9, 10, 16], "src": [0, 2, 4, 8], "code": [0, 6], "output": [0, 3, 5, 7], "result": [0, 1], "One": 0, "can": [0, 1, 2, 4, 7, 8, 10, 11, 12, 13, 14, 15, 17, 18], "also": [0, 6, 10], "option": [0, 1, 2, 3, 4, 5, 7, 9, 11, 12, 13, 14, 15, 17, 18], "provid": [0, 4, 5, 7, 8, 11, 12, 13, 14, 15, 17, 18], "nondefault": 0, "valu": [0, 1, 2, 3, 4, 7, 10, 11, 12, 13, 14, 15, 17, 18], "electron": [0, 3], "mass": 0, "dalton": 0, "ratio": [0, 4], "For": [0, 3, 5, 7, 13, 17], "readabl": [0, 3, 11, 12, 13, 14, 15, 17, 18], "conveni": [0, 3], "we": [0, 3, 11, 12, 13, 14, 15, 17, 18], "few": [0, 3], "abbrevi": [0, 3], "throughout": [0, 3, 5], "z": [0, 1, 2, 3, 4, 9, 16], "number": [0, 1, 2, 3, 4, 5, 7, 15, 18], "an": [0, 1, 3, 4, 5, 7, 11, 12, 13, 18], "sym": [0, 3], "symbol": [0, 1, 2, 3, 4, 5, 9, 16], "e": [0, 3, 8], "g": [0, 3, 7, 8], "h": [0, 1, 2, 3, 4, 7, 10], "hydrogen": [0, 3], "he": [0, 3], "helium": [0, 3], "py": [0, 1, 2, 3, 4, 7, 8], "amu2m": 0, "data_dir": [0, 3], "src_dir": [0, 1, 2, 3, 4], "posit": [0, 1, 2, 3, 4, 7], "argument": [0, 1, 2, 3, 4, 7, 11, 12, 13, 14, 15, 17, 18], "inform": [0, 3], "file": [0, 1, 2, 3, 4, 5, 7, 9, 10], "destin": [0, 1, 2, 3, 4, 7], "gener": [0, 1, 2, 3, 4, 5, 6, 7, 9, 10, 17], "sourc": [0, 1, 2, 3, 4, 6, 9], "help": [0, 1, 2, 3, 4, 7], "show": [0, 1, 2, 3, 4, 7], "messag": [0, 1, 2, 3, 4, 7], "exit": [0, 1, 2, 3, 4, 7], "one": [0, 7, 11, 12, 13, 14, 15, 17, 18], "default": [0, 1, 2, 4, 5, 7, 10, 11, 12, 13, 14, 15, 17, 18], "1822": 0, "888486192": 0, "follow": [0, 1, 2, 3, 4, 7, 11, 12, 13, 14, 15, 17, 18], "": [0, 1, 3, 10, 11, 12, 13, 14, 15, 17, 18], "elementnam": [0, 1, 2, 3, 4], "txt": [0, 1, 2, 3, 4, 8], "ciaaw": 0, "isotopemass": 0, "load_el": 0, "hpp": [0, 1, 2], "atomicdata": [0, 3], "add_isotop": 0, "num": 0, "int": [0, 3, 4], "float": [0, 4], "none": [0, 1, 2, 3, 4, 7], "add": [0, 1, 3, 4, 7], "new": [0, 2, 7], "isotop": [0, 9, 16], "list": [0, 1, 2, 3, 4, 5, 7, 10, 11, 12, 13, 14, 15, 17, 18], "element": [0, 1, 2, 3, 7], "paramet": [0, 1, 2, 3, 4, 5, 7], "neutron": 0, "n": [0, 2, 3, 11, 12, 13, 14, 15, 17, 18], "__repr__": [0, 3, 4], "return": [0, 1, 2, 3, 4, 5, 7, 11, 12, 13, 14, 15, 17, 18], "format": [0, 1, 2, 3, 4, 5, 7], "text": [0, 4, 5, 7], "represent": [0, 1, 4], "type": [0, 1, 2, 3, 4, 5, 7, 8], "str": [0, 1, 2, 3, 4, 5, 7], "parse_symbol": [0, 3], "name_fil": [0, 3], "dict": [0, 1, 2, 3, 4, 5, 7], "pars": [0, 1, 2, 3, 4, 5, 7], "given": [0, 1, 2, 3, 4, 5, 7, 11, 12, 13, 14, 15, 17, 18], "them": [0, 1, 3, 10], "exist": [0, 3, 7, 11, 12, 13, 14, 15, 17, 18], "collect": [0, 1, 2, 3, 4, 7], "name": [0, 1, 3, 4, 5, 7, 10, 11, 12, 13, 14, 15, 17, 18], "current": [0, 3, 4, 6, 7], "load": [0, 3], "ad": [0, 3, 5], "here": [0, 3, 6, 10], "_parse_ciaaw_isotop": 0, "iso_fil": 0, "commiss": 0, "abund": 0, "weight": 0, "_parse_ciaww_mass": 0, "mass_fil": 0, "_write_z_from_sym": 0, "out_dir": [0, 3], "z_from_sym": 0, "cpp": [0, 1, 2, 3, 4], "header": [0, 3, 7, 10], "_write_sym_from_z": 0, "sym_from_z": 0, "_write_atoms_averag": 0, "atoms_averag": 0, "_write_atoms_isotop": 0, "atoms_isotop": 0, "main": [0, 1, 2, 3, 4, 7], "arg": [0, 1, 2, 3, 4, 7], "argpars": [0, 1, 2, 3, 4, 7], "namespac": [0, 1, 2, 3, 4, 7], "entri": [0, 1, 2, 3, 4, 7], "point": [0, 1, 2, 3, 4, 7, 10], "info": [0, 3], "command": [0, 1, 2, 3, 4, 7, 9], "line": [0, 1, 2, 3, 4, 5, 7], "parse_arg": [0, 1, 2, 3, 4, 7], "loop": 1, "over": 1, "seri": 1, "basi": [1, 2, 5, 7, 8, 9, 10, 16], "set": [1, 2, 5, 7, 8, 10, 11, 12, 13, 14, 15], "write": [1, 2, 4, 5, 7], "out": [1, 7], "fill": 1, "The": [1, 2, 4, 6, 7, 8, 10, 11, 12, 13, 14, 15, 17, 18], "suitabl": 1, "basissetexchang": [1, 7], "r": [1, 2, 4, 8], "atoms_dir": [1, 2, 4], "basis_set_sourc": 1, "If": [1, 2, 3, 4, 7, 10], "combin": [1, 2, 4], "flag": [1, 2, 4], "recurs": [1, 2, 4, 5], "search": [1, 2, 4, 5], "toggl": [1, 2, 4, 7], "through": [1, 2, 4], "off": [1, 2, 4, 7, 10], "path": [1, 2, 3, 4, 5, 10], "where": [1, 2, 3, 4, 5, 8, 10], "found": [1, 2, 4, 5, 10], "base": [1, 2, 4, 7, 8], "includ": [1, 2, 4, 7, 8], "ar": [1, 2, 4, 6, 7, 10, 11, 12, 13, 14, 15, 17, 18], "must": [1, 2, 3, 4, 7, 10, 11, 12, 13, 14, 15, 17, 18], "present": [1, 2, 4], "befor": [1, 2, 3, 4, 11, 12, 13, 14, 15, 17, 18], "all_basis_set_fil": [1, 2, 7], "basis_set_list": 1, "load_basis_set": 1, "shell": [1, 3, 17, 18], "l": [1, 3], "num_format": 1, "10e": 1, "is_pur": 1, "true": [1, 7, 11, 12, 13, 14, 15, 17, 18], "repres": [1, 2, 3, 4], "add_prim": 1, "exp": 1, "coef": 1, "primit": [1, 17, 18], "expon": 1, "contract": [1, 7], "coeffici": 1, "cxxifi": [1, 4], "center": 1, "tab": [1, 2, 4], "c": [1, 2, 4, 5, 6, 9, 11, 12, 13, 14, 15, 17, 18], "chemist": [1, 11, 12, 13, 17, 18], "string": [1, 2, 4, 7], "_write_basis_fil": 1, "out_fil": [1, 2], "bs_name": [1, 2, 5], "basis_set": [1, 2, 8, 17, 18], "_write_bas": 1, "_parse_bases_gb": 1, "basis_set_filenam": 1, "sym2z": [1, 2, 4], "l2num": 1, "filepath": [1, 2, 4], "basis_set_filepath": 1, "full": [1, 2, 3, 4, 5, 7, 11, 12, 13, 14, 15, 17, 18], "dictionari": [1, 3, 4, 5, 7], "associ": [1, 8, 14, 18], "orbit": [1, 7], "letter": 1, "support": 1, "each": [1, 3, 4, 7, 8, 11, 12, 13, 14, 15, 17, 18], "_parse_bases_nw": 1, "_parse_bas": 1, "nwchem": [1, 7], "specifi": [1, 2, 4, 7], "redirect": 1, "correct": 1, "structur": 1, "kei": [1, 3, 7, 11, 12, 13, 14, 15, 17, 18], "effect": 1, "make": [1, 4, 10], "map": [1, 2, 4, 5], "anoth": 1, "basis_nam": [1, 7], "access": [1, 3, 6, 11, 12, 13, 14, 15, 17, 18], "list_of_shel": 1, "lowercas": [1, 2, 4], "convers": 1, "p": 1, "d": 1, "f": 1, "etc": 1, "correspond": [1, 5, 8, 13], "0": 1, "1": [1, 3, 4], "2": [1, 10], "3": [1, 3, 10], "rais": [1, 2, 4, 7], "runtimeerror": [1, 2, 4, 7], "unsupport": [1, 2, 4], "densiti": [2, 4], "atomic_density_dir": 2, "nwx_atomic_dens": 2, "atomic_dens": 2, "made": [2, 8], "manual": [2, 11, 12, 13, 14, 15, 17, 18], "python": [2, 8], "add_dens": 2, "cmake": 2, "make_square_arr": 2, "spacer": 2, "_write_den_fil": 2, "_write_dens": 2, "_parse_densities_xml": 2, "xml": 2, "deprec": 2, "adopt": 2, "y": [2, 10], "02": 2, "23": 2, "sort": 2, "_parse_densities_dat": 2, "dat": 2, "_parse_dens": 2, "extens": [2, 4, 5, 7], "contain": [3, 5, 6, 7, 8], "ground": 3, "state": [3, 10], "configur": [3, 9], "periodict": 3, "object": 3, "generate_ptable_config": 3, "nist": 3, "atomicion": 3, "atomic_configur": 3, "atomconfig": 3, "lmax": 3, "nmax": 3, "7": 3, "i_to_lchar": 3, "lchar_to_i": 3, "properti": 3, "config_ful": 3, "index": 3, "config": 3, "reduc": 3, "per": 3, "np": 3, "nd": 3, "tupl": [3, 7], "repr": 3, "self": 3, "twice": 3, "same": 3, "two": 3, "allow": 3, "via": 3, "either": [3, 7, 8], "parse_config_str": 3, "sconf": 3, "length": 3, "parse_nist_config": 3, "ip_fil": 3, "ioniz": 3, "energi": 3, "spectra": 3, "databas": 3, "other": [3, 5, 11, 12, 13, 14, 15, 17, 18], "similar": 3, "should": [3, 10], "element_nam": 3, "atomic_config": 3, "close": 3, "core": 3, "ne": 3, "3s2": 3, "_write_config": 3, "electronic_configur": 3, "molecul": [4, 9, 16], "ang2au": 4, "molecule_dir": 4, "angstrom": 4, "unit": [4, 10], "8897161646320724": 4, "load_molecul": 4, "add_atom": 4, "cart": 4, "cartesian": 4, "coordin": 4, "indent": 4, "charact": [4, 5], "_parse_molecules_xyz": 4, "xyz": 4, "file_nam": 4, "_parse_molecul": 4, "_write_load_molecul": 4, "mol": 4, "all": [4, 5, 7, 10], "helper": 5, "variou": 5, "find_fil": 5, "source_root": 5, "fals": [5, 7, 11, 12, 13, 14, 15, 17, 18], "find": 5, "root": 5, "relev": 5, "bool": [5, 7], "whether": [5, 11, 12, 13, 14, 15, 17, 18], "subdirectori": 5, "respect": [5, 11, 12, 13, 14, 15, 17, 18], "lookup_extens": 5, "bse": [5, 7], "identifi": [5, 7, 11, 12, 13, 18], "sanitize_basis_nam": 5, "sanit": 5, "exampl": [5, 7], "replac": 5, "word": 5, "hyphen": 5, "plu": 5, "sign": 5, "underscor": 5, "desanitize_basis_nam": 5, "desanit": 5, "_star": 5, "asterisk": 5, "unsanit": 5, "write_warn": 5, "fout": 5, "script_nam": 5, "prefix": 5, "warn": 5, "sai": 5, "wa": [5, 7, 11, 12, 13, 14, 15, 17, 18], "auto": 5, "overwritten": 5, "next": 5, "time": 5, "class": [5, 10, 11, 12, 13, 14, 15, 17, 18], "io": 5, "textiobas": 5, "open": 5, "addit": [5, 8, 11, 12, 13, 14, 15, 17, 18], "everi": 5, "comment": 5, "quick": 5, "dirti": 5, "fix": 5, "non": 5, "mai": [5, 10], "differ": [5, 7], "packag": 6, "refer": 6, "version": [6, 8, 10], "form": 6, "some": [6, 8], "download": [6, 7, 8], "those": 6, "generate_atomicinfo": 6, "generate_basi": [6, 8], "generate_dens": 6, "generate_elec_config": 6, "generate_molecul": 6, "helper_fxn": 6, "scrape_bs": [6, 8], "web": 7, "scraper": 7, "exchang": [7, 8], "o": 7, "outformat": 7, "optimize_gener": 7, "optim": 7, "bsebasissetscrap": 7, "base_url": 7, "http": [7, 10], "www": 7, "org": 7, "user_ag": 7, "nwchemex": [7, 10], "email": 7, "uncontract_gener": 7, "uncontract_seg": 7, "uncontract_spdf": 7, "make_gener": 7, "header_toggl": 7, "add_filt": 7, "metadata_kei": 7, "metadata": 7, "filter": 7, "updat": 7, "valid": [7, 11, 12, 13, 14, 15, 17, 18], "scrape": 7, "alreadi": 7, "append": 7, "match": 7, "exactli": 7, "when": [7, 10], "multipl": 7, "guarante": 7, "least": 7, "necessarili": 7, "howev": [7, 10], "appli": 7, "sequenti": 7, "so": 7, "famili": 7, "popl": [7, 8], "dun": [7, 8], "role": 7, "optri": 7, "onli": [7, 10, 11, 12, 13, 14, 15, 17, 18], "thei": [7, 10], "have": [7, 8, 11, 12, 13, 14, 15, 17, 18], "retriev": [7, 11, 12, 13, 14, 15, 17, 18], "member": 7, "filtered_basis_set": 7, "filtered_metadata": 7, "desir": [7, 8, 13], "download_basis_set": 7, "singl": 7, "left": 7, "empti": 7, "get": [7, 10], "comma": 7, "separ": 7, "could": 7, "obtain": 7, "clean": 7, "download_valid_basis_set": 7, "avail": [7, 11, 12, 13, 14, 15, 17, 18], "download_valid_format": 7, "get_extens": 7, "set_head": 7, "request": [7, 11, 12], "descript": 7, "who": 7, "ping": 7, "api": [7, 8, 10, 11, 12, 13, 14, 15, 17, 18], "send": 7, "share": 7, "set_default_format": 7, "get_filtered_basis_set": 7, "doe": [7, 11, 12, 13, 14, 15, 17, 18], "chang": 7, "validate_basis_set_nam": 7, "against": 7, "invalid": 7, "validate_format_nam": 7, "_create_param": 7, "_write_basis_set": 7, "basis_data": 7, "dot": 7, "page": [8, 10], "detail": [8, 10], "decis": 8, "most": [8, 10], "librari": [8, 10], "autogener": [8, 11, 12, 13, 14, 18], "util": 8, "data_manag": [8, 9], "reference_data": 8, "internet": 8, "Of": [8, 10], "particular": [8, 10, 13, 17], "note": [8, 10], "limit": 8, "select": 8, "master": 8, "branch": 8, "m": 8, "venv": 8, "data_venv": 8, "bin": 8, "activ": [8, 10], "pip": 8, "instal": [8, 9], "requir": 8, "onc": 8, "been": 8, "properli": [8, 10], "place": 8, "consult": 8, "document": [8, 10], "more": [8, 10], "usag": 8, "generated_data": 8, "larger": 8, "ahlrich": 8, "depend": 9, "design": 9, "submodul": 9, "modul": [9, 10], "nwx": [9, 16], "sto": [9, 16], "3g": [9, 16], "chemcach": 10, "cmaiz": 10, "its": 10, "build": 10, "system": 10, "instruct": 10, "usual": 10, "b": 10, "build_dir": 10, "dcmake_install_prefix": 10, "dcmake_toolchain_fil": 10, "toolchain": 10, "target": 10, "parallel": 10, "user": 10, "just": 10, "bbuild": 10, "you": 10, "want": 10, "your": 10, "sure": 10, "nwx_module_path": 10, "both": 10, "python_execut": 10, "python3_execut": 10, "interpret": 10, "virtual": 10, "environ": 10, "which": [10, 17], "python3": 10, "recogn": [10, 11, 12, 13, 14, 15, 17, 18], "build_test": 10, "truth": 10, "enabl": 10, "test": 10, "build_doc": 10, "variabl": 10, "defin": [10, 11, 12, 14, 15, 18], "nwx_cxx_api_doc": 10, "only_build_doc": 10, "process": 10, "skip": 10, "aspect": 10, "asid": 10, "chemcache_cxx_api": 10, "build_pybind11_pybind": 10, "On": 10, "built": 10, "enable_experimental_featur": 10, "function": [10, 11, 12, 13, 14, 15, 17, 18], "pre": 10, "releas": 10, "These": 10, "minimum": 10, "14": 10, "reli": 10, "17": 10, "standard": 10, "work": 10, "ani": 10, "compliant": 10, "gcc": 10, "9": 10, "x": 10, "newer": 10, "certain": 10, "featur": [10, 11, 12, 18], "home": 10, "build_pybind11_pythonbind": 10, "develop": 10, "section": [10, 11, 12, 13, 14, 15, 17, 18], "under": 10, "normal": 10, "circumst": 10, "ignor": 10, "primarili": 10, "complet": 10, "url": 10, "github": 10, "com": 10, "simul": 10, "see": 10, "inherit": 10, "catchorg": 10, "control": 10, "No": [11, 12, 13, 14, 15, 17, 18], "citat": [11, 12, 13, 14, 15, 17, 18], "satisfi": [11, 12, 13, 14, 15, 17, 18], "simd": [11, 12, 13, 14, 15, 17, 18], "unsign": [11, 12, 13, 14, 15, 17, 18], "long": [11, 12, 13, 14, 15, 17, 18], "accept": [11, 12, 13, 14, 15, 17, 18], "call": [11, 12, 13, 14, 15, 17, 18], "change_input": [11, 12, 13, 14, 15, 17, 18], "pass": [11, 12, 13, 14, 15, 17, 18], "summari": [11, 12, 13, 14, 15, 17, 18], "column": [11, 12, 13, 14, 15, 17, 18], "initi": [11, 12, 13, 14, 15, 17, 18], "A": [11, 12, 13, 14, 15, 17, 18], "human": [11, 12, 13, 14, 15, 17, 18], "what": [11, 12, 13, 14, 15, 17, 18], "subsect": [11, 12, 13, 14, 15, 17, 18], "head": [11, 12, 13, 14, 15, 17, 18], "within": [11, 12, 13, 14, 15, 17, 18], "inord": [11, 12, 13, 14, 15, 17, 18], "opaqu": [11, 12, 13, 14, 15, 17, 18], "influenc": [11, 12, 13, 14, 15, 17, 18], "memoiz": [11, 12, 13, 14, 15, 17, 18], "domain": [11, 12, 13, 14, 15, 17, 18], "restrict": [11, 12, 13, 14, 15, 17, 18], "criteria": [11, 12, 13, 14, 15, 17, 18], "obei": [11, 12, 13, 14, 15, 17, 18], "deem": [11, 12, 13, 14, 15, 17, 18], "tabul": [11, 12, 13, 14, 15, 17, 18], "how": [11, 12, 13, 14, 15, 17, 18], "comput": [11, 12, 13, 14, 15, 17, 18], "subset": [11, 12, 13, 14, 15, 17, 18], "advanc": [11, 12, 13, 14, 15, 17, 18], "instanc": [12, 13], "std": [12, 13, 14, 15], "pair": 12, "moleculefromstr": 13, "__cxx11": [13, 14, 15], "basic_str": [13, 14, 15], "char": [13, 14, 15], "char_trait": [13, 14, 15], "alloc": [13, 14, 15], "callback": [13, 17], "symbolfromz": 14, "zfromsymbol": 15, "determin": 17, "molecularbasisset": 17, "aobasisset": 17, "atomicbasisset": [17, 18], "contractedgaussian": [17, 18], "doubl": [17, 18]}, "objects": {"": [[6, 0, 0, "-", "data_management"]], "data_management": [[0, 0, 0, "-", "generate_atomicinfo"], [1, 0, 0, "-", "generate_basis"], [2, 0, 0, "-", "generate_densities"], [3, 0, 0, "-", "generate_elec_configs"], [4, 0, 0, "-", "generate_molecules"], [5, 0, 0, "-", "helper_fxns"], [7, 0, 0, "-", "scrape_bse"]], "data_management.generate_atomicinfo": [[0, 1, 1, "", "AtomicData"], [0, 3, 1, "", "_parse_ciaaw_isotopes"], [0, 3, 1, "", "_parse_ciaww_mass"], [0, 3, 1, "", "_write_atoms_average"], [0, 3, 1, "", "_write_atoms_isotope"], [0, 3, 1, "", "_write_sym_from_z"], [0, 3, 1, "", "_write_z_from_sym"], [0, 3, 1, "", "main"], [0, 3, 1, "", "parse_args"], [0, 3, 1, "", "parse_symbols"]], "data_management.generate_atomicinfo.AtomicData": [[0, 2, 1, "", "__repr__"], [0, 2, 1, "", "add_isotope"]], "data_management.generate_basis": [[1, 1, 1, "", "Shell"], [1, 3, 1, "", "_parse_bases"], [1, 3, 1, "", "_parse_bases_gbs"], [1, 3, 1, "", "_parse_bases_nw"], [1, 3, 1, "", "_write_bases"], [1, 3, 1, "", "_write_basis_files"], [1, 3, 1, "", "main"], [1, 3, 1, "", "parse_args"]], "data_management.generate_basis.Shell": [[1, 2, 1, "", "add_prim"], [1, 2, 1, "", "cxxify"]], "data_management.generate_densities": [[2, 3, 1, "", "_parse_densities"], [2, 3, 1, "", "_parse_densities_dat"], [2, 3, 1, "", "_parse_densities_xml"], [2, 3, 1, "", "_write_den_files"], [2, 3, 1, "", "_write_densities"], [2, 3, 1, "", "main"], [2, 3, 1, "", "make_square_arr"], [2, 3, 1, "", "parse_args"]], "data_management.generate_elec_configs": [[3, 1, 1, "", "AtomicData"], [3, 5, 1, "", "LMAX"], [3, 5, 1, "", "NMAX"], [3, 3, 1, "", "_write_configs"], [3, 5, 1, "", "i_to_lchar"], [3, 5, 1, "", "lchar_to_i"], [3, 3, 1, "", "main"], [3, 3, 1, "", "parse_args"], [3, 3, 1, "", "parse_config_str"], [3, 3, 1, "", "parse_nist_configs"], [3, 3, 1, "", "parse_symbols"]], "data_management.generate_elec_configs.AtomicData": [[3, 2, 1, "", "__repr__"], [3, 4, 1, "", "config"], [3, 4, 1, "", "config_full"]], "data_management.generate_molecules": [[4, 1, 1, "", "Molecule"], [4, 3, 1, "", "_parse_molecules"], [4, 3, 1, "", "_parse_molecules_xyz"], [4, 3, 1, "", "_write_load_molecules"], [4, 3, 1, "", "main"], [4, 3, 1, "", "parse_args"]], "data_management.generate_molecules.Molecule": [[4, 2, 1, "", "__repr__"], [4, 2, 1, "", "add_atom"], [4, 2, 1, "", "cxxify"]], "data_management.helper_fxns": [[5, 3, 1, "", "desanitize_basis_name"], [5, 3, 1, "", "find_files"], [5, 3, 1, "", "lookup_extension"], [5, 3, 1, "", "sanitize_basis_name"], [5, 3, 1, "", "write_warning"]], "data_management.scrape_bse": [[7, 1, 1, "", "BSEBasisSetScraper"], [7, 3, 1, "", "_write_basis_set"], [7, 3, 1, "", "main"], [7, 3, 1, "", "parse_args"]], "data_management.scrape_bse.BSEBasisSetScraper": [[7, 2, 1, "", "_create_params"], [7, 2, 1, "", "add_filter"], [7, 2, 1, "", "download_basis_set"], [7, 2, 1, "", "download_valid_basis_sets"], [7, 2, 1, "", "download_valid_formats"], [7, 2, 1, "", "get_extension"], [7, 2, 1, "", "get_filtered_basis_sets"], [7, 2, 1, "", "set_default_format"], [7, 2, 1, "", "set_header"], [7, 2, 1, "", "validate_basis_set_name"], [7, 2, 1, "", "validate_format_name"]]}, "objtypes": {"0": "py:module", "1": "py:class", "2": "py:method", "3": "py:function", "4": "py:property", "5": "py:data"}, "objnames": {"0": ["py", "module", "Python module"], "1": ["py", "class", "Python class"], "2": ["py", "method", "Python method"], "3": ["py", "function", "Python function"], "4": ["py", "property", "Python property"], "5": ["py", "data", "Python data"]}, "titleterms": {"data_manag": [0, 1, 2, 3, 4, 5, 6, 7], "generate_atomicinfo": 0, "usag": [0, 1, 2, 3, 4, 7], "modul": [0, 1, 2, 3, 4, 5, 7, 11, 12, 13, 14, 15, 16, 17, 18], "content": [0, 1, 2, 3, 4, 5, 7, 9], "class": [0, 1, 3, 4, 7], "function": [0, 1, 2, 3, 4, 5, 7], "generate_basi": 1, "generate_dens": 2, "generate_elec_config": 3, "attribut": 3, "generate_molecul": 4, "helper_fxn": 5, "submodul": [6, 11, 12, 13, 14, 15, 17, 18], "scrape_bs": 7, "chemcach": [8, 9], "design": 8, "gener": 8, "sourc": 8, "file": 8, "api": [9, 16], "instal": 10, "command": 10, "configur": 10, "option": 10, "depend": 10, "requir": 10, "cmake": 10, "c": 10, "compil": 10, "doxygen": 10, "python": 10, "other": 10, "simd": 10, "catch2": 10, "atom": [11, 12, 14, 18], "abund": 11, "weight": 11, "mass": [11, 12], "pleas": [11, 12, 13, 14, 15, 17, 18], "cite": [11, 12, 13, 14, 15, 17, 18], "properti": [11, 12, 13, 14, 15, 17, 18], "type": [11, 12, 13, 14, 15, 17, 18], "input": [11, 12, 13, 14, 15, 17, 18], "quick": [11, 12, 13, 14, 15, 17, 18], "refer": [11, 12, 13, 14, 15, 17, 18], "detail": [11, 12, 13, 14, 15, 17, 18], "descript": [11, 12, 13, 14, 15, 17, 18], "id": [11, 12, 18], "result": [11, 12, 13, 14, 15, 17, 18], "isotop": 12, "nwx": 13, "molecul": [13, 17], "nwchemex": 13, "string": 13, "symbol": [14, 15], "from": [14, 15], "z": [14, 15], "number": 14, "sto": [17, 18], "3g": [17, 18], "molecular": 17, "basi": [17, 18], "set": [17, 18]}, "envversion": {"sphinx.domains.c": 3, "sphinx.domains.changeset": 1, "sphinx.domains.citation": 1, "sphinx.domains.cpp": 9, "sphinx.domains.index": 1, "sphinx.domains.javascript": 3, "sphinx.domains.math": 2, "sphinx.domains.python": 4, "sphinx.domains.rst": 2, "sphinx.domains.std": 2, "sphinx.ext.todo": 2, "sphinx": 60}, "alltitles": {"data_management.generate_atomicinfo": [[0, "module-data_management.generate_atomicinfo"]], "Usage": [[0, "usage"], [1, "usage"], [2, "usage"], [3, "usage"], [4, "usage"], [7, "usage"]], "Module Contents": [[0, "module-contents"], [1, "module-contents"], [2, "module-contents"], [3, "module-contents"], [4, "module-contents"], [5, "module-contents"], [7, "module-contents"]], "Classes": [[0, "classes"], [1, "classes"], [3, "classes"], [4, "classes"], [7, "classes"]], "Functions": [[0, "functions"], [1, "functions"], [2, "functions"], [3, "functions"], [4, "functions"], [5, "functions"], [7, "functions"]], "data_management.generate_basis": [[1, "module-data_management.generate_basis"]], "data_management.generate_densities": [[2, "module-data_management.generate_densities"]], "data_management.generate_elec_configs": [[3, "module-data_management.generate_elec_configs"]], "Attributes": [[3, "attributes"]], "data_management.generate_molecules": [[4, "module-data_management.generate_molecules"]], "data_management.helper_fxns": [[5, "module-data_management.helper_fxns"]], "data_management": [[6, "module-data_management"]], "Submodules": [[6, "submodules"], [11, "submodules"], [12, "submodules"], [13, "submodules"], [14, "submodules"], [15, "submodules"], [17, "submodules"], [18, "submodules"]], "data_management.scrape_bse": [[7, "module-data_management.scrape_bse"]], "ChemCache Design": [[8, "chemcache-design"]], "Generation of Source Files": [[8, "generation-of-source-files"]], "ChemCache": [[9, "chemcache"]], "Contents:": [[9, null]], "APIs:": [[9, null]], "Installation": [[10, "installation"]], "Install Command": [[10, "install-command"]], "Configuration Options": [[10, "configuration-options"]], "Dependencies": [[10, "dependencies"]], "Required Dependencies": [[10, "required-dependencies"]], "CMake": [[10, "cmake"]], "C++ Compiler": [[10, "c-compiler"]], "Optional Dependencies": [[10, "optional-dependencies"]], "Doxygen": [[10, "doxygen"]], "Python": [[10, "python"]], "Other Dependencies": [[10, "other-dependencies"]], "SimDE": [[10, "simde"]], "Catch2": [[10, "catch2"]], "Atom": [[11, "atom"]], "Atoms with Abundance-Weighted Mass": [[11, "atoms-with-abundance-weighted-mass"]], "Please Cite": [[11, "please-cite"], [12, "please-cite"], [13, "please-cite"], [14, "please-cite"], [15, "please-cite"], [17, "please-cite"], [18, "please-cite"]], "Property Types": [[11, "property-types"], [12, "property-types"], [13, "property-types"], [14, "property-types"], [15, "property-types"], [17, "property-types"], [18, "property-types"]], "Module Inputs": [[11, "module-inputs"], [12, "module-inputs"], [13, "module-inputs"], [14, "module-inputs"], [15, "module-inputs"], [17, "module-inputs"], [18, "module-inputs"]], "Quick Reference": [[11, "quick-reference"], [12, "quick-reference"], [13, "quick-reference"], [14, "quick-reference"], [15, "quick-reference"], [17, "quick-reference"], [18, "quick-reference"]], "Detailed Descriptions": [[11, "detailed-descriptions"], [12, "detailed-descriptions"], [13, "detailed-descriptions"], [14, "detailed-descriptions"], [15, "detailed-descriptions"], [17, "detailed-descriptions"], [18, "detailed-descriptions"]], "Atom ID": [[11, "atom-id"], [12, "atom-id"], [18, "atom-id"]], "Module Results": [[11, "module-results"], [12, "module-results"], [13, "module-results"], [14, "module-results"], [15, "module-results"], [17, "module-results"], [18, "module-results"]], "Atom Isotope": [[12, "atom-isotope"]], "Atoms with Isotope Mass": [[12, "atoms-with-isotope-mass"]], "NWX Molecules": [[13, "nwx-molecules"]], "NWChemEx Molecules": [[13, "nwchemex-molecules"]], "String": [[13, "string"]], "Symbol from Z": [[14, "symbol-from-z"]], "Atomic Symbol from Atomic Number": [[14, "atomic-symbol-from-atomic-number"]], "Z": [[14, "z"]], "Z from Symbol": [[15, "z-from-symbol"]], "Symbol": [[15, "symbol"]], "Modules API": [[16, "modules-api"]], "sto-3g": [[17, "sto-3g"]], "Molecular Basis Set": [[17, "molecular-basis-set"]], "Molecule": [[17, "molecule"]], "sto-3g atomic basis": [[18, "sto-3g-atomic-basis"]], "sto-3g atomic basis set": [[18, "sto-3g-atomic-basis-set"]]}, "indexentries": {"atomicdata (class in data_management.generate_atomicinfo)": [[0, "data_management.generate_atomicinfo.AtomicData"]], "__repr__() (data_management.generate_atomicinfo.atomicdata method)": [[0, "data_management.generate_atomicinfo.AtomicData.__repr__"]], "_parse_ciaaw_isotopes() (in module data_management.generate_atomicinfo)": [[0, "data_management.generate_atomicinfo._parse_ciaaw_isotopes"]], "_parse_ciaww_mass() (in module data_management.generate_atomicinfo)": [[0, "data_management.generate_atomicinfo._parse_ciaww_mass"]], "_write_atoms_average() (in module data_management.generate_atomicinfo)": [[0, "data_management.generate_atomicinfo._write_atoms_average"]], "_write_atoms_isotope() (in module data_management.generate_atomicinfo)": [[0, "data_management.generate_atomicinfo._write_atoms_isotope"]], "_write_sym_from_z() (in module data_management.generate_atomicinfo)": [[0, "data_management.generate_atomicinfo._write_sym_from_z"]], "_write_z_from_sym() (in module data_management.generate_atomicinfo)": [[0, "data_management.generate_atomicinfo._write_z_from_sym"]], "add_isotope() (data_management.generate_atomicinfo.atomicdata method)": [[0, "data_management.generate_atomicinfo.AtomicData.add_isotope"]], "data_management.generate_atomicinfo": [[0, "module-data_management.generate_atomicinfo"]], "main() (in module data_management.generate_atomicinfo)": [[0, "data_management.generate_atomicinfo.main"]], "module": [[0, "module-data_management.generate_atomicinfo"], [1, "module-data_management.generate_basis"], [2, "module-data_management.generate_densities"], [3, "module-data_management.generate_elec_configs"], [4, "module-data_management.generate_molecules"], [5, "module-data_management.helper_fxns"], [6, "module-data_management"], [7, "module-data_management.scrape_bse"]], "parse_args() (in module data_management.generate_atomicinfo)": [[0, "data_management.generate_atomicinfo.parse_args"]], "parse_symbols() (in module data_management.generate_atomicinfo)": [[0, "data_management.generate_atomicinfo.parse_symbols"]], "shell (class in data_management.generate_basis)": [[1, "data_management.generate_basis.Shell"]], "_parse_bases() (in module data_management.generate_basis)": [[1, "data_management.generate_basis._parse_bases"]], "_parse_bases_gbs() (in module data_management.generate_basis)": [[1, "data_management.generate_basis._parse_bases_gbs"]], "_parse_bases_nw() (in module data_management.generate_basis)": [[1, "data_management.generate_basis._parse_bases_nw"]], "_write_bases() (in module data_management.generate_basis)": [[1, "data_management.generate_basis._write_bases"]], "_write_basis_files() (in module data_management.generate_basis)": [[1, "data_management.generate_basis._write_basis_files"]], "add_prim() (data_management.generate_basis.shell method)": [[1, "data_management.generate_basis.Shell.add_prim"]], "cxxify() (data_management.generate_basis.shell method)": [[1, "data_management.generate_basis.Shell.cxxify"]], "data_management.generate_basis": [[1, "module-data_management.generate_basis"]], "main() (in module data_management.generate_basis)": [[1, "data_management.generate_basis.main"]], "parse_args() (in module data_management.generate_basis)": [[1, "data_management.generate_basis.parse_args"]], "_parse_densities() (in module data_management.generate_densities)": [[2, "data_management.generate_densities._parse_densities"]], "_parse_densities_dat() (in module data_management.generate_densities)": [[2, "data_management.generate_densities._parse_densities_dat"]], "_parse_densities_xml() (in module data_management.generate_densities)": [[2, "data_management.generate_densities._parse_densities_xml"]], "_write_den_files() (in module data_management.generate_densities)": [[2, "data_management.generate_densities._write_den_files"]], "_write_densities() (in module data_management.generate_densities)": [[2, "data_management.generate_densities._write_densities"]], "data_management.generate_densities": [[2, "module-data_management.generate_densities"]], "main() (in module data_management.generate_densities)": [[2, "data_management.generate_densities.main"]], "make_square_arr() (in module data_management.generate_densities)": [[2, "data_management.generate_densities.make_square_arr"]], "parse_args() (in module data_management.generate_densities)": [[2, "data_management.generate_densities.parse_args"]], "atomicdata (class in data_management.generate_elec_configs)": [[3, "data_management.generate_elec_configs.AtomicData"]], "lmax (in module data_management.generate_elec_configs)": [[3, "data_management.generate_elec_configs.LMAX"]], "nmax (in module data_management.generate_elec_configs)": [[3, "data_management.generate_elec_configs.NMAX"]], "__repr__() (data_management.generate_elec_configs.atomicdata method)": [[3, "data_management.generate_elec_configs.AtomicData.__repr__"]], "_write_configs() (in module data_management.generate_elec_configs)": [[3, "data_management.generate_elec_configs._write_configs"]], "config (data_management.generate_elec_configs.atomicdata property)": [[3, "data_management.generate_elec_configs.AtomicData.config"]], "config_full (data_management.generate_elec_configs.atomicdata property)": [[3, "data_management.generate_elec_configs.AtomicData.config_full"]], "data_management.generate_elec_configs": [[3, "module-data_management.generate_elec_configs"]], "i_to_lchar (in module data_management.generate_elec_configs)": [[3, "data_management.generate_elec_configs.i_to_lchar"]], "lchar_to_i (in module data_management.generate_elec_configs)": [[3, "data_management.generate_elec_configs.lchar_to_i"]], "main() (in module data_management.generate_elec_configs)": [[3, "data_management.generate_elec_configs.main"]], "parse_args() (in module data_management.generate_elec_configs)": [[3, "data_management.generate_elec_configs.parse_args"]], "parse_config_str() (in module data_management.generate_elec_configs)": [[3, "data_management.generate_elec_configs.parse_config_str"]], "parse_nist_configs() (in module data_management.generate_elec_configs)": [[3, "data_management.generate_elec_configs.parse_nist_configs"]], "parse_symbols() (in module data_management.generate_elec_configs)": [[3, "data_management.generate_elec_configs.parse_symbols"]], "molecule (class in data_management.generate_molecules)": [[4, "data_management.generate_molecules.Molecule"]], "__repr__() (data_management.generate_molecules.molecule method)": [[4, "data_management.generate_molecules.Molecule.__repr__"]], "_parse_molecules() (in module data_management.generate_molecules)": [[4, "data_management.generate_molecules._parse_molecules"]], "_parse_molecules_xyz() (in module data_management.generate_molecules)": [[4, "data_management.generate_molecules._parse_molecules_xyz"]], "_write_load_molecules() (in module data_management.generate_molecules)": [[4, "data_management.generate_molecules._write_load_molecules"]], "add_atom() (data_management.generate_molecules.molecule method)": [[4, "data_management.generate_molecules.Molecule.add_atom"]], "cxxify() (data_management.generate_molecules.molecule method)": [[4, "data_management.generate_molecules.Molecule.cxxify"]], "data_management.generate_molecules": [[4, "module-data_management.generate_molecules"]], "main() (in module data_management.generate_molecules)": [[4, "data_management.generate_molecules.main"]], "parse_args() (in module data_management.generate_molecules)": [[4, "data_management.generate_molecules.parse_args"]], "data_management.helper_fxns": [[5, "module-data_management.helper_fxns"]], "desanitize_basis_name() (in module data_management.helper_fxns)": [[5, "data_management.helper_fxns.desanitize_basis_name"]], "find_files() (in module data_management.helper_fxns)": [[5, "data_management.helper_fxns.find_files"]], "lookup_extension() (in module data_management.helper_fxns)": [[5, "data_management.helper_fxns.lookup_extension"]], "sanitize_basis_name() (in module data_management.helper_fxns)": [[5, "data_management.helper_fxns.sanitize_basis_name"]], "write_warning() (in module data_management.helper_fxns)": [[5, "data_management.helper_fxns.write_warning"]], "data_management": [[6, "module-data_management"]], "bsebasissetscraper (class in data_management.scrape_bse)": [[7, "data_management.scrape_bse.BSEBasisSetScraper"]], "_create_params() (data_management.scrape_bse.bsebasissetscraper method)": [[7, "data_management.scrape_bse.BSEBasisSetScraper._create_params"]], "_write_basis_set() (in module data_management.scrape_bse)": [[7, "data_management.scrape_bse._write_basis_set"]], "add_filter() (data_management.scrape_bse.bsebasissetscraper method)": [[7, "data_management.scrape_bse.BSEBasisSetScraper.add_filter"]], "data_management.scrape_bse": [[7, "module-data_management.scrape_bse"]], "download_basis_set() (data_management.scrape_bse.bsebasissetscraper method)": [[7, "data_management.scrape_bse.BSEBasisSetScraper.download_basis_set"]], "download_valid_basis_sets() (data_management.scrape_bse.bsebasissetscraper method)": [[7, "data_management.scrape_bse.BSEBasisSetScraper.download_valid_basis_sets"]], "download_valid_formats() (data_management.scrape_bse.bsebasissetscraper method)": [[7, "data_management.scrape_bse.BSEBasisSetScraper.download_valid_formats"]], "get_extension() (data_management.scrape_bse.bsebasissetscraper method)": [[7, "data_management.scrape_bse.BSEBasisSetScraper.get_extension"]], "get_filtered_basis_sets() (data_management.scrape_bse.bsebasissetscraper method)": [[7, "data_management.scrape_bse.BSEBasisSetScraper.get_filtered_basis_sets"]], "main() (in module data_management.scrape_bse)": [[7, "data_management.scrape_bse.main"]], "parse_args() (in module data_management.scrape_bse)": [[7, "data_management.scrape_bse.parse_args"]], "set_default_format() (data_management.scrape_bse.bsebasissetscraper method)": [[7, "data_management.scrape_bse.BSEBasisSetScraper.set_default_format"]], "set_header() (data_management.scrape_bse.bsebasissetscraper method)": [[7, "data_management.scrape_bse.BSEBasisSetScraper.set_header"]], "validate_basis_set_name() (data_management.scrape_bse.bsebasissetscraper method)": [[7, "data_management.scrape_bse.BSEBasisSetScraper.validate_basis_set_name"]], "validate_format_name() (data_management.scrape_bse.bsebasissetscraper method)": [[7, "data_management.scrape_bse.BSEBasisSetScraper.validate_format_name"]]}})