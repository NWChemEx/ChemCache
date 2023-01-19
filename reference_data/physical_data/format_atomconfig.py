#!/usr/bin/env python3
# Author Kevin Gasperich
import re
import sys

# current file generated with:
# format_atomconfig.py ElementNames.txt NIST-ATOMICION.txt > NIST-ATOMCONFIGS.txt

LMAX = 3
i_to_lchar = 'spdfgh'[:LMAX+1]
lchar_to_i = {c: i for i, c in enumerate(i_to_lchar)}
BRACKETS = '{}'


def get_symbols(fpath):
    """
    input:
      fpath, str: path to file containing one line per element
                  elements should be in order, starting from 0
                  each line should contain:
                      atomic_number  element_symbol  element_name

    return:
      ztos, list: map from atomic number to tuple with (symbol, name)
      stoz, dict: map from symbol to atomic number
    """
    ztos = []
    stoz = {}
    with open(fpath, 'r') as f:
        for line in f:
            z = int(line.split()[0])
            assert(z == len(ztos))
            ztos.append(tuple(line.split()[1:]))
            stoz[ztos[-1][0]] = z
    return ztos, stoz


def getconfigs(conf_path, elem_path):
    """
    return:
      header: list of strings (header lines from input file)
              (input assumed to have 5 lines before data, no lines after)
      ztoconf, dict: map from atomic number to tuple representing 
                     grd. state atomic config, separated by l
                     i.e. (N_s, N_p, N_d, N_f)
    """
    header = []
    ztoconf = {0: (0, 0, 0, 0)}
    ztos, stoz = get_symbols(elem_path)
    n_header_lines = 5
    with open(conf_path, 'r') as f:
        for _, line in zip(range(n_header_lines), f):
            header.append(line)
        for line in f:
            z_s, name, conf_s, _ = line.split()
            z = int(z_s)  # atomic number
            conf_s = conf_s.strip('[').split(']')  # separate core if present
            # get core config
            ccore = ztoconf[stoz[conf_s[0]]] if len(
                conf_s) == 2 else (0,)*(LMAX+1)
            # get string representing remaining elec config
            cval_s = conf_s[-1]
            revcval_s = cval_s[::-1]  # reverse config (easier to parse)
            shells = re.findall(r'(\d*[a-z]\d)', revcval_s)
            cval = [0]*(LMAX+1)  # valence config
            for shellrev in shells:
                shell = shellrev[::-1]  # put back in correct order
                matches = re.match(r'(?P<n>\d)(?P<l>\w)(?P<e>\d*)', shell)
                # 1-elec shells have implicit 'e'
                nelec = int(matches['e']) if matches['e'] else 1
                cval[lchar_to_i[matches['l']]] += nelec
            # add sum contributions and add to dict
            ztoconf[z] = tuple(sum(nc) for nc in zip(cval, ccore))
    return header, ztoconf


if __name__ == "__main__":
    assert len(sys.argv) == 3, "Usage: format_atomconfig.py elementfile nistfile"

    elemfile = sys.argv[1]
    nistfile = sys.argv[2]

    header, ztoconf = getconfigs(nistfile, elemfile)
    ztos, stoz = get_symbols(elemfile)
    print(''.join(header), end='')
    for z, conf in ztoconf.items():
        print(f'{z:<3d}  {ztos[z][0]:<2s}   {ztos[z][1]:<13s}   {conf}'.
              replace('(', BRACKETS[0]).
              replace(')', BRACKETS[1]))
