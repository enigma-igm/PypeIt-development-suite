"""
Check the background-exposure grouping used by jwst_run_redux for nod subtraction.

Usage: python3 -m pypeitdev.jwst.test_bkg_indices
"""
import os

from pypeitdev.jwst.jwst_targets import jwst_targets


def bkg_indices_from(basenames):
    # Same logic as in jwst_run_redux
    nexp = len(basenames)
    exp_sets = [b.rsplit('_', 1)[0] for b in basenames]
    return [tuple(j for j in range(nexp) if exp_sets[j] == exp_sets[i] and j != i) for i in range(nexp)]


def basenames_from(exp_list):
    return [os.path.basename(f).split('_nrs2')[0] for f in exp_list[1]]


def main():
    # One set of three must reproduce the old hardcoded list
    one_set = ['jw04233003001_03101_00002', 'jw04233003001_03101_00003', 'jw04233003001_03101_00004']
    assert bkg_indices_from(one_set) == [(1, 2), (0, 2), (0, 1)]

    # Two sets of three for program 4713 source 6733
    exp_list, _, _ = jwst_targets('4713', 'G140M', '4713_6733')
    names = basenames_from(exp_list)
    bkg = bkg_indices_from(names)
    for name, b in zip(names, bkg):
        print(name, '-> background', [names[j] for j in b])
    assert bkg == [(1, 2), (0, 2), (0, 1), (4, 5), (3, 5), (3, 4)]
    print('OK')


if __name__ == '__main__':
    main()
