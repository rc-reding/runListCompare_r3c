#!/usr/bin/env python
# Run from pipeline root directory

import os
import shutil
import hashlib
import subprocess


cwd = os.getcwd()

def run(cmd, cwd=cwd):
    return subprocess.run(cmd,
                          cwd=cwd,
                          shell=True,
                          check=True,
                          universal_newlines=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)

def test_all_cd():
    # PhyML, seq_keep=0, varsite_keep=0
    out_path = 'tests/output/cd'
    if os.path.isdir(out_path):
        shutil.rmtree(out_path)
    else:
        os.makedirs(out_path)
    run('python runListCompare.py tests/data/cd/cd.ini')
    with open('tests/output/cd/align_positions.txt', 'rb') as fh:
        contents = fh.read()
    assert hashlib.md5(contents).hexdigest() == '83da02b151226c443430c5c54c43ccea'


def test_all_ec():
    # IQtree, seq_keep=0.7, varsite_keep=0.7, trailing empty lines in seqlist
    out_path = 'tests/output/ec'
    if os.path.isdir(out_path):
        shutil.rmtree(out_path)
    else:
        os.makedirs(out_path)
    run('python runListCompare.py tests/data/ec/ec.ini')
    with open('tests/output/ec/align_positions.txt', 'rb') as fh:
        contents = fh.read()
    assert hashlib.md5(contents).hexdigest() == '811f7f54d6d5a5867bdeb718f1ac6e9a'
