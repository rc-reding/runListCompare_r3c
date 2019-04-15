#!/usr/bin/env python
# Run from pipeline root directory

import subprocess

def test_all():
	subprocess.check_output('python runListCompare.py tests/data/test.ini', shell=True)
