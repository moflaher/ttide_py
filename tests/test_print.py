import copy
import sys
from io import StringIO

import numpy as np
import pytest

import base as bmod
from ttide.t_tide import t_tide


def compare_string2file(string, fname):
    with open(bmod.testdir + 'data/print/' + fname, 'r') as fl:
        tdata = fl.read().replace('\r\n', '\n')
    string = string.replace('\r\n', '\n')
    # The above .replace() calls makes sure that line-endings are
    # consistent
    assert string == tdata, ("Test failed on file '%s'" % fname)

@pytest.mark.parametrize(
    "kwargs,fname",
    bmod.cases,
    ids=[fname for _, fname in bmod.cases],
)
def test_print_cases(kwargs, fname):
    kwargs = copy.deepcopy(kwargs)

    np.random.seed(29034230)

    out_old = sys.stdout
    stdout = StringIO()

    try:
        sys.stdout = stdout
        t_tide(**kwargs)
    finally:
        sys.stdout = out_old

    compare_string2file(stdout.getvalue(), fname)

