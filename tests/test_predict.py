import numpy as np
from ttide.t_tide import t_tide
import base as bmod
import copy
import pytest

cases = copy.deepcopy(bmod.cases)

t = np.arange(0, 60, 0.3)


def compare_vec2file(x0, fname):
    x1 = np.loadtxt(bmod.testdir + 'data/predict/' + fname)
    if len(x1) == 2 * len(x0):
        x1 = x1.view(complex)
    assert (np.abs(x0 - x1) < 1e-2).all(), ("Test failed on file '%s'" % fname)

@pytest.mark.parametrize(
    "kwargs,fname",
    cases,
    ids=[fname for _, fname in cases],
)
def test_predict_cases(kwargs, fname):
    kwargs = copy.deepcopy(kwargs)
    kwargs['out_style'] = None

    out = t_tide(**kwargs)
    xout = out.t_predic(t)

    compare_vec2file(xout, fname)

# if __name__ == '__main__':

#     ###
#     # This block generates the output files.
#     # USE WITH CAUTION!
#     for tst in gen_predict_tests(make_data=True):
#         pass

#     # for f, vec, fname in gen_predict_tests():
#     #     f(vec, fname)
