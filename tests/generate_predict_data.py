import copy

import numpy as np

import base as bmod
from ttide.t_tide import t_tide


cases = copy.deepcopy(bmod.cases)
t = np.arange(0, 60, 0.3)


def gen_predict_data():
    for kwargs, fname in cases:
        kwargs['out_style'] = None
        out = t_tide(**kwargs)
        xout = out.t_predic(t)

        np.savetxt(
            bmod.testdir + 'data/predict/' + fname,
            xout.view(float),
            fmt='%0.5f',
        )


if __name__ == '__main__':
    gen_predict_data()