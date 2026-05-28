import sys
from io import StringIO

import numpy as np

import base as bmod
from ttide.t_tide import t_tide


def gen_print_data():
    out_old = sys.stdout

    for kwargs, fname in bmod.cases:
        np.random.seed(29034230)

        sys.stdout = stdout = StringIO()
        t_tide(**kwargs)
        sys.stdout = out_old

        with open(bmod.testdir + 'data/print/' + fname, 'w') as fl:
            fl.write(stdout.getvalue())


if __name__ == '__main__':
    gen_print_data()