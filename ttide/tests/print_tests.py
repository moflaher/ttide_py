from ttide.t_tide import t_tide
import ttide.tests.base as bmod
from io import StringIO
import sys
import numpy as np


def compare_string2file(string, fname):
    with open(bmod.testdir + 'data/print/' + fname, 'r') as fl:
        tdata = fl.read().replace('\r\n', '\n')
    string = string.replace('\r\n', '\n')
    # The above .replace() calls makes sure that line-endings are
    # consistent
    assert string == tdata, ("Test failed on file '%s'" % fname)


def gen_print_tests(make_data=False):
    """A generator function for the 'print' tests.

    If make_data=True, this will write out the data files.
    """
    out_old = sys.stdout
    for kwargs, fname in bmod.cases:
        np.random.seed(29034230)
        # Redirect stdout to a StringIO instance
        sys.stdout = stdout = StringIO()
        t_tide(**kwargs)
        sys.stdout = out_old
        if make_data:
            with open(bmod.testdir + 'data/print/' + fname, 'w') as fl:
                fl.write(stdout.getvalue())
            yield None
        else:
            yield compare_string2file, stdout.getvalue(), fname


# if __name__ == '__main__':

#     # ###
#     # # This block generates the output files.
#     # # USE WITH CAUTION!
#     # for t in gen_print_tests(make_data=True):
#     #     pass

#     ###
#     # This block is for joining files for comparison to the old testdata.out file.
#     # It is not identical, but nearly so...
#     outstr = ''
#     for kwargs, fl in cases:
#         with open(testdir + 'data/' + fl, 'r') as fl:
#             outstr += fl.read()
#             outstr += '\n' + '*' * 80 + '\n\n'
#     with open(testdir + 'data/testdata.out2', 'w') as fl:
#         fl.write(outstr)
