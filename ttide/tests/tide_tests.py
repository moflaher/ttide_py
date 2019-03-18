from ttide.tests import base
import ttide as tt


def test_errcalc():
    """
    test different values for errcalc
    should not fail
    """
    errcalc_opts = ("cboot", "linear", "none")

    for errcalc in errcalc_opts:
        tt.t_tide(base.ein, dt=1., synth=0, ray=0.5, out_style=None, errcalc=errcalc)


def test_all():
    test_errcalc()


if __name__ == '__main__':
    test_all()
