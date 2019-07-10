from ttide.tests import base
import ttide as tt
import logging
logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def test_errcalc():
    """
    test different values for errcalc
    should not fail
    """
    errcalc_opts = ("cboot", "linear", "none")


    con_cboot = None

    msg_format = """

    {} should be the same for all errcalc values.
    But:
    errcalc={}, yields mean {}={}
    and
    errcalc={}, yields mean {}={}
    """


    for errcalc in errcalc_opts:
        con = tt.t_tide(base.ein, dt=1., synth=0, ray=0.5, out_style=None, errcalc=errcalc)

        print(con["tidecon"])

        if con_cboot is None:
            con_cboot = con
            amp_mean_cboot = con_cboot["tidecon"][:, 0].mean()
            pha_mean_cboot = con_cboot["tidecon"][:, 2].mean()
        else:
            # amplitude and phase should not change
            amp_mean = con["tidecon"][:, 0].mean()
            pha_mean = con["tidecon"][:, 2].mean()




            param = "amplitude"
            msg = msg_format.format(param, errcalc, param, amp_mean, errcalc_opts[0], param, amp_mean_cboot)
            assert amp_mean == amp_mean_cboot, msg

            param = "phase"
            msg = msg_format.format(param, errcalc, param, pha_mean, errcalc_opts[0], param, pha_mean_cboot)
            assert pha_mean == pha_mean_cboot, msg


def test_errcalc_same():
    """
    Make sure constituents are the same for different errcalc_opts
    """

    errcalc_opts = ("cboot", "linear", "none")

    for errcalc in errcalc_opts:
        tc = tt.t_tide(base.ein, dt=1., synth=0, ray=0.5, errcalc=errcalc)

        logger.debug(tc)

        raise Exception


def test_all():
    test_errcalc_same()
    test_errcalc()


if __name__ == '__main__':
    test_all()
