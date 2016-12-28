from . import t_utils as tu
from .t_predic import t_predic


class TTideCon(dict):
    """The ttide consituents class.
    """

    def t_predic(self, time):
        return t_predic(time,
                        names=self['nameu'], freq=self['fu'],
                        tidecon=self['tidecon'], lat=self['lat'],
                        ltype=self['ltype'], synth=self['synth'])

    __call__ = t_predic

    def pandas_style(self, to_file=None):
        outstr = tu.pandas_style(self)
        if to_file is None:
            return outstr
        elif isinstance(to_file, file):
            to_file.write(outstr)
        else:
            with open(to_file, 'w') as fl:
                fl.write(outstr)

    def classic_style(self, to_file=None):
        outstr = tu.classic_style(self)
        if to_file is None:
            return outstr
        elif isinstance(to_file, file):
            to_file.write(outstr)
        else:
            with open(to_file, 'w') as fl:
                fl.write(outstr)
