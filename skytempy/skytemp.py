# coding: utf-8

from __future__ import print_function, division

import numpy as np, argparse as ap
from scipy.interpolate import CubicSpline
from astropy.table import Table
import astropy.units as u
from astropy_healpix import HEALPix


__all__ = ['SkyTemp',]


def proc_args():
    pars = ap.ArgumentParser(description="Sky temperature calculator")
    pars.add_argument('fits', help="Path to fits file containing temperature"
                      "map")
    pars.add_argument('gl', type=float, help="Galactic Longitude in degrees")
    pars.add_argument('gb', type=float, help="Galactic Latitude in degrees")
    pars.add_argument('-f', '--freq', default=408, type=float,
                      help="Desired frequency for the temperature in MHz "
                      "(default: 408)")
    pars.add_argument('-s', '--spind', default=-2.6, type=float,
                      help="Preferred spectral index for the sky (include "
                      "sign; default: -2.6)")

    return(vars(pars.parse_args()))


class SkyTemp:
    """
    Class defining the radio temperature of the sky at a given position.
    Designed to use the Remazeilles et al. 2014 reprocessing maps, or
    others with similar structure (Healpix pixellation).

    """

    def __init__(self, gl, gb, fits):
        """
        Parameters
        ----------
        gl - float, Galactic longitude
        gb - float, Galactic latitude
        fits - str, path to fits file containing map

        Useful attributes
        -----------------
        temp408 - float, the temperature of the sky position at 408 MHz

        """

        self.gl = gl
        self.gb = gb
        self.fits_file = fits
        self.table, self.meta = self.read_fits()
        self.flat_tab = self.table.flatten()
        self.temp408 = float(self.interp())

    def read_fits(self):
        t = Table.read(self.fits_file)
        return(t['TEMPERATURE'], t.meta)

    def get_pix(self, gl, gb):
        hp = HEALPix(nside=self.meta["NSIDE"],
                     order=self.meta["ORDERING"])
        return(hp.lonlat_to_healpix(gl*u.deg, gb*u.deg))

    def interp(self):
        gl = self.gl
        gb = self.gb
        # Want a few nearby values separated by beamsize
        beamsize = self.meta["BEAMSIZE"]/60
        lvals = np.array([gl-beamsize, gl, gl+beamsize])
        bvals = np.array([gb-beamsize, gb, gb+beamsize])
        xvals = self.get_pix(lvals, np.array([gb, gb, gb]))
        yvals = self.get_pix(np.array([gl, gl, gl]), bvals)
        temps_l = self.flat_tab[xvals]
        temps_b = self.flat_tab[yvals]

        cs_l = CubicSpline(lvals, temps_l)
        cs_b = CubicSpline(bvals, temps_b)
        return(np.mean([cs_l(gl), cs_b(gb)]))

    def get_temp(self, freq, spec_ind=-2.6):
        """
        Derive the temperature of the sky at the given frequency.

        Parameters
        ----------
        freq - float or numpy.ndarray of floats, the desired frequency/ies
            for conversion
        spec_ind - float, the spectral index of the sky temperature
            (Default: -2.6)

        Returns
        -------
        float or numpy.ndarray of floats (with the shape of `freq`)
            representing the temperature at the given frequency/ies

        """

        return(self.temp408*(freq/408)**spec_ind)


def main(fits, gl, gb, freq=408, spind=-2.6):
    print("gl      gb     freq  temp")
    temp = SkyTemp(gl, gb, fits).get_temp(freq, spind)
    print("{:<7.2f} {:<6.2f} {:<5.0f} {:.3f}".format(gl, gb, freq, temp))


if __name__ == "__main__":
    main(**proc_args())
