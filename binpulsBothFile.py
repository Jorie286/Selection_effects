# equivalent file of a binpuls____.c file

# NOTE: some of these files could be combined, I am keeping them separate for now so that they mirror the .c files better.

import numpy as np
import math
import uconst
import mymath
import gal2eq
import gal_cart
import fortran_functions
import seed_from_time
import misc_functions
import pulsar
import survey
import pulsar_survey_functions

NS_USED = 5 # the actual number of surveys used

def main(argc, *argv):
    a, e, m1, m2, porb, ce ,p_mer, a_mer, tf, l, p1, pdot, dt, tb, b, step, x, y, z, vx, vy, vz, q, r, area, e0, porb0
    """
    Main body of the code.

    Input:
        a,
        e,
        m1,
        m2,
        porb,
        ce,
        p_mer,
        a_mer,
        tf,
        l,
        p1,
        pdot,
        dt,
        tb,
        b,
        step,
        x,
        y,
        z,
        vx,
        vy,
        vz,
        q,
        r,
        area,
        e0,
        porb0,

    Returns:

    """
    if argc!=4:
        printf("Wrong usage, you have to pass: input output survey\n")
        break

    with open(argv[1], "r") as in_file:
        data = in_file.read()

    with open(argv[2], "w") as out_file: # put this at the end of the function to output the results to a file

