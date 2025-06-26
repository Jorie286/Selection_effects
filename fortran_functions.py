# equivalent file of fortran_functions.h and fortran_functions.c

def kggalmod(r, z):
    """
    Kuijken & Gilmore MNRAS 239 571-603 (89)
    returns the derivatives of the galactic potential

    Input:
        r, radial coordinate ???
        z, cartisian z coordinate
    """
    # masses in solar masses
    mdisk, mnuc, mbulge = 1.45e11, 9.3e9, 1e10
    beta = [0.4, 0.5, 0.1]

    # distances in Kpc
    h = [0.325, 0.090, 0.125]
    a=2.4
    b_disk, b_nuc, b_bulge = 5.5, 0.25, 1.5

    # gravitational constant time Msun
    gmsun = 4.498502167e-12 # kpc**3 Myr**-2

    # temporary variables for speed
    shz1=sqrt(h(1)*h(1)+z*z)
    shz2=sqrt(h(2)*h(2)+z*z)
    shz3=sqrt(h(3)*h(3)+z*z)
    b1shz=beta(1)*shz1
    b2shz=beta(2)*shz2
    b3shz=beta(3)*shz3

    # disk/ halo
    dpdr=(mdisk*r/(b_disk*b_disk + r*r + (a + b1shz + b2shz + b3shz)**2)**1.5)

    dpdz=(mdisk*(beta(1)*z/shz1 + beta(2)*z/shz2 +beta(3)*z/shz3)*(a + b1shz + b2shz +b3shz)/(b_disk*b_disk + r*r + (a + b1shz + b2shz + b3shz)**2)**1.5)

    # nucleus
    rbz = (z*z+b_nuc*b_nuc+r*r)**-1.5
    dpdz=dpdz+mnuc*z*rbz
    dpdr=dpdr+mnuc*r*rbz

    # bulge
    rbz = (z*z+b_bulge*b_bulge+r*r)**-1.5
    dpdz=(dpdz+mbulge*z*rbz)*gmsun
    dpdr=(dpdr+mbulge*r*rbz)*gmsun
