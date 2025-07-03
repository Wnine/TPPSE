import galois
import numpy
import time

from sympy import symbols, div, GF



def tppse_setup(ashe,Lmax,n,roots):
    average = n//(n//16)
    EQx = [None] * average
    Sx = [None] * average
    Kx = [None] * average
    x = symbols('x')
    for i in range(average):
        pwi_poly = ashe.create_poly_sympy(Lmax - 1, roots)
        rploy = ashe.random_poly_sympy(Lmax - 1)
        Qx, _ = div(rploy * x ** (2*128), pwi_poly, domain = GF(ashe.p))
        Qx = [int(coeff) for coeff in Qx.all_coeffs()]
        Qx = galois.Poly(Qx, field=ashe.GF)
        enc_pwi = ashe.enc_poly(Qx)
        EQx.append(enc_pwi)
        si = ashe.hash_K(i)
        ki = ashe.hash_K(i)
        Sx[i]= (-si % ashe.p)* Qx
        Kx[i]= ki * Qx
    sumSx = numpy.sum(Sx)
    sumKx = numpy.sum(Kx)
    encSumSx = ashe.enc_poly(sumSx)
    encSumKx = ashe.enc_poly(sumKx)
    EDB = {"EQx":EQx,"encSumSx":encSumSx,"encSumKx":encSumKx,"K":ashe.K,}
    return EDB

db = Database()
ashe = Ashe(128)
tppse_setup(ashe,db.L,db.Nw_arr[i],DB[w])


