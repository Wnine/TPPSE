import random,time

import galois
import numpy
from sympy import symbols, div, GF



def tppse_setup(ashe,L,NW,roots):
    average = NW//(NW//16)
    EQx = [None] * average
    Sx = [None] * average
    Kx = [None] * average
    x = symbols('x')
    for i in range(average):
        pwi_poly = ashe.create_poly_sympy(L - 1, roots)
        rploy = ashe.random_poly_sympy(L - 1)
        Qx, _ = div(rploy * x ** (2*128), pwi_poly, domain = GF(ashe.p))
    #    print(Qx)
        Qx = [int(coeff) for coeff in Qx.all_coeffs()]
        Qx = galois.Poly(Qx, field=ashe.GF)
        EQx[i]= ashe.enc_poly(Qx)
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

def tppse_search(ashe, EDB,NW):
    # 1.令牌生成算法
    m1 = random.randint(1, ashe.p)
    TQ= []
    average = NW//(NW//16)
    for i in range(average):
        ei = random.randint(1, ashe.p)
        si = ashe.hash_K(i)
        ki = ashe.hash_K(i)
        vi = (m1 * (ei + ki) + si) % ashe.p
        TQ.append((vi,ei))
    # 2.搜索算法
    EQx = EDB["EQx"]
    vQx,vQx_counter,eQx,eQx_counter =[],[],[],[]
    for i in range(average):
        vQwix = ashe.hom_scal_poly(TQ[i][0], EQx[4])
        eQwix = ashe.hom_scal_poly(TQ[i][1], EQx[4])
        vQx.append(vQwix[0])
        vQx_counter.append(vQwix[1])
        eQx.append(eQwix[0])
        eQx_counter.append(eQwix[1])
    EPSigma = ashe.hom_add_poly(EDB["encSumSx"],ashe.sum_add_poly(vQx,vQx_counter))
    EPSigma_pie=ashe.hom_add_poly(EDB["encSumKx"],ashe.sum_add_poly(eQx,eQx_counter))
    # 3.解密算法
    PSigma1_arr, PSigma1 = ashe.dec_poly(EPSigma)
    PSigma1 = PSigma1*(ashe.GF(m1)**(-1))
    id_domain = [i+1 for i in range(2**17)]
    length = 128
    padding_length = length - len(id_domain) % length
    id_domain += [0] * padding_length
    for i in range(0, len(id_domain), length):
        ashe.multi_points_fft(PSigma1_arr,id_domain[i:i+length])

db = Database()
ashe = Ashe(128)
search_arr = []
EDB = tppse_setup(ashe, db.L,db.Nw_arr[i],roots_arr)
tppse_search(ashe, EDB,db.Nw_arr[i])
