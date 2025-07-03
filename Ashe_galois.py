import hashlib, random
from collections import Counter
from numpy import iterable

import galois
from gmpy2 import mpz
from sympy import ntt, intt, Poly, symbols, GF, isprime, primitive_root  # 符号计算库
from Crypto.Random import get_random_bytes
import numpy as np
from sympy.utilities.iterables import ibin
from sympy.utilities.misc import as_int


class Ashe:
    # 初始化
    def __init__(self, flag):
        if flag == 128:
            self.p = 850705917302346158658436518579420528641 #5*(2**127)+1
            self.g =  17
        if flag == 256:
            self.p = 364250417292324200797399107529409794995451282322276160234714826826044553475976593409 #3*(2**276)+1
            self.g = 19
        self.GF = galois.GF(self.p)
        self.K = get_random_bytes(16)

    # 定义：带密钥的哈希函数
    def hash_K(self, i):
        i_bytes = i.to_bytes(16, byteorder='big')
        data = i_bytes + self.K
        hash_value = hashlib.sha256(data).hexdigest()
        hash_int = int(hash_value, 16) % self.p
        return hash_int

    # 定义：加密一个多项式 ,随机生成一个标识符i只用在setup.py中
    def enc_poly(self, plain_poly):
        # 预先分配位置
        length = len(plain_poly.coefficients())
        enc_poly_arr = [0] * length
        counter = [Counter({})] * length
        for i, coeff in enumerate(plain_poly.coefficients()):
            identity = random.randint(1, 1000000)
            enc_coeff = (int(coeff) + self.hash_K(identity)) % self.p
            enc_poly_arr[i] = enc_coeff
            counter[i] = Counter({identity: 1})
        return galois.Poly(enc_poly_arr, field=self.GF), np.array(counter)

    def enc_poly_sympy(self, plain_poly):
        # 预先分配位置
        length = len(plain_poly.all_coeffs())
        enc_poly_arr = [0] * length
        counter = [Counter({})] * length
        for i, coeff in enumerate(plain_poly.all_coeffs()):
            identity = random.randint(1, 1000000)
            enc_coeff = (int(coeff) + self.hash_K(identity)) % self.p
            enc_poly_arr[i] = enc_coeff
            counter[i] = Counter({identity: 1})
        return galois.Poly(enc_poly_arr, field=self.GF), np.array(counter)

# 输入两个sympy.Poly 多项式  只用在VSSEsetup.py中
    def sympy_NTT(self, rx_poly, px_poly):
        if len(rx_poly.all_coeffs()) & (len(rx_poly.all_coeffs()) - 1):
            print("输入的多项式系数不是2的次幂")
        length_rx = len(rx_poly.all_coeffs())
        m = 2 ** (length_rx.bit_length())
        rxpoly_arr = [mpz(i) for i in rx_poly.all_coeffs()] + [mpz(0)] * (m - length_rx)
        pxpoly_arr = [mpz(i) for i in px_poly.all_coeffs()] + [mpz(0)] * (m - length_rx)
        ntt_rx = self.sympy_ntt(rxpoly_arr)
        ntt_pxpoly_arr = self.sympy_ntt(pxpoly_arr)
        c_arr = [a * b % self.p for a, b in zip(ntt_rx, ntt_pxpoly_arr)]
        rx_mul_px = self.sympy_ntt(c_arr, inverse=True)
        return galois.Poly(rx_mul_px[:-1], field=self.GF)

    def sympy_ntt(self, seq, inverse=False):
        """Utility function for the Number Theoretic Transform"""
        if not iterable(seq):
            raise TypeError("Expected a sequence of integer coefficients "
                            "for Number Theoretic Transform")
        p = self.p
        # 将输入序列的每个元素转换为模p的余数
        a = [as_int(x) % p for x in seq]
        n = len(seq)
        if n < 1:
            return a
        b = n.bit_length() - 1
        if n & (n - 1):
            b += 1
            n = 2 ** b
        if (p - 1) % n:
            raise ValueError("Expected prime modulus of the form (m*2**k + 1)")
        a += [0] * (n - len(a))
        for i in range(1, n):
            j = int(ibin(i, b, str=True)[::-1], 2)
            if i < j:
                a[i], a[j] = a[j], a[i]
        pr = self.g
        rt = pow(pr, (p - 1) // n, p)
        if inverse:
            rt = pow(rt, p - 2, p)
        w = [1] * (n // 2)
        for i in range(1, n // 2):
            w[i] = w[i - 1] * rt % p
        h = 2
        while h <= n:
            hf, ut = h // 2, n // h
            for i in range(0, n, h):
                for j in range(hf):
                    u, v = a[i + j], a[i + j + hf] * w[ut * j]
                    a[i + j], a[i + j + hf] = (u + v) % p, (u - v) % p
            h *= 2
        if inverse:
            rv = pow(n, p - 2, p)
            a = [x * rv % p for x in a]

        return a

    # 定义：一个加密多项式 + 一个加密多项式
    def hom_add_poly(self, encf, encg):
        return encf[0] + encg[0], np.add(encf[1], encg[1])

    # 实现n个加密多项式的相加
    def sum_add_poly(self, encf_arr, counter_arr):
        return np.sum(encf_arr), np.sum(counter_arr, axis=0)

    # 定义：一个标量 * 一个加密多项式
    def hom_scal_poly(self, k, en_poly):
        counter = en_poly[1]
        for c in counter:
            for key in c:
                c[key] *= k
        return k * en_poly[0], counter

    # 定义：解密一个多项式
    def dec_poly(self, encpoly):
        poly_arr = []
        cipher_poly = encpoly[0]
        counter = encpoly[1]
        for i, coeff in enumerate(cipher_poly.coefficients()):
            coeff_to_counter = counter[i]
            hash_acc = 0
            for identity in coeff_to_counter:
                hash_acc = hash_acc + (coeff_to_counter[identity] * self.hash_K(identity) % self.p) % self.p
            plain_coeff = (int(coeff) + hash_acc) % self.p
            poly_arr.append(plain_coeff)
        return poly_arr, galois.Poly(poly_arr, field=self.GF)

    # 多项式在多点计算的 O(nlogn)算法
    def multi_points_fft(self, vals, domain):
        if len(vals) == 1:
            return vals
        L = self.multi_points_fft(vals[::2], domain[::2])
        R = self.multi_points_fft(vals[1::2], domain[::2])
        o = [0 for i in vals]
        for i, (x, y) in enumerate(zip(L, R)):
            y_times_root = y * domain[i] % self.p
            o[i] = (x + y_times_root) % self.p
            o[i + len(L)] = (x - y_times_root) % self.p
        return o

    def horner_evaluation(self, poly, x):
        coeffs = [int(i) for i in poly.coefficients()]

        result = 0
        for coeff in reversed(coeffs):
            result = ((result * x) % self.p + coeff) % self.p
        return result

    def random_enc_poly(self, L):
        random_poly = galois.Poly.Random(L, None, self.GF)
        return self.enc_poly(random_poly)

    def create_poly_sympy(self, degree, roots):
        n = len(roots)
        coefficients = [1]
        for i in range(n):
            coefficients[0] *= -roots[i] % self.p  # 根据韦达定理，计算多项式系数
            for j in range(i, 0, -1):
                coefficients[j] = (coefficients[j - 1] - (roots[i] * coefficients[j]) % self.p) % self.p
            coefficients.append(1)
        poly_arr = coefficients[::-1] + [0] * (degree - len(roots))
        return Poly(poly_arr, symbols('x'), domain=GF(self.p))

    def random_poly_sympy(self, L):
        return Poly([random.randint(1, self.p) for _ in range(L + 1)], symbols('x'), domain=GF(self.p))
