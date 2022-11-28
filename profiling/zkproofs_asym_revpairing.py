from charm.toolbox.pairinggroup import PairingGroup, ZR, G1, G2, GT, pair
import random
import time
import functools
import sys


group = PairingGroup('BN254')
f1, g1, h1 = [group.random(G2) for i in range(3)]
g2 = group.random(G1)
eg1g2 = pair(g2, g1)
ef1g2 = pair(g2, f1)
eh1g2 = pair(g2, h1)
iden = g1 ** 0
invg2 = g2 ** (-1)
invef1g2 = ef1g2 ** (-1)
inveg1g2 = eg1g2 ** (-1)
f1.initPP()
g1.initPP()
h1.initPP()
g2.initPP()
ef1g2.initPP()
eg1g2.initPP()
eh1g2.initPP()
invg2.initPP()
invef1g2.initPP()
inveg1g2.initPP()

#### Other globals ####
sigs = []
comms = []
phi = []
_rands = []

#### Utils ####


def sz(base64_input):
    b_padded = base64_input.split(str.encode(":"))[1]
    pad_size = b_padded.count(str.encode("="))
    b_len_without_pad = len(b_padded)-4
    byte_len = (b_len_without_pad * 3)/4 + (3-pad_size)-1
    bit_len = byte_len * 8
    return bit_len

#### Commitments ####


def commit(m, _r):
    return (g1**m) * (h1**_r)


def open(c, m, _r):
    return c == (g1**m) * (h1**_r)


def test_commitment():
    m = 3
    _r = group.random(ZR)
    c = commit(m, _r)
    cserial = group.serialize(c)
    cprime = group.deserialize(cserial)
    return open(cprime, m, _r)


#### BLS signatures ####

def blskeygen():
    _sk = group.random(ZR)
    pk = g2 ** _sk
    return _sk, pk


def blssign(m, _sk):
    mhash = group.hash(m, type=G2)
    sigma = mhash ** _sk
    return sigma


def blsverify(sigma, m, pk):
    mhash = group.hash(m, type=G2)
    return pair(mhash, pk) == pair(sigma, g2)


def test_bls():
    _sk, pk = blskeygen()
    m = b'Hello'
    sigma = blssign(m, _sk)
    return blsverify(sigma, m, pk)


#### BB signatures ####

def bbkeygen():
    _sk = group.random(ZR)
    pk = g2 ** _sk
    return _sk, pk


def bbsign(m, _sk):
    sigma = g1 ** (1 / (m + _sk))
    return sigma


def bbverify(sigma, m, pk):
    return pair(pk * (g2 ** m), sigma) == pair(g2, g1)


def test_bb():
    _sk, pk = bbkeygen()
    m = 5
    sigma = bbsign(m, _sk)
    return bbverify(sigma, m, pk)


#### BBS+ signatures ####

def bbspluskeygen():
    _sk = group.random(ZR)
    pk = g2 ** _sk
    return _sk, pk


def bbsplussign(m, _sk):
    c, r = group.random(ZR, 2)
    A = (f1 * (g1 ** m) * (h1 ** r)) ** (1 / (c + _sk))
    return (A, c, r)


def bbsplusverify(sigma, m, pk):
    A, c, r = sigma
    return pair(pk * (g2 ** c), A) == pair(g2, f1 * (g1 ** m) * (h1 ** r))


def bbsplussign_commitment(C, _sk):
    c, rdash = group.random(ZR, 2)
    A = (f1 * (h1 ** rdash) * C) ** (1 / (c + _sk))
    return (A, c, rdash)


def bbsplussign_randomise(sigma, r):
    A, c, rdash = sigma
    rdoubledash = rdash + r
    return (A, c, rdoubledash)


def test_bbsplus():
    _sk, pk = bbspluskeygen()
    m = 5
    sigma = bbsplussign(m, _sk)
    return bbsplusverify(sigma, m, pk)


def test_bbsplus_commitment():
    _sk, pk = bbspluskeygen()
    m = 5
    r = group.random(ZR)
    C = commit(m, r)
    sigma = bbsplussign_commitment(C, _sk)
    sigmasdoubledash = bbsplussign_randomise(sigma, r)
    return bbsplusverify(sigmasdoubledash, m, pk)


#### PoK1: PoK of a signature on a committed message ####
# PK{(rho, r, v): C = g^rho h^r and V = g^{v/(x+rho)}}

def pok1commit(V):
    _s, _t, _m = group.random(ZR, 3)
    a = (pair(g2, V) ** (-_s)) * (eg1g2 ** _t)
    D = (g1 ** _s) * (h1 ** _m)
    return _s, _t, _m, a, D


def pok1chal():
    return group.random(ZR)


def pok1resp(_s, _t, _m, _rho, _r, _v, c):
    zrho = _s - (_rho * c)
    zv = _t - (_v * c)
    zr = _m - (_r * c)
    return zrho, zv, zr


def pok1verify(C, V, a, D, zrho, zv, zr, y, c):
    cond1 = (D == (C ** c) * (h1 ** zr) * (g1 ** zrho))
    cond2 = (a == (pair(y, V) ** c) *
             (pair(g2, V) ** (-zrho)) * (eg1g2 ** (zv)))
    return cond1 and cond2


def pok1(C, V, _rho, _r, _v, y):
    _s, _t, _m, a, D = pok1commit(V)
    c = pok1chal()
    zrho, zv, zr = pok1resp(_s, _t, _m, _rho, _r, _v, c)
    return pok1verify(C, V, a, D, zrho, zv, zr, y, c)


def pok1nizkproof(C, V, _rho, _r, _v, y):
    stmt = (g1, h1, g2, C, V)
    _s, _t, _m, a, D = pok1commit(V)

    c = group.hash((stmt, (a, D)), type=ZR)

    zrho, zv, zr = pok1resp(_s, _t, _m, _rho, _r, _v, c)
    return c, (zrho, zv, zr)


def pok1nizkverify(C, V, y, pf):
    stmt = (g1, h1, g2, C, V)
    c, (zrho, zv, zr) = pf
    verif = ((C ** c) * (h1 ** zr) * (g1 ** zrho), (pair(y, V) ** c)
             * (pair(g2, V) ** (-zrho)) * (eg1g2 ** (zv)))
    return c == group.hash((stmt, verif), type=ZR)


#### Set membership proof (C commits some rho in set phi) ####

def genphi(n):
    return [group.random(ZR) for i in range(n)]


def gencomms(phi):
    n = len(phi)
    _shift = random.randint(0, n-1)
    _rands = []
    comms = []
    for i in range(n):
        j = (i + _shift) % n
        _rand = group.random(ZR)
        comm = commit(phi[j], _rand)
        _rands.append(_rand)
        comms.append(comm)
    return _shift, _rands, comms


def verfsigs(phi):
    _sk, pk = bbkeygen()
    sigs = []
    for m in phi:
        sig = bbsign(m, _sk)
        sigs.append(sig)
    return pk, sigs

def set_membership_nizkproof(comms, phi, y, i, _rands, sigs, _shift, n):
    C = comms[i]

    # Prover identifies the randomness and the signature corresponding to the i^th commitment
    _r = _rands[i]
    _j = (i + _shift) % n
    _sig = sigs[_j]
    _rho = phi[_j]

    # Prover sends a blinded signature on the shifted message
    _v = group.random(ZR)
    V = _sig ** _v

    # NIZK proof
    return V, pok1nizkproof(C, V, _rho, _r, _v, y)

def set_membership_nizkverify(comms, V_pfs, y, i):
    V, pf = V_pfs[i]
    C = comms[i]

    return pok1nizkverify(C, V, y, pf)


def test_set_membership(universal=False):
    n = 1000
    print("Set membership (given C, does it commit some value in set Phi of %s values):" % n)
    phi = genphi(n)
    _shift, _rands, comms = gencomms(phi)

    # Verifier sends fresh public key and signatures on all elements of phi
    group.InitBenchmark()
    group.StartBenchmark(["RealTime"])
    y, sigs = verfsigs(phi)
    group.EndBenchmark()
    print(" - Generating %s verifier signatures (sec):" %
          n, group.GetBenchmark("RealTime"))

    # Verifier decides to check i^th commitment in comms
    i = random.randint(0, n-1)

    # NIZK proof
    group.StartBenchmark(["RealTime"])
    V, pf = set_membership_nizkproof(comms, phi, y, i, _rands, sigs, _shift, n)
    group.EndBenchmark()
    print(" - [Online] Generating NIZK proofs for a single element (sec):",
          group.GetBenchmark("RealTime"))

    # Size of NIZK proof
    c, (zrho, zv, zr) = pf
    print(" - [Online] Size of NIZK proof for a single element (bits):",
        sz(group.serialize(V)) + sz(group.serialize(c)) + sz(group.serialize(zrho)) + sz(group.serialize(zv)) + sz(group.serialize(zr)))

    # NIZK verify
    group.StartBenchmark(["RealTime"])
    C = comms[i]
    status = pok1nizkverify(C, V, y, pf)
    group.EndBenchmark()
    print(" - [Online] Verifying NIZK proofs for a single element (sec):",
          group.GetBenchmark("RealTime"))

    if not universal:
        return status

    # Generating multiple NIZKs
    group.StartBenchmark(["RealTime"])
    V_pfs = []
    for i in range(len(comms)):
        V, pf = set_membership_nizkproof(comms, phi, y, i, _rands, sigs, _shift, n)
        V_pfs.append((V,pf))
    group.EndBenchmark()
    print(" - Generating %d NIZK proofs (sec):" % (n),
          group.GetBenchmark("RealTime"))

    # Verifying multiple NIZKs
    group.StartBenchmark(["RealTime"])
    for i in range(len(comms)):
        currstatus = set_membership_nizkverify(comms, V_pfs, y, i)
        status = status and currstatus
    group.EndBenchmark()
    print(" - Verifying %d NIZK proofs (sec):" % n,
        group.GetBenchmark("RealTime"))

    return status

#### POK2: PoK of a representation of a group element ####
# PK{(rho,r): C = g^rho h^r}


def pok2commit():
    _r1, _r2 = group.random(ZR, 2)
    a = (g1 ** _r1) * (h1 ** _r2)
    return _r1, _r2, a


def pok2resp(_r1, _r2, _rho, _r, c):
    z1 = _r1 - (_rho * c)
    z2 = _r2 - (_r * c)
    return z1, z2


def pok2nizkproof(C, _rho, _r):
    stmt = (g1, h1, C)
    _r1, _r2, a = pok2commit()
    c = group.hash((stmt, a), type=ZR)
    z1, z2 = pok2resp(_r1, _r2, _rho, _r, c)
    return c, (z1, z2)


def pok2nizkverify(C, pf):
    stmt = (g1, h1, C)
    c, (z1, z2) = pf
    verif = (C ** c) * (g1 ** z1) * (h1 ** z2)
    return c == group.hash((stmt, verif), type=ZR)


def pok2nizkproofs(comms, phi, _shift, _rands):
    n = len(comms)
    pfs = []
    for i, comm in enumerate(comms):
        _j = (i + _shift) % n
        _rho = phi[_j]
        _r = _rands[i]
        pf = pok2nizkproof(comm, _rho, _r)
        pfs.append(pf)
    return pfs


def pok2nizkverifies(comms, pfs):
    status = True
    for i, comm in enumerate(comms):
        pf = pfs[i]
        status = status and pok2nizkverify(comm, pf)
    return status


def test_pok2():
    _rho, _r = group.random(ZR, 2)
    C = commit(_rho, _r)
    pf = pok2nizkproof(C, _rho, _r)
    return pok2nizkverify(C, pf)

#### PoK3: PoK of a BBS+ signature on a given message ####
# PK{(s1, s2, c, r, d1, d2): B1 = g^s1 h^s2 and B1^c = g^d1 h^d2 and
# e(B2, y)e(f, f)^(-1) = e(B2, f)^{-c}e(g, y)^s2 e(g, f)^d2 e(g, f)^rho e(h, f)^r}


def pok3nizkproof(B1, B2, y, eg1y, rho, _s1, _s2, _c, _r, _d1, _d2):
    invB1 = B1 ** (-1)
    eB2y = pair(y, B2)
    inveB2g2 = pair(invg2, B2)
    egstar = eB2y * (invef1g2) * (inveg1g2 ** (rho))

    # Now proving PK{(s1, s2, c, r, d1, d2):
    #                  B1     = g^s1 h^s2             AND
    #                  iden   = invB1^c g^d1 h^d2     AND
    #                  egstar = inveB2g2^c eg1y^s2 eg1g2^d2 eh1g2^r }

    stmt = (f1, g1, h1, g2, invB1, inveB2g2,
            eg1y, eg1g2, eh1g2, B1, iden, egstar)

    # Commit (each ri corresponds to the respective known discrete log)
    _r1, _r2, _r3, _r4, _r5, _r6 = group.random(ZR, 6)
    com1 = (g1 ** _r1) * (h1 ** _r2)
    com2 = (invB1 ** _r3) * (g1 ** _r5) * (h1 ** _r6)
    com3 = (inveB2g2 ** _r3) * (eg1y ** _r2) * (eg1g2 ** _r6) * (eh1g2 ** _r4)

    # Challenge
    c = group.hash((stmt, (com1, com2, com3)), type=ZR)

    # Response
    z1 = _r1 - c * _s1
    z2 = _r2 - c * _s2
    z3 = _r3 - c * _c
    z4 = _r4 - c * _r
    z5 = _r5 - c * _d1
    z6 = _r6 - c * _d2

    return c, (z1, z2, z3, z4, z5, z6)


def pok3nizkverify(B1, B2, y, eg1y, rho, pf):
    invB1 = B1 ** (-1)
    eB2y = pair(y, B2)
    inveB2g2 = pair(invg2, B2)
    egstar = eB2y * (invef1g2) * (inveg1g2 ** (rho))

    # Now verifying PK{(s1, s2, c, r, d1, d2):
    #                  B1     = g^s1 h^s2             AND
    #                  iden   = invB1^c g^d1 h^d2     AND
    #                  egstar = inveB2g2^c eg1y^s2 eg1g2^d2 eh1g2^r }

    stmt = (f1, g1, h1, g2, invB1, inveB2g2,
            eg1y, eg1g2, eh1g2, B1, iden, egstar)

    c, (z1, z2, z3, z4, z5, z6) = pf

    v1 = (B1 ** c) * (g1 ** z1) * (h1 ** z2)
    v2 = (iden) * (invB1 ** z3) * (g1 ** z5) * (h1 ** z6)
    v3 = (egstar ** c) * (inveB2g2 ** z3) * \
        (eg1y ** z2) * (eg1g2 ** z6) * (eh1g2 ** z4)

    return c == group.hash((stmt, (v1, v2, v3)), type=ZR)


def test_pok3():
    rho = group.random(ZR)
    _x, y = bbspluskeygen()
    _A, _c, _r = bbsplussign(rho, _x)
    _s1, _s2 = group.random(ZR, 2)
    _d1 = _c * _s1
    _d2 = _c * _s2
    B1 = commit(_s1, _s2)
    B2 = _A * (g1 ** _s2)
    eg1y = pair(y, g1)
    pf = pok3nizkproof(B1, B2, y, eg1y, rho, _s1, _s2, _c, _r, _d1, _d2)
    return pok3nizkverify(B1, B2, y, eg1y, rho, pf)

#### Reverse set membership proof (rho is committed by some C in set phi') ####
# PK{(r): g^rho h^r \in Phi'}

def verfsigs_commitment(comms):
    _sk, pk = bbspluskeygen()
    sigs = []
    for comm in comms:
        sig = bbsplussign_commitment(comm, _sk)
        sigs.append(sig)
    return pk, sigs

def reverse_set_membership_nizkproof(phi, y, eg1y, j, _rands, sigs, _shift, n):
    rho = phi[j]

    # Prover identifies the index for the commitment committing rho
    _i = (j - _shift) % n
    _sig = sigs[_i]
    _rand = _rands[_i]

    # NIZK proof
    _sigdoubledash = bbsplussign_randomise(_sig, _rand)
    _A, _c, _r = _sigdoubledash
    _s1, _s2 = group.random(ZR, 2)
    _d1 = _c * _s1
    _d2 = _c * _s2
    B1 = commit(_s1, _s2)
    B2 = _A * (g1 ** _s2)
    return B1, B2, pok3nizkproof(B1, B2, y, eg1y, rho, _s1, _s2, _c, _r, _d1, _d2)

def reverse_set_membership_nizkverify(phi, B1B2_pfs, y, eg1y, j):
    B1, B2, pf = B1B2_pfs[j]
    rho = phi[j]

    return pok3nizkverify(B1, B2, y, eg1y, rho, pf)


def test_reverse_set_membership(universal=False):
    n = 1000
    print("Reverse set membership (given rho, is it committed by some commitment in set Phi' of %s commitments):" % n)
    phi = genphi(n)
    _shift, _rands, comms = gencomms(phi)

    # Prover proves knowledge of the committed values
    group.InitBenchmark()
    group.StartBenchmark(["RealTime"])
    pfs = pok2nizkproofs(comms, phi, _shift, _rands)
    status = pok2nizkverifies(comms, pfs)
    group.EndBenchmark()
    print(" - Proving knowledge of %s committed values (sec):" %
          n, group.GetBenchmark("RealTime"))

    # Verifier sends fresh public key and signatures on all elements of comms
    group.InitBenchmark()
    group.StartBenchmark(["RealTime"])
    y, sigs = verfsigs_commitment(comms)
    group.EndBenchmark()
    print(" - Generating %s verifier signatures (sec):" %
          n, group.GetBenchmark("RealTime"))

    # Verifier decides to check j^th element in phi
    j = random.randint(0, n-1)

    # NIZK proof generation
    group.StartBenchmark(["RealTime"])
    eg1y = pair(y, g1)
    # eg1y.initPP() - benefits of pre-computation do not show up for single element calcs
    B1, B2, pf = reverse_set_membership_nizkproof(phi, y, eg1y, j, _rands, sigs, _shift, n)
    group.EndBenchmark()
    print(" - [Online] Generating NIZK proof for a single element (sec):",
          group.GetBenchmark("RealTime"))

    # Size of the proof
    c, (z1, z2, z3, z4, z5, z6) = pf
    print(" - [Online] Size of NIZK proof for a single element (bits):", sz(group.serialize(B1)) + sz(group.serialize(B2)) + sz(group.serialize(c)) + sz(group.serialize(z1)) +
          sz(group.serialize(z2)) + sz(group.serialize(z3)) + sz(group.serialize(z4)) + sz(group.serialize(z5)) + sz(group.serialize(z6)))

    # NIZK proof verification
    group.StartBenchmark(["RealTime"])
    rho = phi[j]
    eg1y = pair(y, g1)
    #eg1y.initPP() - the benefits of pre-computation do not show up for single element calcs
    status = pok3nizkverify(B1, B2, y, eg1y, rho, pf)
    group.EndBenchmark()
    print(" - [Online] Verifying NIZK proof for a single element (sec):",
          group.GetBenchmark("RealTime"))

    if not universal:
        return status

    # Generating multiple NIZKs
    group.StartBenchmark(["RealTime"])
    B1B2_pfs = []
    eg1y = pair(y, g1)
    y.initPP()
    eg1y.initPP()
    for j in range(len(phi)):
        B1, B2, pf = reverse_set_membership_nizkproof(phi, y, eg1y, j, _rands, sigs, _shift, n)
        B1B2_pfs.append((B1,B2,pf))
    group.EndBenchmark()
    print(" - Generating %d NIZK proofs (sec):" % (n),
          group.GetBenchmark("RealTime"))

    # Verifying multiple NIZKs
    group.StartBenchmark(["RealTime"])
    eg1y = pair(y, g1)
    y.initPP()
    eg1y.initPP()
    for j in range(len(phi)):
        currstatus = reverse_set_membership_nizkverify(phi, B1B2_pfs, y, eg1y, j)
        status = status and currstatus
    group.EndBenchmark()
    print(" - Verifying %d NIZK proofs (sec):" % n,
        group.GetBenchmark("RealTime"))

    return status


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv)==2 else None
    print("Curve: %s (security equivalent to 3048 bits RSA security level)\n" %
          group.groupType())
    if mode == "all":
        print("Commitments ok?:", test_commitment())
        print("BLS signatures ok?:", test_bls())
        print("BB signatures ok?:", test_bb())
        print("BBS+ signatures ok?:", test_bbsplus())
        print("BBS+ signatures on committed messages ok?:", test_bbsplus_commitment())
    print("Set membership ok?:", test_set_membership(mode=="all"), '\n')
    if mode == "all":
        print("Knowledge of representation of group element ok?:", test_pok2())
        print("Knowledge of BBS+ signature ok?:", test_pok3())
    print("Reverse set membership ok?:", test_reverse_set_membership(mode=="all"))
