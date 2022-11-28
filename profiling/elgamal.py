from charm.toolbox.pairinggroup import ZR
from globals import group, g1

from parutils import rank, serialize_bcast

def elgamal_keygen():
    sk = group.random(ZR)
    pk = g1 ** sk
    return sk, pk

def elgamal_encrypt(pk, m):
    r  = group.random(ZR)
    return (g1 ** r, m * (pk ** r))

def elgamal_encrypt_with_randomness(pk, m,r):
    return (g1 ** r, m * (pk ** r))

def elgamal_decrypt(sk,c):
    c1, c2 = c 
    return c2 / (c1 ** sk)

def elgamal_mult(c1, c2):
    return (c1[0]*c2[0], c1[1]*c2[1])

def elgamal_div(c1, c2):
    return (c1[0]/c2[0], c1[1]/c2[1])

def elgamal_exp(c1, a):
    return (c1[0] ** a, c1[1] ** a)

def elgamal_reencrypt(pk, c):
    c_iden = elgamal_encrypt(pk, g1 ** 0)
    return elgamal_mult(c, c_iden)

def elgamal_reencrypt_with_randomness(pk, c ,r):
    c_iden = elgamal_encrypt_with_randomness(pk, g1 ** 0,r)
    return elgamal_mult(c, c_iden)

def elgamal_keygen_par():
    if rank == 0:
        elg_sk, elg_pk = elgamal_keygen()
    else:
        elg_sk, elg_pk = None, None
    elg_sk, elg_pk = serialize_bcast(elg_sk, elg_pk)
    return elg_sk, elg_pk
