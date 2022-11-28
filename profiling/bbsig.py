import random
from charm.toolbox.pairinggroup import ZR, pair

from globals import group, g1, f2, eg1f2
from parutils import rank, serialize_bcast

def bbkeygen():
    _sk = group.random(ZR)
    pk = f2 ** _sk
    return _sk, pk

def bbsign(m, _sk):
    sigma = g1 ** (1 / (m + _sk))
    return sigma

def bbverify(sigma, m, pk):
    return pair(sigma, pk * (f2 ** m)) == eg1f2

def bbbatchverify(sigmas, ms, pk):
    # Choose random delta_i. Batch verification is thus verifying the following:
    #    prod_i e(sigma_i, pk * (f2 ** m_i))**delta_i = eg1f2 ** delta_i
    # which reduces to:
    #    prod_i (e(sigma_i, pk) ** delta_i) * (e(sigma_i, f2 ** m_i) ** delta_i) = prod_i eg1f2 ** delta_i
    #
    # The above can be efficiently checked by making only O(1) pairing computations:
    #    e(prod_i sigma_i ** delta_i, pk) * e(prod_i sigma_i ** (m_i * delta_i), f2) = eg1f2 ** (sigma_i delta_i)
    #
    # Ref: Anna Lisa Ferrara, Matthew Green, Susan Hohenberger, ``Practical Short Signature Batch Verification'', https://eprint.iacr.org/2008/015.pdf

    deltas = [random.getrandbits(80) for _ in range(len(sigmas))]

    sigma_delta_prod = g1 ** 0
    sigma_mdelta_prod = g1 ** 0
    delta_sum = 0

    for i in range(len(sigmas)):
        sigma_delta = sigmas[i] ** deltas[i]
        sigma_mdelta = sigma_delta ** ms[i]
        sigma_delta_prod = sigma_delta_prod * (sigma_delta)
        sigma_mdelta_prod = sigma_mdelta_prod * (sigma_mdelta)
        delta_sum = delta_sum + deltas[i]

    return pair(sigma_delta_prod, pk) * pair(sigma_mdelta_prod, f2) == eg1f2 ** delta_sum

def bbkeygen_par():
    if rank == 0:
        _sk, pk = bbkeygen()
    else:
        _sk, pk = None, None
 
    _mysk, mypk = serialize_bcast(_sk, pk)
    return _mysk, mypk