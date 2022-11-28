from charm.toolbox.pairinggroup import ZR

from globals import group
from misc import timed

beaver_a_shares = []
beaver_b_shares = []
beaver_c_shares = []
beaver_callcount = 0

#### Simple alpha-out-of-alpha secret sharing scheme

def sharerand(alpha):
    shares = [0]*alpha
    for a in range(alpha):
        shares[a] = group.random(ZR)
    return shares

def share(val, alpha):
    shares = [0]*alpha
    sum_shares = group.init(ZR, 0)
    for a in range(alpha-1):
        shares[a] = group.random(ZR)
        sum_shares = sum_shares + shares[a]
    shares[alpha-1] = val - sum_shares
    return shares

def reconstruct(shares):
    val = group.init(ZR, 0)
    for share in shares:
        val = val + share
    return val

def sharemult(shares1, shares2, alpha):
    global beaver_callcount
    dshares, eshares = [], []
    for a in range(alpha):
        dshares_a = shares1[a] - beaver_a_shares[beaver_callcount][a]
        eshares_a = shares2[a] - beaver_b_shares[beaver_callcount][a]
        dshares.append(dshares_a)
        eshares.append(eshares_a)

    d = reconstruct(dshares)
    e = reconstruct(eshares)

    multshares = []
    for a in range(alpha):
        multshare_a = (
            d * beaver_b_shares[beaver_callcount][a] + 
            e * beaver_a_shares[beaver_callcount][a] + 
            beaver_c_shares[beaver_callcount][a]
        )
        if a == 0:
            multshare_a += d*e
        multshares.append(multshare_a)

    beaver_callcount += 1

    return multshares

@timed
def gen_beaver_triples(myn, alpha):
    global beaver_a_shares, beaver_b_shares, beaver_c_shares
    for i in range(myn):
        beaver_a = group.random(ZR)
        beaver_b = group.random(ZR)
        beaver_c = beaver_a * beaver_b
        beaver_a_shares.append(share(beaver_a, alpha))
        beaver_b_shares.append(share(beaver_b, alpha))
        beaver_c_shares.append(share(beaver_c, alpha))

if __name__ == "__main__":
    alpha = 2
    gen_beaver_triples(1, alpha)
    shares1 = share(group.init(ZR, 3), alpha)
    shares2 = share(group.init(ZR, 5), alpha)
    multshares = sharemult(shares1, shares2, alpha)
    mul = reconstruct(multshares)
    print("Mul", mul)