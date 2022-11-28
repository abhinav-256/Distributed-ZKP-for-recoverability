from charm.toolbox.pairinggroup import PairingGroup, ZR, G1, G2, GT, pair
import random
import time
import functools
import sys
import io
import os
import pickle
from math import log2
from mpi4py import MPI
from collections.abc import Iterable

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()
strios = []

#### Parallel utils ####

def init_strios():
    for p in range(nprocs):
        strio = io.StringIO()
        strios.append(strio)
    strios[rank].write("******** Worker %s **********\n" % rank)

def print_strios():
    strios[rank].write("---------------------------\n")
    print(strios[rank].getvalue())

def close_strios():
    strios[rank].close()

def pprint(*x):
    initprint = True
    if nprocs == 1:
        print(*x)
    else:
        print(*x, file=strios[rank])

def get_count(items):
    """ Get the size of each sub-task. """
    ave, res = divmod(len(items), nprocs)
    return [ave + 1 if p < res else ave for p in range(nprocs)]

def parsplit(items):
    """ Split the list of items to a list of sublists such that
    the i^th sublist contains the items intended for the i^th
    process. """
    count = get_count(items)

    itemsplits = []
    last = 0
    for p in range(nprocs):
        itemsplit = items[last:last+count[p]]
        itemsplits.append(itemsplit)
        last = last + count[p]

    return itemsplits

def parunsplit(itemsplits):
    """ Unsplit a list of sublists to a single list. """
    items = []
    for itemsplit in itemsplits:
        items.extend(itemsplit)
    return items

def serialize_wrapper(item):
    """ If item is a group element, call group.serialize; otherwise
    serialize using pickle. """
    if item.__class__.__name__ == "Element":
        return group.serialize(item)
    elif isinstance(item, Iterable):
        return [serialize_wrapper(iitem) for iitem in item]
    else:
        return item

def deserialize_wrapper(sitem):
    """ Try to deserialize the serialized item as a group element.
    If it fails, deserialize using pickle. """
    try:
        item = group.deserialize(sitem)
    except Exception:
        if isinstance(sitem, Iterable):
            item = [deserialize_wrapper(siitem) for siitem in sitem]
        else:
            item = sitem
    return item

def serialize_bcast(*items):
    """ Serialize the items and broadcast to each process. """
    if rank == 0:
        sitems = []
        for item in items:
            sitems.append(serialize_wrapper(item))
    else:
        sitems = None

    sitems = comm.bcast(sitems, root=0)

    myitems = []
    for sitem in sitems:
        myitem = deserialize_wrapper(sitem)
        myitems.append(myitem)

    return myitems[0] if len(myitems) == 1 else myitems

def serialize_scatter(*itemlists):
    """ Serialize each item in each itemlist and scatter them to each process. """
    start = time.time()
    if rank == 0:
        sitemlists = []
        psitemlists = []
        for i, itemlist in enumerate(itemlists):
            t1 = time.time()
            sitemlist = [serialize_wrapper(item) for item in itemlist]
            t2 = time.time()
            # pprint("            Serialised itemlist %s of length %s:" % (i, len(sitemlist)), t2 - t1, "s")
            psitemlist = parsplit(sitemlist)
            t3 = time.time()
            # pprint("            Parsplit itemlist %s:" % i, t3 - t2, "s")
            sitemlists.append(sitemlist)
            psitemlists.append(psitemlist)
            t4  = time.time()
            # pprint("            Append itemlist %s:" % i, t4 - t3, "s")
        t5 = time.time()
        zipped_psitemlists = zip(*psitemlists)
        t6 = time.time()
        # pprint("        Zip itemlists:", t6 - t5, "s")
    else:
        zipped_psitemlists = None
    barrier_start = time.time()
    comm.barrier()
    # pprint("        Serialisation barrier:", time.time() - barrier_start, "s")
    # pprint("        Total serialisation and zipping:", time.time() - start, "s")

    scatter_start = time.time()
    mysitemlists = comm.scatter(zipped_psitemlists, root=0)
    scatter_end = time.time()
    # pprint("        Actual MPI scatter:", scatter_end - scatter_start, "s")

    deserialisation_start = time.time()
    myitemlists = []
    for i, mysitemlist in enumerate(mysitemlists):
        t1 = time.time()
        myitemlist = [deserialize_wrapper(sitem) for sitem in mysitemlist]
        # pprint("            Deserialising myitemlist %s of length %s:" % (i, len(myitemlist)), time.time() - t1, "s")
        myitemlists.append(myitemlist)
    barrier_start = time.time()
    comm.barrier()
    # pprint("        Deserialisation barrier:", time.time() - barrier_start, "s")
    # pprint("        Total deserialisation:", time.time() - deserialisation_start, "s")

    return myitemlists[0] if len(myitemlists) == 1 else myitemlists

def serialize_gather(*myitemlists):
    """ Serialize each item in myitemlist and gather at the root. """
    smyitemlists = []
    for myitemlist in myitemlists:
        smyitemlist = [serialize_wrapper(item) for item in myitemlist]
        smyitemlists.append(smyitemlist)
    zipped_smyitemlist = list(zip(*smyitemlists))
    zipped_psitemlist = comm.gather(zipped_smyitemlist, root=0)

    if rank == 0:
        zitemlist = parunsplit(zipped_psitemlist)
        sitemlists = list(zip(*zitemlist)) # Actually this is the way to "unzip" each item
        itemlists = []
        for sitemlist in sitemlists:
            itemlist = [deserialize_wrapper(sitem) for sitem in sitemlist]
            itemlists.append(itemlist)
        return itemlists[0] if len(itemlists) == 1 else itemlist
    else:
        return None

#### Other utils ####

def statusstr(status):
    return "(ok)" if status else "(not ok)"

def expand_scientific_notation(numstr):
    was_neg = False
    if not ("e" in numstr):
        return numstr
    if numstr.startswith('-'):
        numstr = numstr[1:]
        was_neg = True 
    str_vals = str(numstr).split('e')
    coef = float(str_vals[0])
    exp = int(str_vals[1])
    return_val = ''
    if int(exp) > 0:
        return_val += str(coef).replace('.', '')
        return_val += ''.join(['0' for _ in range(0, abs(exp - len(str(coef).split('.')[1])))])
    elif int(exp) < 0:
        return_val += '0.'
        return_val += ''.join(['0' for _ in range(0, abs(exp) - 1)])
        return_val += str(coef).replace('.', '')
    if was_neg:
        return_val='-'+return_val
    return return_val

def fmt(num):
    num_truncated_precision = f'{num:.2}'
    return expand_scientific_notation(num_truncated_precision)

#### Group initialisation ####

def generators():
    f1, g1, h1 = [group.random(G1) for i in range(3)]
    f2 = group.random(G2)
    ef1f2 = pair(f1, f2)
    eg1f2 = pair(g1, f2)
    eh1f2 = pair(h1, f2)
    invf1 = f1 ** (-1)
    invg1 = g1 ** (-1)
    invh1 = h1 ** (-1)
    invf2 = f2 ** (-1)
    invef1f2 = ef1f2 ** (-1)
    inveg1f2 = eg1f2 ** (-1)
    inveh1f2 = eh1f2 ** (-1)
    iden = g1 ** 0
    return (f1, g1, h1, f2, ef1f2, eg1f2, eh1f2, invf1, invg1, invh1, invf2, invef1f2, inveg1f2, inveh1f2, iden)

def initgens():
    genserials = [group.serialize(gen) for gen in generators()] if rank == 0 else None
    genserials = comm.bcast(genserials, root=0)
    gens = [group.deserialize(genserial) for genserial in genserials]
    for gen in gens:
        gen.initPP()
    return gens

group = PairingGroup('BN254')
genstart = time.time()
f1, g1, h1, f2, ef1f2, eg1f2, eh1f2, invf1, invg1, invh1, invf2, invef1f2, inveg1f2, inveh1f2, iden = initgens()
genend = time.time()

#### Commitments ####

def commit(m, _r):
    return (g1**m) * (h1**_r)


def open(c, m, _r):
    return c == (g1**m) * (h1**_r)

#### NEW #######BB signatures #### #### REMOVED

#### BBS+ signatures####

def bbspluskeygen():
    _sk = group.random(ZR)
    pk = f2 ** _sk
    return _sk, pk

def bbsplussign(m, _sk):
    c, r = group.random(ZR, 2)
    S = (f1 * (g1 ** m) * (h1 ** r)) ** (1 / (c + _sk))
    return (S, c, r)

def bbsplusverify(sigma, m, pk):
    S, c, r = sigma
    return pair(S, pk * (f2 ** c)) == pair(f1 * (g1 ** m) * (h1 ** r), f2)

def bbsplussign_commitment(C, _sk):
    c, rdash = group.random(ZR, 2)
    S = (f1 * (h1 ** rdash) * C) ** (1 / (c + _sk))
    return (S, c, rdash)

def bbsplussign_randomise(sigma, r):
    S, c, rdash = sigma
    rdoubledash = rdash + r
    return (S, c, rdoubledash)

def bbsplusbatchverify(sigmas, ms, pk):
    # Consider that the signature sigma can be extracted as S, c, r = sigma.
    # Choose random delta_i. Batch verification is thus verifying the following:
    #    prod_i e(S_i, pk )**delta_i * e((S_i ** c_i)( g1 ** (-m_i))(h1 ** (-r_i)), f2) ** delta_i = prod_i ef1f2 ** delta_i
    #
    # The above can be efficiently calculated using a single pairing computation:
    #    e(prod_i S_i ** delta_i, pk) * e(prod_i ((S_i ** c_i)( g1 ** (-m_i))(h1 ** (-r_i))) ** delta_i, f2) = ef1f2 ** (sum_i delta_i)
    #
    # Ref: Anna Lisa Ferrara, Matthew Green, Susan Hohenberger, ``Practical Short Signature Batch Verification'', https://eprint.iacr.org/2008/015.pdf

    deltas = [random.getrandbits(80) for _ in range(len(sigmas))]

    S_delta_prod = g1 ** 0
    Sgh_delta_prod = g1 ** 0
    delta_sum = 0
    q = group.order()

    for i in range(len(sigmas)):
        S, c, r = sigmas[i]
        m = ms[i]
        S_delta = S ** deltas[i]
        Sgh_delta = ((S ** c) * (g1 ** (q-m)) * (h1 ** (q-r))) ** deltas[i]
        S_delta_prod = S_delta_prod * (S_delta)
        Sgh_delta_prod = Sgh_delta_prod * (Sgh_delta)
        delta_sum = delta_sum + deltas[i]

    return pair(S_delta_prod, pk) * pair(Sgh_delta_prod, f2) == ef1f2 ** delta_sum

#### Generating the set and its commitments ####

def genvals(n): ## generating values
    genvals_start = time.time()
    vals = [group.random(ZR) for _ in range(n)]
    genvals_end = time.time()
    t_genvals = genvals_end - genvals_start
    return vals, t_genvals

def gencomms(vals): ## commiting values
    gencomms_start = time.time()
    n = len(vals)
    _shift = random.randint(0, n-1)
    _rands = []
    comms = []
    for i in range(n):
        j = (i + _shift) % n
        _rand = group.random(ZR)
        comm = commit(vals[j], _rand)
        _rands.append(_rand)
        comms.append(comm)
    gencomms_end = time.time()
    t_gencomms = gencomms_end - gencomms_start
    return _shift, _rands, comms, t_gencomms

#### NEW##################Set membership proof (C commits some rho in set vals) ####  #### REMOVED


###
##NEW##test_set_membership REMOVED
####

#### Reverse set membership proof (rho is committed by some C in set vals') ####

def PoKcommproof(C, _t, _r):               #### RENAME
    """ PoK of the opening of the commitment, i.e.:
     PK{(t, r): C = g^t h^r} """
    stmt = (g1, h1, C)

    # Commit
    _r1, _r2 = group.random(ZR, 2)
    a = (g1 ** _r1) * (h1 ** _r2)

    # Challenge
    c = group.hash(stmt + (a,), type=ZR)

    # Response
    z1 = _r1 - (_t * c)
    z2 = _r2 - (_r * c)

    return c, (z1, z2)

def PoKcommverify(C, pf):
    stmt = (g1, h1, C)
    c, (z1, z2) = pf
    verif = (C ** c) * (g1 ** z1) * (h1 ** z2)
    return c == group.hash(stmt +(verif,), type=ZR)

def get_ivals(vals, _shift):
    if rank == 0:
        getivals_start = time.time()
        n = len(vals)
        ivals = []
        for i in range(n):
            _j = (i + _shift) % n
            m = vals[_j]
            ivals.append(m)
        getivals_end = time.time()
        t_getivals = getivals_end - getivals_start
        pprint("   * [T:root] Gathering per-core committed messages:", fmt(t_getivals), "s")
    else:
        ivals = None
        t_getivals = 0
    return ivals, t_getivals

def PoKcommproofs(comms, ivals, _rands):
    n = len(comms)
    pr = []
    for i, comm in enumerate(comms):
        _t = ivals[i]
        _r = _rands[i]
        pf = PoKcommproof(comm, _t, _r)
        pr.append(pf)
    return pr

def PoKcommverifies(comms, pr):
    status = True
    for i, comm in enumerate(comms):
        pf = pr[i]
        status = status and PoKcommverify(comm, pf)
    return status

def PoKcommproofs_par(comms, ivals, _rands):
    map_start = time.time()
    mycomms, myivals, _myrands = serialize_scatter(comms, ivals, _rands)
    map_end = time.time()
    t_map = map_end - map_start
    pprint("      * [T] Scattering comms, ivals, _rands to worker nodes:", fmt(t_map), "s")
    pprint("   * [T] Distributing data to worker nodes:", fmt(t_map), "s")

    genproofs_start = time.time()
    mypr = PoKcommproofs(mycomms, myivals, _myrands)
    genproofs_end = time.time()
    t_genproofs = genproofs_end - genproofs_start
    pprint("   * [T] Generating %s NIZK PoKs of commitments:" % len(mypr), fmt(t_genproofs), "s")

    verproofs_start = time.time()
    status = PoKcommverifies(mycomms, mypr)
    verproofs_end = time.time()
    t_verproofs = verproofs_end - verproofs_start
    pprint("   * [A] Verifying %s NIZK PoKs of commitments:" % len(mypr), fmt(t_verproofs), "s", statusstr(status))

    return status, t_map, t_genproofs, t_verproofs

def verfsigs_commitment(comms, _sk):
    sigs = []
    for comm in comms:
        sig = bbsplussign_commitment(comm, _sk)
        sigs.append(sig)
    return sigs

def verify_verfsigs_commitment(ivals, y, sigs, _rands):
    rhos = []
    rsigs = []
    for i, sig in enumerate(sigs):
        _t = ivals[i]
        _r = _rands[i]
        rsigs.append(bbsplussign_randomise(sigs[i], _r))
        rhos.append(_t)
    return bbsplusbatchverify(rsigs, rhos, y)

def verfsigs_commitment_par(comms, ivals, _rands):
    map_start = time.time()
    if rank == 0:
        _sk, pk = bbspluskeygen()
    else:
        _sk, pk = None, None
    bcast_start = time.time()
    _mysk, mypk = serialize_bcast(_sk, pk)
    bcast_end = time.time()
    pprint("      = [A] Broadcasting sk,pk to worker nodes", fmt(bcast_end - bcast_start), "s")
    scatter_start = time.time()
    mycomms, myivals, _myrands = serialize_scatter(comms, ivals, _rands)
    scatter_end = time.time()
    pprint("      = [A] Scattering comms, ivals, _rands to worker nodes", fmt(scatter_end - scatter_start), "s")
    map_end = time.time()
    t_map = map_end - map_start
    pprint("   * [A] Distributing data to worker nodes:", fmt(t_map), "s")

    # Generating signatures at the verifier
    verfsigs_start = time.time()
    mysigs = verfsigs_commitment(mycomms, _mysk)
    verfsigs_end = time.time()
    t_gensigs = verfsigs_end - verfsigs_start
    pprint("   * [A] Generating %s verifier BBS+ signatures:" % len(mysigs), fmt(t_gensigs), "s")

    # Verifying signatures at the prover
    verify_verfsigs_start = time.time()
    status = verify_verfsigs_commitment(myivals, mypk, mysigs, _myrands)
    verify_verfsigs_end = time.time()
    t_versigs = verify_verfsigs_end - verify_verfsigs_start
    pprint("   * [T] Verifying %s verifier BBS+ signatures:" % len(mysigs), fmt(t_versigs), "s", statusstr(status))

    # Gather signatures at the root process
    gathersigs_start = time.time()
    sigs = serialize_gather(mysigs)
    gathersigs_end = time.time()
    t_gathersigs = gathersigs_end - gathersigs_start
    pprint("   * [T] Gathering verifier BBS+ signatures at the root:", fmt(t_gathersigs), "s")

    return status, sigs, pk, t_map, t_gensigs, t_versigs, t_gathersigs

def pipokproof(B1, B2, y, eg1y, rho, _s1, _s2, _c, _r, _d1, _d2):    
    invB1 = B1 ** (-1)
    eB2y = pair(B2, y)
    inveB2f2 = pair(B2, invf2)
    egstar = eB2y * (invef1f2) * (inveg1f2 ** (rho))

    # Now proving PK{(s1, s2, c, r, d1, d2):
    #                  B1     = g1^s1 h1^s2             AND
    #                  iden   = invB1^c g1^d1 h1^d2     AND
    #                  egstar = inveB2f2^c eg1y^s2 eg1f2^d2 eh1f2^r }

    stmt = (f1, g1, h1, f2, invB1, inveB2f2,
            eg1y, eg1f2, eh1f2, B1, iden, egstar)

    # Commit (each ri corresponds to the respective known discrete log)
    _r1, _r2, _r3, _r4, _r5, _r6 = group.random(ZR, 6)
    com1 = (g1 ** _r1) * (h1 ** _r2)
    com2 = (invB1 ** _r3) * (g1 ** _r5) * (h1 ** _r6)
    com3 = (inveB2f2 ** _r3) * (eg1y ** _r2) * (eg1f2 ** _r6) * (eh1f2 ** _r4)

    # Challenge
    c = group.hash(stmt + (com1, com2, com3), type=ZR)

    # Response
    z1 = _r1 - c * _s1
    z2 = _r2 - c * _s2
    z3 = _r3 - c * _c
    z4 = _r4 - c * _r
    z5 = _r5 - c * _d1
    z6 = _r6 - c * _d2

    return c, (z1, z2, z3, z4, z5, z6)

def pipokverify(B1, B2, y, eg1y, rho, pf):
    invB1 = B1 ** (-1)
    eB2y = pair(B2, y)
    inveB2f2 = pair(B2, invf2)
    egstar = eB2y * (invef1f2) * (inveg1f2 ** (rho))

    # Now verifying PK{(s1, s2, c, r, d1, d2):
    #                  B1     = g^s1 h^s2             AND
    #                  iden   = invB1^c g^d1 h^d2     AND
    #                  egstar = inveB2f2^c eg1y^s2 eg1f2^d2 eh1f2^r }

    stmt = (f1, g1, h1, f2, invB1, inveB2f2,
            eg1y, eg1f2, eh1f2, B1, iden, egstar)

    c, (z1, z2, z3, z4, z5, z6) = pf

    v1 = (B1 ** c) * (g1 ** z1) * (h1 ** z2)
    v2 = (iden) * (invB1 ** z3) * (g1 ** z5) * (h1 ** z6)
    v3 = (egstar ** c) * (inveB2f2 ** z3) * \
        (eg1y ** z2) * (eg1f2 ** z6) * (eh1f2 ** z4)

    return c == group.hash(stmt + (v1, v2, v3), type=ZR)

def get_jsigs_jrands(sigs, _rands, _shift):
    if rank == 0:
        getsigsrands_start = time.time()
        n = len(sigs)
        jsigs, _jrands = [], []
        for j in range(len(sigs)):
            _i = (j - _shift) % n
            jsigs.append(sigs[_i])
            _jrands.append(_rands[_i])
        getsigsrands_end = time.time()
        t_getsigsrands = getsigsrands_end - getsigsrands_start
        pprint("   * [T:root] Gathering per-core verifier signatures:", fmt(t_getsigsrands), "s")
    else:
        jsigs, _jrands = None, None
        t_getsigsrands = 0
    return jsigs, _jrands, t_getsigsrands

def pipokproofs(vals, y, _jrands, jsigs):
    B1B2_pr = []
    eg1y = pair(g1, y)
    eg1y.initPP()
    for j in range(len(vals)):
        rho = vals[j]
        _sig = jsigs[j]
        _rand = _jrands[j]

        # NIZK proof
        _sigdoubledash = bbsplussign_randomise(_sig, _rand)
        _S, _c, _r = _sigdoubledash
        _s1, _s2 = group.random(ZR, 2)
        _d1 = _c * _s1
        _d2 = _c * _s2
        B1 = commit(_s1, _s2)
        B2 = _S * (g1 ** _s2)
        B1B2_pr.append((B1, B2, pipokproof(B1, B2, y, eg1y, rho, _s1, _s2, _c, _r, _d1, _d2)))
    return B1B2_pr

def pipokverifies(vals, B1B2_pr, y):
    status = True
    eg1y = pair(g1, y)
    eg1y.initPP()
    for j in range(len(vals)):
        B1, B2, pf = B1B2_pr[j]
        rho = vals[j]
        status = status and pipokverify(B1, B2, y, eg1y, rho, pf)
    return status

def pipokproofs_par(vals, y, _jrands, jsigs):
    map_start = time.time()
    bcast_start = time.time()
    myy = serialize_bcast(y)
    bcast_end = time.time()
    pprint("      = [T] Broadcasting pk to worker nodes", fmt(bcast_end - bcast_start), "s")
    scatter_start = time.time()
    myvals, _myjrands, myjsigs = serialize_scatter(vals, _jrands, jsigs)
    scatter_end = time.time()
    pprint("      = [T] Scattering vals, jrands, jsigs to worker nodes", fmt(scatter_end - scatter_start), "s")
    map_end = time.time()
    t_map = map_end - map_start
    pprint("   * [T] Distributing data to worker nodes:", fmt(t_map), "s")

    genproofs_start = time.time()
    myB1B2_pr = pipokproofs(myvals, myy, _myjrands, myjsigs)
    genproofs_end = time.time()
    t_genproofs = genproofs_end - genproofs_start
    pprint("   * [T] Generating %s NIZK proofs:" % len(myB1B2_pr), fmt(t_genproofs), "s")

    verproofs_start = time.time()
    status = pipokverifies(myvals, myB1B2_pr, myy)
    verproofs_end = time.time()
    t_verproofs = verproofs_end - verproofs_start
    pprint("   * [A] Verifying %s NIZK proofs:" % len(myB1B2_pr), fmt(t_verproofs), "s", statusstr(status))

    return status, t_map, t_genproofs, t_verproofs

def test_reverse_set_membership(n, mode):
    # Recording the set of plaintexts and corresponding commitments
    if rank == 0: 
        vals, t_genvals = genvals(n)
        _shift, _rands, comms, t_gencomms = gencomms(vals)
    else:
        vals, _shift, _rands, comms = None, None, None, None
        t_genvals, t_gencomms = 0.0, 0.0
    pprint(" - [T:root] Generating %s plaintext values:" % n, fmt(t_genvals), "s")
    pprint(" - [T:root] Generating %s commitments:" % n, fmt(t_gencomms), "s")

    # Prover proves knowledge of the committed values
    pprint(" - [T/A] Generating and verifying %s NIZK PoKs of committed values:" % n)
    ivals, t_getivals = get_ivals(vals, _shift)
    status, t_map, t_genproofs, t_verproofs = PoKcommproofs_par(comms, ivals, _rands)
    tea_pokcomm = t_getivals + t_map + t_genproofs
    tuv_pokcomm = t_map + t_verproofs
    pprint("   * [T/A] Total:", fmt(tea_pokcomm), "/", fmt(tuv_pokcomm), "s")

    # Verifier sends fresh public key and BBS+ signatures on all elements of vals 
    # and prover verifies each signature in batch mode
    pprint(" - [A/T] Generating and verifying %s verifier BBS+ signatures:" % n)
    status1, sigs, y, t_map, t_gensigs, t_versigs, t_gathersigs = verfsigs_commitment_par(comms, ivals, _rands)
    tuv_sigs = t_map + t_gensigs
    tea_sigs = t_map + t_versigs + t_gathersigs
    pprint("   * [A/T] Total:", fmt(tuv_sigs), "/", fmt(tea_sigs), "s")

    if mode == "basic":
        return status1

    # Generating multiple NIZKs
    pprint(" - [T/A] Generating and verifying %d NIZK proofs:" % n)
    jsigs, _jrands, t_getsigsrands = get_jsigs_jrands(sigs, _rands, _shift)
    status2, t_map, t_genproofs, t_verproofs = pipokproofs_par(vals, y, _jrands, jsigs)
    tea_pr = t_getsigsrands + t_map + t_genproofs
    tuv_pr = t_map + t_verproofs
    pprint("   * [T/A] Total:", fmt(tea_pr), "/", fmt(tuv_pr), "s")

    return status1 and status2


if __name__ == "__main__":
    n = int(sys.argv[1])
    mode = sys.argv[2] if len(sys.argv) == 3 else None
    try:
        init_strios()
        gtype = group.groupType()
        pprint("Curve: %s" % gtype)
        pprint("\nInitialising generators:", fmt(genend - genstart), "s")
        #pprint("\nSet membership:")
        #test_set_membership(n, mode)
        comm.barrier()
        pprint("\nReverse set membership:")
        test_reverse_set_membership(n, mode)
    finally:
        print_strios()
        close_strios()
