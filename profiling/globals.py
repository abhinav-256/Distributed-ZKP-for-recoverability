import sys
from mpi4py import MPI
from charm.toolbox.pairinggroup import PairingGroup, G1, G2, GT, ZR, pair
from charm.toolbox.integergroup import RSAGroup

n = int(sys.argv[1])
alpha = int(sys.argv[2])

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

group = PairingGroup('BN254')
pai_group = RSAGroup()
q = group.order()
pai_p, pai_q, pai_n = pai_group.paramgen(secparam=1536)

# choose beta to be roughly (here, exactly) equal to the order 
# of the exponent group
beta = 1#group.order()

def generators():
    f1, g1, h1 = [group.random(G1) for i in range(3)]
    f2 = group.random(G2)
    fT = group.random(GT)
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
    idenT = eg1f2 ** (0)
    return (f1, g1, h1, f2, fT, ef1f2, eg1f2, eh1f2, invf1, invg1, invh1, invf2, invef1f2, inveg1f2, inveh1f2, iden, idenT)

def initgens_par():
    genserials = [group.serialize(gen) for gen in generators()] if rank == 0 else None
    genserials = comm.bcast(genserials, root=0)
    gens = [group.deserialize(genserial) for genserial in genserials]
    for gen in gens:
        gen.initPP()
    return gens

f1, g1, h1, f2, fT, ef1f2, eg1f2, eh1f2, invf1, invg1, invh1, invf2, invef1f2, inveg1f2, inveh1f2, iden, idenT = initgens_par()
