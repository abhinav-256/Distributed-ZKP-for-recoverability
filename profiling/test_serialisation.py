from charm.core.engine.util import objectToBytes, bytesToObject
from charm.toolbox.pairinggroup import PairingGroup, ZR, G1, G2, GT, pair
from collections.abc import Iterable
import sys
import time

def serialize_wrapper(item, group):
    """ If item is a group element, call group.serialize; otherwise
    serialize using pickle. """
    if item.__class__.__name__ == "Element":
        return group.serialize(item)
    elif isinstance(item, Iterable):
        return [serialize_wrapper(iitem, group) for iitem in item]
    else:
        return item



if __name__ == "__main__":
     n = int(sys.argv[1])
     gr = PairingGroup('BN254')

     gen_start = time.time()
     g1s = [gr.random(G1) for i in range(n)]
     gen_end = time.time()
     print("generation of %s G1 elements:" % n, gen_end - gen_start, "s")

     serial2_start = time.time()
     g1s_bytes2 = serialize_wrapper(g1s, gr)
     serial2_end = time.time()
     print("my serialisation of %s G1 elements:" % n, serial2_end - serial2_start, "s")

     serial_start = time.time()
     g1s_bytes = objectToBytes(g1s, gr)
     serial_end = time.time()
     print("serialisation of %s G1 elements:" % n, serial_end - serial_start, "s")

     deserial_start = time.time()
     g1s_dash = bytesToObject(g1s_bytes, gr)
     deserial_end = time.time()
     print("deserialisation of %s G1 elements:" % n, deserial_end - deserial_start, "s")


