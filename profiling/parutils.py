from mpi4py import MPI
from collections.abc import Iterable
import io
import time

from globals import group, pai_group

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()
strios = []

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

def fullname(o):
    klass = o.__class__
    module = klass.__module__
    if module == '__builtin__':
        return klass.__name__ # avoid outputs like '__builtin__.str'
    return module + '.' + klass.__name__

def serialize_wrapper(item):
    """ If item is a group element, call group.serialize; otherwise
    serialize using pickle. """
    fn = fullname(item)
    if fn == "integer.Element":
        return (fn, pai_group.serialize(item))
    elif item.__class__.__name__ == "Element":
        return (fn, group.serialize(item))
    elif isinstance(item, str):
        return (fn, bytes(item, 'utf-8'))
    elif isinstance(item, Iterable):
        return (fn, [serialize_wrapper(iitem) for iitem in item])
    else:
        return (fn, item)

def deserialize_wrapper(sitem):
    """ Try to deserialize the serialized item as a group element.
    If it fails, deserialize using pickle. """
    fn, _sitem = sitem
    if fn == "integer.Element":
        item = pai_group.deserialize(_sitem)
    elif fn == "pairing.Element":
        item = group.deserialize(_sitem)
    elif fn == "builtins.str":
        item = str(_sitem, 'utf-8')
    elif isinstance(_sitem, Iterable):
        item = [deserialize_wrapper(siitem) for siitem in _sitem]
    else:
        item = _sitem
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
