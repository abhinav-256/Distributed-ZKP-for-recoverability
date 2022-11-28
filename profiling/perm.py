import random
from parutils import rank, serialize_scatter, serialize_gather

def gen_rand_perm_par(n):
    if rank == 0:
        firstn = range(n)
        pi = random.sample(firstn, len(firstn))
        reverse_pi =[0]*len(pi)
        for i in range(len(pi)):
            dest = pi[i]
            reverse_pi[dest] = i
    else:
        pi,reverse_pi = None, None
    mypi = serialize_scatter(pi)
    myreverse_pi = serialize_scatter(reverse_pi)
    return mypi, myreverse_pi

def permute_par(items, perm):
    items_perm = zip(items, perm)
    all_items_perm = serialize_gather(items_perm)
    if rank == 0:
        permuted_all_items_perm = sorted(all_items_perm, key = lambda x: x[1])
        permuted_all_items,_ = zip(*permuted_all_items_perm)
    else:
        permuted_all_items = None
    my_permuted_items = serialize_scatter(permuted_all_items)
    return my_permuted_items

if __name__ == "__main__":
    items_dest = [("a", 3), ("b", 0), ("c", 1), ("d", 2)]
    permuted_items = permute_par_old(items_dest)
    print("Permuted items (old)", permuted_items)

    items = ["a", "b", "c", "d"]
    perm = [3,0,1,2]
    permuted_items = permute_par(items, perm)
    print("Permuted items (new)", permuted_items)
