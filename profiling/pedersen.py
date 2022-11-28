from globals import g1, h1

def commit(m, _r):
    return (g1**m) * (h1**_r)


def open(c, m, _r):
    return c == (g1**m) * (h1**_r)
