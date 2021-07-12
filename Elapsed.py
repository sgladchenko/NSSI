#!/usr/bin/env python3

import timeit

def elapsed(function):
    def wrapper(*args, **kwargs):
        start_time = timeit.default_timer()
        res = function(*args, **kwargs)
        elapsedtime = timeit.default_timer() - start_time
        print(f"Elapsed time at {function.__name__}: {elapsedtime:.3f} sec")
        return res
    return wrapper

@elapsed
def dummyFunction(N):
    return sum(i**2 for i in range(N))

class Foo:
    def __init__(self, power):
        self.power = power

    @elapsed
    def powersum(self, N):
        return sum(i**self.power for i in range(N))

if __name__ == "__main__":
    foo = Foo(2)
    pw = foo.powersum(100000)
    dF = dummyFunction(100000)
    print(f"{pw=}, {dF=}")