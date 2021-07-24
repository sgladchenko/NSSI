import timeit

def elapsed(function):
    def wrapper(*args, **kwargs):
        start_time = timeit.default_timer()
        res = function(*args, **kwargs)
        elapsedtime = timeit.default_timer() - start_time
        print(f"Elapsed time at {function.__name__}: {elapsedtime:.3f} sec")
        return res
    return wrapper    