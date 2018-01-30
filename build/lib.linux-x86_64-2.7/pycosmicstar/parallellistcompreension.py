#!/usr/bin/env python
# *-* Coding: UTF-8 *-*

# Copyright (c) Mathieu Blondel 2009
# BSD license

# Based on ideas from Gael Varoquaux
# http://gael-varoquaux.info/blog/wp-content/uploads/2009/11/parallel_py

from __future__ import division, absolute_import

from functools import wraps


def pool_map_seq(func, iterable, chunksize=None, njobs=None):
    return list(map(func, iterable))


def pool_zipped_map_seq(func, iterable, chunksize=None, njobs=None):
    ret = []
    for args in iterable:
        if not isinstance(args, list) and not isinstance(args, tuple):
            args = (args, )
        ret.append(func(*args))
    return ret


def pool_sequentialize(iterable, njobs=None):
    ret = []
    for func, args, kw in iterable:
        ret.append(func(*args, **kw))
    return ret

try:
    import multiprocessing

    def pool_map(func, iterable, chunksize=None, njobs=None):
        """
        func must be an unary function
        """
        pool = multiprocessing.Pool(njobs)
        return pool.map(func, iterable, chunksize)

    def pool_zipped_map(func, iterable, chunksize=None, njobs=None):
        """
        func can be of variable arity and each element in iterable should
        be a tuple of the same length as func's arity
        """
        # FIXME: chunksize is currently ignored
        pool = multiprocessing.Pool(njobs)

        jobs = []
        for args in iterable:
            if not isinstance(args, list) and not isinstance(args, tuple):
                args = (args, )
            jobs.append(pool.apply_async(func, args))

        return [job.get() for job in jobs]

    def pool_parallelize(iterable, njobs=None):
        pool = multiprocessing.Pool(njobs)

        jobs = []
        for func, args, kw in iterable:
            jobs.append(pool.apply_async(func, args, kw))

        return [job.get() for job in jobs]

except ImportError:
    pool_map = pool_map_seq
    pool_zipped_map = pool_zipped_map_seq
    pool_parallelize = pool_sequentialize


def delayed(func):
    @wraps(func)
    def delayed_function(*args, **kw):
        return func, args, kw
    return delayed_function


def parallelized(func):
    @wraps(func)
    def wrapper(iterable, njobs=None):
        return pool_zipped_map(func, iterable, njobs=njobs)
    return wrapper

if __name__ == "__main__":
    from math import sqrt

    sqrtd = delayed(sqrt)
    powd = delayed(pow)

    squares = [1, 4, 9, 16]
    print((pool_parallelize([sqrtd(i) for i in squares], njobs=4)))
    print((pool_parallelize([powd(i, 0.5) for i in squares], njobs=4)))

    sqrtp = parallelized(sqrt)
    powp = parallelized(pow)
    print((sqrtp([i for i in squares], njobs=4)))
    print((powp([(i, 0.5) for i in squares], njobs=4)))
