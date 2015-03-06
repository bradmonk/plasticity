import random


def qux():
    return random.random()


def waldo(J):
    print 'type(J): %r' % (type(J),)
    print 'J.dtype: %r' % (J.dtype,)
    print 'J.shape: %r' % (J.shape,)
    J *= 2.1
