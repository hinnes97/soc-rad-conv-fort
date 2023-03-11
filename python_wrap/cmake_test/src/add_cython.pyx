cdef extern from "add_cython.h":
    cdef void c_add(int *a, int *b, int *c)

def add(int a, int b):
    cdef int c
    c_add(&a, &b, &c)
    return c
