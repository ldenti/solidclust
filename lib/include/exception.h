#ifndef EXCEPTION_H
#define EXCEPTION_H

#define FALSE 0
#define TRUE 1

#define NDEBUG

enum exception_t {
    OK = 0,
    ERR_PARSE,
    ERR_OUT_OF_BOUNDS,
    ERR_FILE,
    ERR_LOGIC,
    ERR_RUNTIME,
    ERR_MALLOC,
    ERR_IO
};

#endif /* EXCEPTION_H */