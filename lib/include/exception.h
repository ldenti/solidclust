#ifndef EXCEPTION_H
#define EXCEPTION_H

#define FALSE 0
#define TRUE 1

#define NS_IN_SEC (1000 * 1000 * 1000)
#define NS_IN_USEC 1000
#define RESOLUTION NS_IN_SEC

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
