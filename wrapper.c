#include "CHOLMOD/include/cholmod.h"
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

cholmod_common * cholmod_allocate()
{
        return (cholmod_common*)calloc(1, sizeof(cholmod_common));
}

void cholmod_release(cholmod_common * ptr)
{
        free(ptr);
}

#define MAKE_ACCESSOR(FIELD, TYPE)                              \
        TYPE cholmod_get_##FIELD(const cholmod_common * common) \
        {                                                       \
                return (TYPE)(common->FIELD);                   \
        }                                                       \
                                                                \
        TYPE cholmod_set_##FIELD(cholmod_common * common, TYPE new_value) \
        {                                                               \
                TYPE old = (TYPE)(common->FIELD);                       \
                common->FIELD = (typeof(common->FIELD))new_value;       \
                return old;                                             \
        }

MAKE_ACCESSOR(print, int)
MAKE_ACCESSOR(print_function, void *)

MAKE_ACCESSOR(dbound, double)
MAKE_ACCESSOR(supernodal_switch, double)
MAKE_ACCESSOR(supernodal, int)
MAKE_ACCESSOR(selected, int)

MAKE_ACCESSOR(itype, int)
MAKE_ACCESSOR(dtype, int)

MAKE_ACCESSOR(status, int)
MAKE_ACCESSOR(fl, double)
MAKE_ACCESSOR(lnz, double)
MAKE_ACCESSOR(anz, double)
MAKE_ACCESSOR(modfl, double)
MAKE_ACCESSOR(malloc_count, size_t)
MAKE_ACCESSOR(memory_usage, size_t)
MAKE_ACCESSOR(memory_inuse, size_t)
MAKE_ACCESSOR(rowfacfl, double)
MAKE_ACCESSOR(aatfl, double)
MAKE_ACCESSOR(blas_ok, int)

#ifdef __cplusplus
}
#endif
