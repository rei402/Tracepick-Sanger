/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "thal_parameters.h"

static unsigned char oligo[26];

unsigned char *input_address(void) { return oligo; }

int prepare(void) {
    thal_parameters params = {0};
    thal_results result = {0};
    set_default_thal_parameters(&params);
    int status = get_thermodynamic_values(&params, &result);
    thal_free_parameters(&params);
    return status;
}

double structure_tm(int hairpin) {
    thal_args args;
    thal_results result = {0};
    set_thal_oligo_default_args(&args);
    args.type = hairpin ? thal_hairpin : thal_any;
    thal(oligo, oligo, &args, THL_FAST, &result, 0);
    if (result.msg[0] != '\0') __builtin_trap();
    return result.no_structure ? 0.0 : result.temp;
}
