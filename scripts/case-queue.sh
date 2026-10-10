#!/usr/bin/env bash
# Bounded independent-case queue. Worker owns its case directory and result file.
# Wait for every worker even after failure, so feedback is complete before tar.
lpbf_run_case_queue() {
    local limit="$1" worker="$2" pid failed=0
    shift 2
    [[ "$limit" =~ ^[1-9][0-9]*$ ]] || return 2
    local -a pending=()
    local selection
    for selection in "$@"; do
        "$worker" "$selection" &
        pending+=("$!")
        if [ "${#pending[@]}" -ge "$limit" ]; then
            pid="${pending[0]}"
            if wait "$pid"; then :; else failed=1; fi
            pending=("${pending[@]:1}")
        fi
    done
    for pid in "${pending[@]}"; do
        if wait "$pid"; then :; else failed=1; fi
    done
    return "$failed"
}
