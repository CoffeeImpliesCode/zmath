# Run the complete unit suite over cross-test's version-pinned Tier 1/2 matrix.
# This is correctness coverage, not corpus replay or profiling.

cross_artifact() {
    local suffix=
    case "$1" in
        *-windows*) suffix=.exe ;;
        wasm32-*|wasm64-*) suffix=.wasm ;;
    esac
    printf '%s/%s/bin/zmath-test%s\n' "$CROSS_PREFIX" "$1" "$suffix"
}

cross_build_command() {
    local target=$1
    local cpu=baseline
    case "$target" in
        x86-*|x86_64-*) cpu=baseline+fma ;;
    esac
    CROSS_COMMAND=(
        "${ZIG_COMMAND[@]}" build test-build
        "-Dtarget=$target" "-Dcpu=$cpu" "-j$CROSS_BUILD_JOBS"
        --prefix "$CROSS_PREFIX/$target"
        --cache-dir "$CROSS_CACHE_DIR/$target"
        --summary all --color off
        "${CROSS_BUILD_ARGS[@]}"
    )
}
