#!/bin/bash
set -euxo pipefail

REPO_ROOT="$(realpath "$(dirname "${BASH_SOURCE[0]}")"/..)"
BUILD_DIR="$(mktemp -d)"
INSTALL_DIR="$(mktemp -d)"
DOWNSTREAM_BUILD_DIR="$(mktemp -d)"

cleanup() {
    /bin/rm -rf "$BUILD_DIR" "$INSTALL_DIR" "$DOWNSTREAM_BUILD_DIR"
}
trap cleanup EXIT TERM INT

cmake \
    -DBUILD_SHARED_LIBS=ON \
    -Denable_doubledouble=ON \
    -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR" \
    -S "$REPO_ROOT" \
    -B "$BUILD_DIR"

cmake --build "$BUILD_DIR"
ctest --test-dir "$BUILD_DIR" --output-on-failure
cmake --install "$BUILD_DIR"

cmake \
    -DCMAKE_PREFIX_PATH="$INSTALL_DIR" \
    -S "$REPO_ROOT/integration/downstreamConsumerProg" \
    -B "$DOWNSTREAM_BUILD_DIR"

cmake --build "$DOWNSTREAM_BUILD_DIR"
ctest --test-dir "$DOWNSTREAM_BUILD_DIR" --output-on-failure
