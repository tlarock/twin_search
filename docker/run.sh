#!/usr/bin/env bash
# Helper to run a command inside the twin_search Linux dev container.
#
# Layout inside the container:
#   /work/twin_search   <- this repository (read-write)
#   /work/argparse     <- morrisfranken/argparse (read-only)
#   /work/discreture   <- mraggi/discreture   (read-only)
#   /build             <- out-of-tree build directories (read-write)
#
# Usage: docker/run.sh <command...>
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PARENT="$(dirname "$REPO")"
BUILDS="${TWIN_SEARCH_BUILD_DIR:-$REPO/../.twin_search_builds}"
mkdir -p "$BUILDS"

exec docker run --rm -it \
    -v "$REPO":/work/twin_search \
    -v "$PARENT/argparse":/work/argparse:ro \
    -v "$PARENT/discreture":/work/discreture:ro \
    -v "$BUILDS":/build \
    -w /work/twin_search \
    twin-search-dev:latest "$@"
