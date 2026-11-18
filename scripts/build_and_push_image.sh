#!/usr/bin/env bash
# Build, test and (after confirmation) publish the pipeline's container image.
#
# Usage:  scripts/build_and_push_image.sh <version>      e.g. 2.0.0
#
# The image repository (e.g. adamd3/panscript) is read from process.container
# in nextflow.config. Images are always version-tagged; :latest is never used.
# After pushing, set process.container to the new tag in nextflow.config and
# re-capture tests/baseline (see CLAUDE.md, "Changing the software environment").
set -euo pipefail

VERSION=${1:-}
if [[ -z "$VERSION" ]]; then
    echo "Usage: $0 <version>   (e.g. 2.0.0)" >&2
    exit 1
fi

cd "$(dirname "$0")/.."
CURRENT=$(grep -E "^process\.container\s*=" nextflow.config | sed -E "s/.*'([^']+)'.*/\1/")
REPO=${CURRENT%:*}
IMAGE="$REPO:$VERSION"

if docker manifest inspect "$IMAGE" > /dev/null 2>&1; then
    echo "ERROR: $IMAGE already exists on the registry; choose a new version." >&2
    exit 1
fi

echo "== Building $IMAGE (linux/amd64)"
docker build --platform linux/amd64 -t "$IMAGE" .

echo "== Smoke test"
docker run --rm --platform linux/amd64 "$IMAGE" bash -c \
    'kallisto version && fastp --version && python -c "import pandas, umap, matplotlib" && Rscript -e "library(DESeq2)" > /dev/null'

echo "== Unit tests inside the image"
docker run --rm --platform linux/amd64 -u "$(id -u):$(id -g)" -e HOME=/tmp \
    -v "$PWD":/repo -w /repo "$IMAGE" python -m unittest discover -s tests/unit

read -r -p "Push $IMAGE to the registry? [y/N] " answer
if [[ "$answer" =~ ^[Yy]$ ]]; then
    docker push "$IMAGE"
    echo
    echo "Pushed $IMAGE. Next:"
    echo "  1. set process.container = '$IMAGE' in nextflow.config"
    echo "  2. run the test profile and re-capture tests/baseline (tests/README.md)"
else
    echo "Not pushed. The image is available locally as $IMAGE."
fi
