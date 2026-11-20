#!/usr/bin/env bash
# Publish the built manual in docs/build to the gh-pages branch on origin,
# which GitHub Pages serves (Settings → Pages → Deploy from a branch →
# gh-pages, / (root)). Build the manual first:
#
#     julia --project=docs -e 'using Pkg; Pkg.instantiate(); include("docs/make.jl")'
#
# Then run docs/publish.sh from anywhere in the repository. An optional
# argument is added to the commit message. The branch is created on the first
# run, and each run replaces its contents with the new build. The working
# copy is not touched: the branch is checked out in a temporary worktree.
set -euo pipefail

root=$(git rev-parse --show-toplevel)
build="$root/docs/build"
if [ ! -f "$build/index.html" ]; then
    echo "There is no built manual in docs/build. Build it first:" >&2
    echo "    julia --project=docs -e 'using Pkg; Pkg.instantiate(); include(\"docs/make.jl\")'" >&2
    exit 1
fi

source=$(git -C "$root" rev-parse --short HEAD)
if ! git -C "$root" diff --quiet HEAD -- src docs/src docs/make.jl; then
    source="$source with uncommitted changes"
fi

tmp=$(mktemp -d)
site="$tmp/site"
cleanup() {
    git -C "$root" worktree remove --force "$site" >/dev/null 2>&1 || true
    rm -rf "$tmp"
}
trap cleanup EXIT

if git -C "$root" fetch -q origin gh-pages 2>/dev/null; then
    git -C "$root" worktree add -q -B gh-pages "$site" FETCH_HEAD
else
    git -C "$root" worktree add -q --detach "$site"
    git -C "$site" checkout -q --orphan gh-pages
fi
# On a new orphan branch every file counts as staged, so the removal is forced.
git -C "$site" rm -rqf --ignore-unmatch .
cp -R "$build"/. "$site"/
touch "$site/.nojekyll"   # serve the files as they are, without Jekyll
git -C "$site" add -A
if git -C "$site" diff --cached --quiet; then
    echo "gh-pages already holds this build."
    exit 0
fi
git -C "$site" commit -q -m "Manual built from $source" ${1:+-m "$1"}
git -C "$site" push -q origin gh-pages
echo "Published the manual built from $source to gh-pages."
