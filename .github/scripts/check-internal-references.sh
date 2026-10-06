#!/usr/bin/env bash
# Refuses references to the private project-management tooling in this public
# repository: backlog item identifiers and the paths of the files that hold
# them. A comment in the code is read by everyone who clones the repository,
# and an identifier is a pointer nobody outside the project can follow.
#
#   check-internal-references.sh            whole tracked tree (CI)
#   check-internal-references.sh --staged   lines added by the staged commit
#                                           (pre-commit hook)
#
# History is not scanned: cleaning the tree does not rewrite it, and this
# check exists to keep new references out, not to audit old commits.
set -euo pipefail

pattern='\b(?:AI|BKL)-[0-9]{3}\b|sestante|backlog\.(?:done\.)?json|\.coding_agent|daily_focus'
self='.github/scripts/check-internal-references.sh'
# The ignore lists name the agent working folder precisely to keep it out of
# the repository and the package, which is the opposite of a leak.
excludes=(":(exclude)$self" ":(exclude).gitignore" ":(exclude).Rbuildignore")

if [[ "${1:-}" == "--staged" ]]; then
  hits=$(git diff --cached -U0 --no-color -- . "${excludes[@]}" \
    | perl -ne 'if (/^\+\+\+ b\/(.*)/) { $f = $1; next } print "$f: $_" if /^\+/ && /'"$pattern"'/')
else
  hits=$(git grep -n -I -P "$pattern" -- . "${excludes[@]}" || true)
fi

if [[ -n "$hits" ]]; then
  echo "Internal project-management references found:"
  echo "$hits"
  echo
  echo "Backlog identifiers and tooling paths stay out of this public repository."
  echo "Describe the change in its own terms instead of pointing to a tracker."
  exit 1
fi
echo "No internal project-management references."
