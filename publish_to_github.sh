#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
REMOTE="https://github.com/sanatat/TrickleBed-Dynamics.git"
command -v git >/dev/null || { echo "Git is required" >&2; exit 1; }
if ! git config user.name >/dev/null || ! git config user.email >/dev/null; then
  echo 'Configure Git name and email first with: git config --global user.name "Your Name" and git config --global user.email "your-email"' >&2
  exit 1
fi
if [[ ! -d .git ]]; then git init -b main; fi
if ! git remote get-url origin >/dev/null 2>&1; then
  git remote add origin "$REMOTE"
fi
if [[ "$(git remote get-url origin)" != "$REMOTE" ]]; then
  echo "Existing origin differs from $REMOTE; refusing to overwrite" >&2
  exit 1
fi
git add -A
if ! git diff --cached --quiet; then
  git commit -m "Publish TrickleBed-Dynamics v2.9 research code and results"
fi
git push -u origin main
