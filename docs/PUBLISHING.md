# Publish under sanatat on GitHub

The prepared archive intentionally has no `.git` history. Create a new, **empty** repository at https://github.com/new named `TrickleBed-Dynamics` (do not initialize it with a README, gitignore, or license). Pick public/private visibility yourself.

From this extracted folder, run:

```bash
cd ~/Downloads/TrickleBed-Dynamics-GitHub
# Optional: set your Git author identity if it is not configured globally
# git config --global user.name "YOUR NAME"
# git config --global user.email "YOUR COMMIT EMAIL"
git init -b main
git add .
git commit -m "Initial reproducible research release: TrickleBed-Dynamics v2.9"
git remote add origin https://github.com/sanatat/TrickleBed-Dynamics.git
git push -u origin main
```

GitHub may request authentication via browser, a credential helper, or a personal access token; a GitHub account password alone is not accepted for HTTPS Git pushes. Never paste tokens into a public chat or commit them. If the branch was already initialized, omit `git init` and use the existing remote/branch instead.

**Read the README and choose a license before making the repository public.** A public repository without a license is source-visible but not generally open source.
