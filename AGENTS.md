# Waiting for GitHub Actions push CI

After a push, use the destination branch and the full commit SHA delivered to that branch. After merging a PR, use the base branch and the resulting merge commit SHA. List matching runs with `gh run list --event push --branch <branch> --commit <sha> --json databaseId,name,status,conclusion`.

For each applicable push run, use `gh run watch <run-id> --interval 300 --exit-status` with its `databaseId`. The command refreshes every five minutes until the run finishes and returns a nonzero exit status on failure. If no matching push run exists, report that for the branch and SHA.
