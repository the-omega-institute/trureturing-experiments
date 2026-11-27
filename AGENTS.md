# Waiting for GitHub Actions CI

When waiting for a CI run, use `gh run watch <run-id> --interval 300 --exit-status` with its GitHub Actions run ID. The command refreshes every five minutes until the run finishes and returns a nonzero exit status on failure.
