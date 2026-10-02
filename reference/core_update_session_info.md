# Persist the session environment

Stores ssEnv into \`.pkgglobalenv\` (in-memory cache, always) and
optionally onto disk as \`session_info.rds\` and
\`\<YYYY-MM-DD\>\_session_info.rds\` inside \`ssEnv\$session_folder\`.

## Usage

``` r
core_update_session_info(ssEnv, save_to_disk = TRUE)
```

## Arguments

- ssEnv:

  list. Session environment.

- save_to_disk:

  logical. If TRUE (default, backward-compatible) writes the session to
  disk as well as in-memory. If FALSE, only updates the in-memory
  cache - fast path for worker bodies inside \`foreach %dorng%\`.

## Value

ssEnv, invisibly.

## Details

\*\*Hot-path callers (inside foreach loops) MUST pass \`save_to_disk =
FALSE\`\*\* to avoid hammering the disk with one 15 MB rds write per
gene per worker. The on-disk persistence is meant for end-of-job
snapshots, not per-iteration state syncs: the split exists because
per-iteration writes were a measured performance regression.
