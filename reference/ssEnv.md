# ssEnv

Internal session environment object persisted between SEMseeker analysis
steps. Stores runtime parameters (result folder paths, technology flag,
alpha threshold, etc.) set by
[`core_init_env`](https://corsaro-lab.github.io/SEMseeker/reference/core_init_env.md)
and retrieved by `core_get_session_info()`.

## Usage

``` r
ssEnv
```

## Format

An `environment` containing named slots for session-level analysis
parameters. Users should not modify this object directly; use
[`core_init_env`](https://corsaro-lab.github.io/SEMseeker/reference/core_init_env.md)
and `core_set_env_variable()` instead.
