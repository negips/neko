# Appendix {#appendix}

## Environmental variable reference {#env-var}

Name                    | Description                                                   | Default value
----                    | -----------                                                   | -------------
`NEKO_AUTOTUNE`         | Force Ax auto-tuning strategy (``'1D'``,``'KSTEP'``)          | Unset
`NEKO_LOG_LEVEL`        | Log verbosity level (integer > 0, default: 1)                 | Unset
`NEKO_GS_STRTGY`        | Gather-scatter device MPI sync. strategy (0 < integer < 5 )   | Unset

### Logging level details

A number of logging levels are supported.

- `NEKO_LOG_LEVEL=0`   : Quiet mode, minimal logging during execution.
- `NEKO_LOG_LEVEL=1`   : Default information mode, adding step informations.
- `NEKO_LOG_LEVEL=2`   : Verbose mode, logging extra details.
- `NEKO_LOG_LEVEL=10`  : Debug mode.