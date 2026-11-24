# Architecture Decision Records

Short records of significant design decisions: the context, the decision, and
its consequences. A record is never edited to change its meaning. A decision
that is reversed gets a new record that supersedes the old one.

| # | Title | Status |
|---|---|---|
| [0001](0001-lossless-integer-hdf5-store.md) | Store captures as lossless integer ADC codes in HDF5 | Accepted |
| [0002](0002-chronological-ordering.md) | Order captures by header timestamp, not file number | Accepted |
| [0003](0003-modular-optional-stages.md) | Processing stages are modular and individually optional | Accepted (not yet implemented) |
| [0004](0004-template-construction.md) | Pulse template built from data on a sub-sample grid with phase equalisation | Accepted |

New records: copy the structure of an existing one, take the next number.
