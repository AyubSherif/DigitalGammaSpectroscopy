# ADR-0002: Order captures by header timestamp, not file number

* **Status:** Accepted
* **Date:** 2026-11-24

## Context

The file numbers in `raw_data_<N>.txt` don't follow acquisition time. In the
2021 source run, consecutive captures in time carry numbers like
3293 → 3430 → 1 → 402. Files 0–99 by number are spread over 1.74 h of the 3 h
run. The numbering isn't string sort order either; its origin is unknown. Every
header carries a `Creation` timestamp with microsecond resolution.

## Decision

* The header timestamp is the authoritative order. `build_index` sorts by
  `(created, file_number)`, and stores are written in that order.
* "The first N captures" (`--limit N`) means the first N **in time**.
* The original file number is kept in the store for traceability.

## Consequences

* Time-dependent effects (gain drift, rate changes, a later live stream) are
  analysed in the correct order.
* The development subset (first 100 in time) covers a continuous 3.7 minutes.
  To sample the whole run (for example for drift studies), a future option can
  pick evenly spaced captures.
* Timestamps are naive local time. They are stored as if UTC
  (`time_s`) for arithmetic only; they are not true UTC.
