# fastq-purge — Rust rewrite design

Status: DRAFT for review · Profiled: 2026-11-04 · Scope: full rewrite of the `fastq-purge` package as a standalone Rust binary.

## 1. Motivation & measured baseline

Current state (Python 3.11, `src/fastq_purge/main.py`, 869 lines, single file):

| Phase (5M-read synthetic set, 20% engineered duplicates) | Time | Per record |
|---|---|---|
| exact: build set from undetermined R1 (`build_undetermined_set`) | 2.93 s | 0.59 µs |
| exact: scan assigned files (`pool_task`) | 3.96 s | 0.79 µs |
| exact: re-read undetermined, write purged (`write_purged_fastq`) | 2.21 s | 0.44 µs |
| **exact total (single process)** | **9.1 s** | ~1.8 µs |
| approx: build bloom from undetermined (`build_bloom_filter`) | 7.96 s | 1.59 µs |
| approx: filter + write (`process_target_file`, gzip L7) | 30.46 s | 6.09 µs |
| **approx total** | **38.4 s** | ~7.7 µs |

Memory: Python set of name-strings ≈ 90–120 B/name → **~9–12 GB at 100 M names** (the practical RAM wall for NovaSeqX-scale undetermined files).

Microbench (per record): `sha256(pickle.dumps(name))` = 0.53 µs vs builtin `hash()` = 0.022 µs; rbloom `add`/`contains` with the current sha256(pickle) hash ≈ 1.1 µs each — the hash function, not the set logic, dominates the approx path.

### Problems beyond speed

1. **The released `--method approx` CLI path is broken.** `main()` passes the lane→pair **dict** to `build_bloom_filter()`, which iterates it as if it were a list of files → `AttributeError: 'str' object has no attribute 'name'` (reproduced 2026-11-04). Approx has only ever worked via direct function calls.
2. **The two methods do opposite things.** exact purges the **undetermined** file (removes undetermined reads whose names appear in assigned files — the documented intent). approx purges the **assigned** file (removes assigned reads whose full names appear in the undetermined bloom). Output files, counts, and semantics differ.
3. **The comparison key differs per method**: exact uses `name.split(" ")[0]` (name0 = instrument:run:lane:tile:x:y), approx uses the **full** `record.name` (including the index sequence). On real Illumina data name0 is globally unique per cluster, so the difference is usually small — but on data with duplicate name0s (multi-run concatenations, some platforms) the methods diverge.
4. **Fragile concurrency**: `mp.set_start_method("fork", force=True)`, three interchangeable threading backends (`loky`/`mp_pool`/`mp_manager`), global module variables shared across forked workers, results shipped as pickled sets over process boundaries, 1 M-record intersection-update batching.
5. **Three passes over the data per undetermined lane**, and in multi-lane runs the assigned files are re-scanned **once per undetermined lane** (main loop: `for un_lane in undetermined_path: build set → scan ALL assigned → write`).
6. Bloom sized for 1 B items by default (`--max-items 1e9`): a 2.4 GB bit array even when the run has 10 M reads.

### What a rewrite buys

- 1 pass structure per concern, 1 concurrency model (threads), 0 pickling of sets.
- ~5–10× less memory (u64 hashes vs strings; ~2 GB at 100 M names vs ~10–12 GB).
- 2–5× wall time at production scale (I/O-bound floor is gzip; per-record CPU is 10–50× faster).
- Both methods work, do the same thing, with documented (or explicitly flagged) differences.
- Single static binary; no Python/deps to ship.

## 2. Goals / non-goals

**Goals**

- CLI drop-in replacement (same flags, same file-discovery rules, same output naming; see §7 checklist).
- exact + approx membership, both purging **undetermined** (canonical direction, §3).
- Correctness-first: differential-tested against the current Python exact path.
- Memory bounded and documented (§6).

**Non-goals**

- Python bindings (this is a standalone tool, not a library).
- Changing the on-disk FASTQ/BCL formats.
- Parallel decompression *within* one gzip stream (across files/lanes is enough; revisit if profiling shows a wall).
- Index-file-aware purging (R1/R2/I1/I2 collapsing stays as a file-selection rule only).

## 3. Canonical semantics (decision points)

The rewrite must decide the three divergences above. Recommendation:

| # | Decision | Recommended | Rationale |
|---|---|---|---|
| D1 | Purge direction | **Purge undetermined** in both methods (remove undetermined reads whose name appears in any assigned file) | Matches README intent and the exact method's behaviour; approx's current direction is an accident of the broken path. |
| D2 | Comparison key | **name0** (`record.name` up to first space) in both methods | name0 identifies the cluster (lane/tile/x/y); index sequences are not part of identity. exact already uses it. |
| D3 | Preserve approx's old behaviour? | Optional `--purge-assigned` flag (default off) that inverts D1, for anyone relying on it | Cheap to implement (same pipeline, swap roles); keeps the door open. |

Each decision must be called out in the changelog; D1/D2 change observable output for edge-case inputs.

**D4** (recorded in the §8 M4 entry, decided 2026-11-06 with the M4 implementation): the §6 bloom sizing "k = 2 at 19.17 bits/item" is corrected to the optimal `k = round(m/n·ln 2)` at the same (rbloom-parity) bit size. Not a semantic decision — by the M4 output-identity proof the approx purged outputs remain byte-identical to exact for any `--fpr`; D4 only fixes the structural FPR the data structure actually delivers.

**D10** (legacy mode; documented as M3 deviation (6); found by the 2026-11-10 rollout
re-verification): the master Python's legacy mode accumulates the undetermined set across
lanes via a module global — `undetermined_set = set()` at `main.py` L37,
`build_undetermined_set` adds each read's name0 to it and returns the global, and main's
per-lane call rebinds its local name to the same object. The per-lane `Number of
undetermined reads` and `Found {N} duplicates in the undetermined file.` lines are therefore
cumulative across lanes (rollout run: 10000/20000/30000/40000 over four lanes). Same bug
family as D9 (the sheet-mode branch). The Rust port keeps fresh per-lane sets (M3
semantics) and does not reproduce the cumulative count lines. The cumulative set cannot
change the purged outputs: a record purged from a lane's file carries a name from that
file, hence in that lane's own set by construction, so the removal predicate
(`name0 not in assigned ∩ S`) reduces to `name0 not in assigned` under both semantics —
the accumulation affects the reported counts only (both rollout runs confirm identical
outputs; the degenerate run was a 100% self-match, Kept 0 on both sides).

## 4. Architecture

```
fastq-purge/
└── (Rust workspace)
    ├── Cargo.toml
    └── src/
        ├── main.rs        # entry point: parse → validate → dispatch
        ├── cli.rs         # clap definitions + validation (§7 parity rules)
        ├── fastq/
        │   ├── mod.rs
        │   ├── read.rs    # streaming FASTQ(.gz) reader: name0 + raw record
        │   └── write.rs   # gzip writer (level 7), single- and paired-stream
        ├── membership/
        │   ├── mod.rs
        │   ├── exact.rs   # HashSet<u64> of xxhash64(name0)
        │   └── bloom.rs   # bitvec + xxh128 (two 64-bit halves)
        └── pipeline.rs    # passes, rayon orchestration, progress reporting
```

Dependencies: `clap` (derive), `flate2` (miniz_oxide backend; `read::MultiGzDecoder` for multi-member input, `write::GzEncoder` level 7 for output), `xxhash-rust` (feature `xxh64` — the crate has **no default features**), later `rayon`, `bitvec`; `anyhow` for the bin. M1 uses a ~70-line minimal level-gated logger instead of `tracing` (added only if M3/M4 need structured logging).

Layout note (M2): the crate is a **lib** (`src/lib.rs` exposing `cli`, `discovery`, `fastq`, `logging`, `membership`, `pipeline`) plus two bins — `fastq-purge-rs` (the tool) and `fastq-dump` (dev/verification: `names`/`records` dumps for differential testing, kept for M6) — and `examples/inflate_bench.rs` (dev-only inflate throughput probe).

### Key data structures

```rust
// fastq/read.rs
pub struct NameCursor { /* position of next record in the decompressed stream */ }

/// Extract name0 (bytes up to the first b' ') from the header line,
/// then advance to the start of the next record.
pub trait NameReader {
    fn next_name0(&mut self) -> Option<Result<&[u8], FastqError>>;
}

// write.rs: paired writer that consumes (r1_bytes, r2_bytes) in lockstep.

// membership/exact.rs
pub struct ExactSet { inner: HashSet<u64> }
impl ExactSet {
    pub fn insert_name(&mut self, name0: &[u8]) {
        self.inner.insert(xxhash_h64(name0, 0));
    }
    pub fn contains_name(&self, name0: &[u8]) -> bool {
        self.inner.contains(&xxhash_h64(name0, 0))
    }
}

// membership/bloom.rs
pub struct Bloom64 {
    bits: BitVec,          // size per §6
    k: u32,                // 2 (one per 64-bit half of xxh128)
}
impl Bloom64 {
    pub fn new(n_items: u64, fpr: f64) -> Self;  // m, k per §6 formulas
    pub fn insert_name(&mut self, name0: &[u8]);  // xxh128 → two bit probes
    pub fn contains_name(&self, name0: &[u8]) -> bool;
}
```

Name extraction detail (the hot primitive): read the decompressed stream in 256 KiB chunks; find the next `\n` (header end), then the first ` ` in the header (name0 end — if absent, name0 = whole header); hash the slice; then skip 3 more `\n`s (FASTQ is strictly 4 lines/record). No per-record allocation, no full-record parse for membership passes. For the write pass, buffer the 4 raw lines as one owned `Vec<u8>` per record (only for kept records).

## 5. Pass structure

Canonical pipeline (D1/D2):

```
pass 1 — build membership          (one thread per undetermined lane, rayon)
  for each undetermined lane L:
      stream undetermined R1 (L):  name0 → insert into set/bloom for L
  (multi-lane: K independent sets, one per undetermined lane — the current
   Python re-scans assigned K times; we scan once and test against all K)

pass 2 — find duplicates           (one thread per assigned file, rayon)
  for each assigned file A (all lanes, R1-only after pair/index collapse):
      stream A: for each name0: for each undetermined lane L with a set:
          if set_L.contains(name0): dups_L.add(name0)     // or count only

pass 3 — write purged undetermined (one thread per undetermined lane, rayon)
  for each undetermined lane L:
      if dups_L empty → log "Found no duplicates..." ; write NOTHING (parity!)
      else stream (R1, R2) pair:
          if !set_of(dups_L).contains(name0): write r1, r2
```

Notes:

- Pass 2 only needs the **set of duplicate names** (not counts) — `dups_L` is a `HashSet<u64>` of hash values (same memory as pass-1 sets, typically far smaller: duplicates ≪ undetermined).
- The Python exact path ships the full duplicate set back from workers and re-tests undetermined against it; we test the same way but the "set" is in-process u64s.
- approx replaces `ExactSet` with `Bloom64` in pass 1; pass 2 membership tests become bloom queries; pass 3 tests `!bloom.contains(name0)`. False positives remove reads that should be kept (documented; `--fpr` controls it).
- Threading: `--threads N` = rayon pool size (file-level parallelism). No fork, no managers, no global state. Gzip decode cost per file is the I/O floor; 4 lanes × 2 files already fills 8 cores at typical N.
- Single-undetermined-lane + single-assigned-file (the common small case) still benefits: no process startup, no pickle, tighter inner loop.

## 6. Hash & memory math

### Hash choice: xxhash

- `xxh64` for name0 identity: ~GB/s/core, 8 B output.
- `xxh128` for the bloom: two independent 64-bit halves → k = 2 probe functions (classic optimal-k for m/n ≈ 19 bits/item is ~13; double-hashing with k=2 at 19 bits/item gives fpr ≈ (e^−19/... ) — **use the standard result**: with k=2 and m = 19.17·n bits, fpr ≈ 2^−(m/n)·... in practice compute empirically and size m with 20% margin; the test plan in §8 includes an empirical FPR check).
- Why not sha256: 30–50× slower than xxhash and provides no benefit for non-adversarial cluster names (current code uses it only because `rbloom` needs *an* int hash and `pickle` was chosen for generality — it costs 0.53 µs/record).

### 64-bit collision analysis (exact mode)

Expected number of colliding name pairs: `n(n−1) / 2^65`

| n (names) | 10 M | 100 M | 1 B |
|---|---|---|---|
| expected colliding pairs | 7×10⁻⁸ | 2.7×10⁻⁴ | 2.7×10⁻² |

At 100 M names there is a 0.027% chance of *any* spurious duplicate in a run. Acceptable for "exact" in practice; if a user needs cryptographic identity, add an optional `--strict` mode that stores name0 strings (memory per §7 table) — do not build it speculatively.

### Memory table (membership, pass-1 structure)

| Structure | Per name | 10 M | 100 M | 1 B |
|---|---|---|---|---|
| Python `set[str]` (current exact) | ~90–120 B | ~1 GB | **~10–12 GB** | infeasible |
| Python `rbloom` @ fpr 1e-4 (current approx, default 1e9 items) | 19.17 bits fixed | 24 MB | 240 MB | 2.4 GB |
| Rust `HashSet<u64>` (hashbrown, load ~0.8) | ~16 B | 0.13 GB | **~2 GB** | ~16 GB |
| Rust `Bloom64` @ fpr 1e-4 | 19.17 bits | 24 MB | **240 MB** | 2.4 GB |

Bloom sizing formula (used in `Bloom64::new`): `m = ⌈n · (−ln fpr) / (ln 2)²⌉` bits → 19.17 bits/item at fpr = 1e-4; `k = 2` (xxh128 halves) with +20% bit margin, or derive k = round(m/n · ln 2) if we switch to more probe functions.

Pass-3 `dups_L` sets: duplicates are a small fraction of undetermined in practice (that's the point of the tool); budget as a fraction of the pass-1 set, monitor in tests.

## 7. CLI parity checklist

Binary name: `fastq-purge` (or `fastq-purge-rs` during development). Flags mirror `parse_args()` (`main.py:42`) exactly:

| Flag | Type | Default | Parity requirement |
|---|---|---|---|
| `--undetermined-path` | file or dir, required | — | dir → **non-recursive** glob `*.fq*` / `*.fastq*`; files whose names contain `purged` are excluded (re-run safety) |
| `--assigned-path` | file/dir, `nargs+`, required | — | dir → **recursive** glob; dedupe; collapse paired/index files: if `re.sub(r"_[IR][0-9]_", "_", name)` was seen, skip with warning `Ignoring '<f>' as it is either a paired or index file`; lane key from `re.search(r"_L[0-9]{3}_", name)` → `L001`; files without it: warning + skip |
| `--output-path` | dir | parent of first undetermined file | create if missing (`mkdir -p` semantics) |
| `--method` | `exact` \| `approx` | `exact` | both must work (approx is currently broken) |
| `--max-items` | int | `1_000_000_000` | approx only; used for bloom sizing — **recommend keeping the default** for parity but logging the actual item count observed |
| `--fpr` | float | `1e-4` | approx only |
| `--threading-method` | `loky` \| `mp_pool` \| `mp_manager` | `mp_pool` | **accept and ignore**; one-time "accepted and ignored" warning in the pipeline phase (NOT in the validation block, to keep the validation log byte-identical to Python) |
| `--threads` | int | `1` | rayon pool size |
| `--log-level` | DEBUG…CRITICAL | `INFO` | same levels; same message content and order on **stdout** (see logging note below) |
| `--version` | | | print version |

Behavioural checklist (from `validate_args` + `main`):

- [ ] Undetermined pairing: key = 3rd `_`-separated token of the filename (e.g. `L001` in `Undetermined_S0_L001_R1_001.fastq.gz`); files sorted; a lane with a single file → `(file, None)` and only that file is purged/written.
- [ ] Output naming: `X.fastq.gz` → `X.purged.fastq.gz`; `X.fastq` → `X.purged.fastq`; written into `--output-path`.
- [ ] exact logs (same strings, same order): `Using exact method`, `Building python set from undetermined reads...`, `Number of undetermined reads: {n} ({mb:.2f} MB)`, `Found {n} duplicates in the undetermined file.`, and when zero duplicates: `Found no duplicates in the undetermined reads. Exiting...` with **no output file created**.
- [ ] approx logs: `Using bloom filter method`, `Building bloom filter...`, per-file `Processed {n} reads` every 1 M records, `Removed {n} reads`, `Kept {n} reads`, `Purged fastq file saved to '{path}'`.
- [ ] Output gzip compression level 7.
- [ ] Invalid path (neither file nor dir) → error log + exit code 1.
- [ ] `--log-level` gates DEBUG-only lines (memory usage, bloom size, item counts).
- [ ] Final `Done!` line.

**Logging contract (verified against the running Python CLI, 2026-11-05)**

The Python CLI configures `RichHandler(rich_tracebacks=True, tracebacks_show_locals=True, markup=True)` — despite `format="%(message)s"` in `basicConfig`, RichHandler renders each record itself: **stdout**, layout `HH:MM:SS LEVEL message  main.py:LINE`, level column fixed width, rich wrapping at console width, timestamp omitted when the second is unchanged, colour when TTY. The Rust CLI reproduces the *message* contract — identical messages, identical order, identical `--log-level` gating, on stdout — as plain bare lines. Deliberate deviation: no rich rendering (wall-clock timestamps are non-deterministic, `main.py:LINE` references Python source that does not exist in the rewrite, wrapping depends on console width). "Path … is neither a file nor a directory" is logged to stdout once before exit 1 (parity); errors with no Python equivalent (clean error for an empty undetermined dir, M1 pipeline bail) go to stderr.

New (documented) behaviours:

- [ ] `--purge-assigned` (D3): inverts direction; off by default.
- [ ] approx + exact produce identical outputs except for bloom false positives (test-asserted on collision-free datasets).
- [ ] Multi-undetermined-lane runs scan assigned files **once** (perf note, not an observable change except speed).

## 8. Test & rollout plan

### Milestone status

- [x] **M1 — CLI + discovery/validation (done 2026-11-05).** Rust workspace at `rust/fastq-purge/` (bin `fastq-purge-rs`, v0.3.0): clap CLI with all flags/defaults/choices; `validate_args` port (undetermined pairing, output naming/dir creation, assigned cleaning + lane grouping); fnmatch port for `pathlib.glob` single-component patterns. 17 unit tests (`cargo test`). Differential validation-log parity vs the Python CLI on the 5 M dataset: **byte-identical message blocks** for (a) single-lane gz run, (b) two-lane run exposing the stale-`lane` log quirk, (c) ungzipped `.fastq` single-read lane (`'None'` lines), (d) re-run with `.purged.*` outputs present (excluded); `--version` rc 0; usage errors rc 2 (both); invalid path → message on stdout + rc 1 (both, no duplicate print); DEBUG-only output-dir line gated identically. Deviations: empty undetermined dir → clean stderr error + rc 1 (Python: IndexError traceback); `--version` prints `fastq-purge-rs 0.3.0` (name/version differ by design); plain log rendering (see §7 note).
  - Findings: the `Assigned files:` per-file log reuses the last undetermined-lane loop variable (main.py:299) — reproduced for parity. `_R1_`/`_R2_` normalisation replaces each segment with the *literal* text `_R[12]_`, which is then re-globbed as an fnmatch class. `re.sub` consumes each matched `_[IR][0-9]_`/`_R[12]_` segment including its trailing underscore, so adjacent segments sharing an underscore are NOT both collapsed/normalised (`_R1_R2_` → `_R2_`/`R[12]_R2_`). The glob patterns `*.fq*` and `*.fastq*` are disjoint (`.fastq` does not contain the substring `.fq`); a name like `x.fq.fastq.gz` matches both and appears twice in `explode_path` output (no dedup — dedup happens downstream). Rust 1.99 removed `impl Display for PathBuf` (use `.display()`).

- [x] **M2 — FASTQ streaming + name0 hashing (done 2026-11-06).** `fastq/read.rs`: line-based port of dnaio 1.2.3's C parser over a growing 256 KiB-chunked buffer (no per-record allocation on the name0 path); `fastq/write.rs`: gzip (level 7) / plain writers, single + paired; `membership/exact.rs`: `ExactSet` = `HashSet<u64>` of xxhash64(name0); dev bin `fastq-dump` (`names`/`records`) for differential testing; crate restructured to lib + bins. 50 tests (17 M1 + 4 membership + 29 reader/writer integration, the full dnaio error/normalisation matrix). Differential on the 5 M dataset (5,000,000 records): name0 dump **byte-identical** to `record.name.split(" ")[0]` via dnaio; full record copy **byte-identical** plain and after gunzip. Single-threaded name pass: 3.6 s Rust (miniz_oxide inflate ≈ 628 MB/s, on par with C zlib ≈ 622 MB/s) vs 2.0 s Python *with pigz multi-threaded inflate* (xopen auto-detects `/usr/bin/pigz`) — the gap is thread count, to be closed in M5 (rayon).
  - Findings (dnaio 1.2.3, all verified empirically; the port reproduces them exactly):
    - A record is strictly 4 lines; header marker `@` **or** `>` (a `>`-header file is mis-parsed by dnaio itself — two-header misdetection, `qualities=None` — not reproduced; Illumina data never uses `>`).
    - CRLF is normalised to LF on read; the `+` line's content (a repeated name) is **dropped** — the writer always emits a bare `+`. `name0` = header minus the leading marker up to the first space byte (0x20); tabs are not separators; no space ⇒ whole name.
    - Error model (exact strings, 1-based physical line numbers): the header `@` check is **deferred until the sequence line has been read** (bad header + EOF ⇒ premature-EOF, not the header error); the `+` check fires when the `+` line is read; `Length of sequence and qualities differ` fires after the qual line at the qual line number; non-ASCII bytes abort as soon as the offending line is read (`Error in FASTQ file at unknown line: Non-ASCII characters found in record.`, fires before truncation errors); premature-EOF quotes the raw bytes of the incomplete record attempt in Python-`repr` style (single quotes, double if the content has a `'` and no `"`).
    - Empty file ⇒ 0 records; a final record without a trailing newline is fine; a trailing blank line is an error (it starts a new record attempt).
    - gzip: **multi-member aware** (Python's gzip module reads concatenated members; flate2's plain `GzDecoder` stops at the first member, so the reader uses `flate2::read::MultiGzDecoder` — which, like `gunzip`, errors on trailing non-gzip data). Python's gzip writer emits mtime=0, but zlib vs miniz-oxide at the same level produce different compressed bytes ⇒ `.gz` outputs are compared after decompression; plain `.fastq` outputs are byte-identical (verified: 317.7 MB vs 312.6 MB raw).
    - **Bug caught by the prefix bisect:** a buffer re-pack (memmove to front) mid-record invalidated the header/sequence line spans, so the deferred `@` check read a stale offset (one buffer position ahead) — a silent name0 corruption or a spurious `expected to start with '@'` error, manifesting only for certain file-size alignments (repro: any prefix ≥ 3000 records of the 5 M plain file). Fix: re-packing is permitted only at record boundaries (`read_line(allow_repack)`; only the header read may re-pack).
    - Parser overhead ≈ 1.1 s for 1.8 GB (≈ 1.6 GB/s): the ASCII check is a second full scan of every line and `rec` (the error-quote buffer) is appended per line; fusing the scans / making the quote lazy are the obvious M5 micro-opts if they show up in profiles.

- [x] **M3 — exact pipeline, single threaded (done 2026-11-06).** `pipeline/mod.rs` (dispatch; one-time compatibility warning when `--threading-method` ≠ `mp-pool` or `--threads` ≠ 1); `pipeline/exact.rs` — the canonical pipeline per §5: pass 1 builds one `ExactSet` per undetermined lane from the R1 name0s; pass 2 streams the assigned files **once** and tests each name0 hash against all K lane sets (Python re-scans the assigned files K times); pass 3 writes the purged undetermined files (paired = zip, stop at shorter, R1 name0 gates the pair). `--purge-assigned` (D3) inverts the direction: assigned reads whose name0 appears in the lane's undetermined set are removed, per-file outputs, warning for assigned lanes without an undetermined counterpart. `membership/exact.rs` gained `hash_name0`, `insert_hash`/`contains_hash` (one hash per name0, tested against all sets), `estimated_bytes`. 5 e2e tests (`tests/pipeline.rs`): paired lane, no-duplicates (no output files + exact log), `--purge-assigned`, single-read lane, gzipped end-to-end. 55 tests total (21 unit + 29 fastq + 5 pipeline).
  - Differential (5 M + 1.1 M mixed datasets): counts **identical** to Python (`Number of undetermined reads`, `Found {N} duplicates`, `Removed`/`Kept`); on the 1.1 M dataset (1 M real records + 100 k synthetic reads whose name0s are absent from the assigned files) the purged R1 and R2 outputs are **byte-identical after decompression** (100 k records kept each side). Wall time on the 5 M set: Rust 17.3 s vs Python 19.2 s single-threaded equivalent (Python's head start is pigz multi-threaded inflate; M5 closes the rest).
  - New log contract (documented in `pipeline/exact.rs`): the Python sequence is per-lane (build → scan *all* assigned → write, ×K); the Rust sequence builds all K sets, scans once, then writes per lane. Same messages, same counts, new order (the pass-2 line set is replaced — deviation 5 below). Deliberate deviations from the Python log: (1) no `memory_usage` DEBUG lines (process/OS memory is not the same quantity and not useful for parity); (2) the `Number of undetermined reads` size figure is `set.len() × 16 B` (hashbrown entry estimate) vs Python's `sys.getsizeof` of a Python `set` — different quantity, same shape; (3) when a lane has 0 duplicates **both** `Found 0 duplicates in the undetermined file.` and `Found no duplicates in the undetermined reads. Exiting...` are printed (matching Python, which prints both too); (4) the `--threading-method`/`--threads` warning wording is new (see §7); (5) pass 2: Rust logs `Scanning assigned files for already-assigned reads...` + per file `Processing target file '{f}'` + `'{f}': {N} new duplicates` (deterministic, §8/M5 log rule); the Python exact logs process-pool internals instead (`Using multiprocessing pool...`, `Purging already assigned reads...`, DEBUG `Creating {N} child processes for lane {lane}`, per file `Process {pid} ('{f}') finished yielding {K} already assigned reads` and `Process finished: {f} with {K} already assigned reads`) — the pool mechanics are not ported (Rust scans on threads; sheet mode deliberately reproduces the branch's pool lines, §8/M7.3). (6) **D10** — the master Python's legacy mode accumulates the undetermined set across lanes via a module global (`main.py` L37 `undetermined_set = set()`; `build_undetermined_set` adds to it and returns the global; main's per-lane binding rebinds to the same object — same bug family as D9): the per-lane `Number of undetermined reads` and `Found {N} duplicates` lines are cumulative (rollout: 10000/20000/30000/40000 over four lanes). The Rust port keeps fresh per-lane sets (M3 semantics) and does not reproduce the cumulative count/Found lines; the cumulative set cannot change the purged outputs (a record purged from a lane's file carries a name from that file, hence in that lane's own set by construction) — see the D10 decision in §3. (7) pass 3: after the `Writing purged fastq files` line, Rust logs `    Removed {R} reads`, `    Kept {W} reads`, `    Purged fastq file saved to '{path}'` (per written pair); the master Python's exact `write_purged_fastq` logs only the `Writing purged fastq files` line — the three writer lines are the master's approx-path `process_target_file` lines (the §7 approx checklist), kept in Rust because exact and approx share one pass-3 writer (the M4 approx contract includes them, matching the master's approx writer).
  - Findings (dnaio paired reader): `dnaio.open(r1, r2)` **validates that the two read names match** — `Reads are improperly paired. Read name '{n1}' in file 1 does not match '{n2}' in file 2.` — and raises `FileFormatError` at the first mismatched pair. Real Illumina paired files share read names (only the `1:N:0`/`2:N:0` suffix differs), so this never fires in production; the Rust pass-3 zip does **not** enforce it (silent zip on mismatched names). Identical behaviour on valid data; the deviation is documented, not implemented.

- [x] **M4 — approx (bloom) pipeline (done 2026-11-06).** `membership/bloom.rs`: `Bloom64` — hand-rolled `Vec<u64>` bit array (deviation from §4's `BitVec`: no new dependency), `m`/`k` per decision **D4** below, Kirsch–Mitzenmacher double hashing `probe_i = (h1 + i·h2) mod m` from two seeded xxh64s, `insert_h`/`contains_h` precomputed-hash fast path (one hash pair per name0, tested against all K lanes). `pipeline/approx.rs`: the canonical 3-pass structure of §5 with **D4 construction** — pass 1: one `Bloom64` per undetermined lane from the R1 name0s (sized by `--max-items`/`--fpr`); pass 2: one scan of the assigned files, each name0 tested against all K lane blooms, matches collected into an **exact** `ExactSet` of duplicates per lane; pass 3: reuses the exact write path (`write_purged_undetermined` in `pipeline/exact.rs`, now `pub(crate)` with a `progress: bool` flag that emits the Python-style `    Processed {n} reads` every 1 M records). `--purge-assigned` (D3) mirrored: per-lane filtering of the assigned files against the lane's undetermined bloom, per-file outputs, per-1 M progress. `pipeline/mod.rs` dispatches `Method::Approx`.
  - **D4 — bloom sizing and the output-identity guarantee** (corrects §6): the §6 combination "k = 2 (xxh128 halves) at m = 19.17·n bits" is inconsistent — 19.17 is the *optimal-k* sizing (rbloom's formula, k ≈ 13); with k = 2 the structural FPR at 19.17 bits/item is `(1 − e^(−2n/m))²` ≈ **1 %**, not 1e-4. M4 keeps the §6 *size* (so the memory table holds: `m = ⌈1.2·n·(−ln fpr)/(ln 2)²⌉` bits, 19.17 bits/item @ fpr 1e-4, 20 % margin) but uses the **optimal** `k = max(1, round(m/n·ln 2))` (16 at 1e6/1e-4 → measured FPR ≈ 1.5×10⁻⁵, better than promised). Consequence for the pipeline: a name enters a lane's dups set only if it **actually appeared in the assigned files**, and an undetermined read is removed only if its name0 is in the dups set — therefore the purged outputs are **byte-identical to the exact method for ANY `--fpr`/`--max-items`** (a structural FP is an assigned name absent from the undetermined file and can never match a pass-3 read; it can only inflate the `Found {K}` count). The only FPR-observable output path is `--purge-assigned`, where the bloom query directly gates removal. Verified empirically below, incl. a saturated-bloom stress test.
  - New log contract (documented in `pipeline/approx.rs`): `Using bloom filter method` → `Building bloom filter...` → per lane (DEBUG `Adding '{und R1}' to bloom filter`, DEBUG `Bloom filter: {bits} bits, {k} probes, sized for {n} items`, then) `Number of undetermined reads: {N} ({MB:.2} MB)` → `Scanning assigned files for already-assigned reads...` → per file `Processing target file '{f}'` → per lane `Found {K} duplicates in the undetermined file.` → the shared pass-3 lines (incl. `Processed {n} reads`, `Found no duplicates... Exiting...` for K = 0) → `Done!`. Deliberate deviations: (1) no `Output file:` / `Reading target fastq file...` lines — those belonged to the old broken direction (filtering assigned files), which D1 inverts; (2) blooms are **per-lane** (Python built one shared bloom across all lanes) — output-identical, and per-lane gives each lane the full bit budget; (3) the DEBUG size line replaces Python's `Bloom filter size:` / `Number of items:` (rbloom internals); (4) `Number of undetermined reads` here is the **record count** (incl. duplicate name0s) and the size figure is the bloom byte count — in exact it is the unique-name count and a set estimate.
  - Differential (1.1 M mixed dataset, §M3): `--method approx --max-items 2000000` (m ≈ 46 M bits, k = 16, achieved FPR ≈ 2×10⁻⁹) → `Found 993958 duplicates` (**exactly** the exact count — zero FPs at this size), `Removed 1000000` / `Kept 100000`, purged R1+R2 **byte-identical after decompression** to both the Python-exact and Rust-exact outputs. Saturation stress `--max-items 1000 --fpr 0.5` (1,732-bit saturated bloom, k = 1): `Found 4849796 duplicates` (maximally inflated by FPs) but `Removed`/`Kept` and the purged outputs **still byte-identical** — the D4 guarantee holds at the worst case.
  - Tests: +3 e2e (`approx_matches_exact_output` asserts byte-identity vs exact on the same dataset — the D4 guarantee as a test; `approx_no_duplicates_writes_nothing`; `approx_purge_assigned_flag`). 62 tests total (25 unit incl. 4 bloom: sizing vs Python f64 reference, degenerate sizing, no false negatives, **empirical FPR ≤ `--fpr`** at 1e-4 with 100k insert + 100k non-member probes + 29 fastq + 8 pipeline).
  - Findings: §6's "xxh128 → k = 2" would have silently delivered a ~1 % structural FPR at the documented sizing (moot for outputs by the D4 proof, but wrong for the `--fpr` promise and for the `--purge-assigned` path); two seeded xxh64s are the natural probe source (≈ same cost as one xxh128). Pass-2 cost note: k = 16 modular probes per name0 per lane — measured negligible (the 1.1 M differential runs in the same wall time as exact).

- [x] **M5 — rayon threading (done 2026-11-06).** `rayon` dependency; main.rs builds the **global pool from `--threads`** (after validation) — `--threads` is now honoured, and the one-time warning in `pipeline/mod.rs` fires only for a non-default `--threading-method` (new wording: `'--threading-method {x}' accepted and ignored (the Rust implementation is thread-based; '--threads {n}' sets the pool size)`). All three passes run on the pool: pass 1 per-lane (build set/bloom), pass 2 per assigned **file**, pass 3 per lane with dups > 0 (`par_bridge`); `--purge-assigned` parallelizes per file.
  - **Deterministic-log rule** (the M5 invariant): workers never print — each returns the lines the sequential run would have emitted for its item (`pipeline::Lines`, aliased `Block`/`ScanBlock` in `pipeline/mod.rs`; on failure: the partial lines collected up to the error point plus the error). The orchestrator flushes blocks in canonical lane/file order, and pass-2 workers return per-lane **partial** `ExactSet`s that the orchestrator merges in file order (`ExactSet::merge` returns the newly-inserted count) — reproducing the sequential per-file `new duplicates` count **exactly** (within-file duplicates dedupe to one; across files the earlier file's copy wins, matching stream-order insert-returns-bool). Consequence: the log is **byte-identical to the M3/M4 sequential run for any `--threads`**, including on failure.
  - Reader micro-opts (the M2-flagged items, §8/M2): the ASCII check is **fused into `read_line`'s `\n` scan** (one pass per line instead of two; the non-ASCII error still fires at exactly the same point — line completion — so every verified error priority is preserved), and the premature-EOF quote is **lazy**: the per-record `rec` copy buffer is replaced by a `rec_start` offset, quoting `buf[rec_start..start)` (set at the record boundary, after the only point a repack may shift the buffer). All 29 dnaio-matrix tests still green; 5 M wall 17.3 s → 17.08 s.
  - Tests: +3 e2e (`two_lanes_cross_lane_dups` — pass 2 tests every assigned name0 against **all** lane sets, so a name from L001's assigned file that is in L002's undetermined set is an L002 duplicate; `parallel_matches_sequential_log_and_output` — `--threads 4` vs default produce **byte-identical stdout and purged outputs** for both methods on a two-lane paired dataset; `threading_warning_only_for_non_default_method`). 65 tests total (25 unit + 29 fastq + 11 pipeline); clippy clean.
  - Performance (release, 16-core machine): **5 M single-lane** (the profiling set): 17.08 s @ `--threads 16` — one file per pass, so no parallelism is possible and the §2 non-goal (no intra-stream parallel decompression) makes the single-stream inflate the floor. **4-lane dataset** (1.1 M und × 4 lanes + 5 M assigned × 4 lanes): **12.05 s @ `--threads 16` (334 % CPU) vs 34.58 s @ `--threads 1` → 2.87×**; parallelism = max(lanes, assigned files) per pass, capped by `--threads` (16 threads, 4 items ⇒ ≤ 4 cores busy per pass). Correctness: all 8 purged outputs (4 lanes × R1/R2) **byte-identical after decompression** across `--threads 16` / `--threads 1` / the M3 single-lane reference. Wall analysis: the §8 acceptance targets (exact ≤ 3 s, approx ≤ 5 s on the 5 M set) assume multi-lane/multi-file workloads; they are unreachable for a single-lane single-file run, where the 312 MB assigned-file inflate (~5.7 s) plus two 62 MB undetermined inflates set the floor regardless of thread count.

- [x] **M6 — differential + perf validation (done 2026-11-06).** All gate items green:
  - **Collision-free reference differential** (gate steps 1+2): dataset with 1 M records/lane, **globally unique** name0s `UND%07d` (matching R1/R2 pair names), assigned = 400 k overlapping `UND%07d` + 600 k unique `ASS%07d` (150 bp). Python exact (`--threads 1`) vs Rust exact: counts identical (`Number of undetermined reads: 1000000`, `Found 400000 duplicates`); purged R1+R2 **byte-identical after decompression** (600 k records each; the first 400 k overlap names correctly removed).
  - **Property tests** (gate step 3): +2 committed unit tests in `membership/bloom.rs` — `empirical_fpr_within_target_10m` (10 M inserts / 10 M non-member probes, measured FPR ≤ `--fpr` 1e-4) and `saturation_increases_found_count` (saturated 1,732-bit bloom ⇒ FPR ≈ 1 vs ≈ 2 FPs comfy — the monotonic-in-`--max-items` mechanism). Binary-level check on the reference dataset: `--max-items 1000000` → `Found 400008` (400000 true + 8 structural FPs); `--max-items 1000 --fpr 0.5` → `Found 1000000` (saturated); **both** `Removed 400000` / `Kept 600000`, purged outputs byte-identical to each other **and** to the exact-mode reference — the D4 guarantee at the binary level.
  - **+2 e2e edge tests**: `all_duplicated_writes_empty_purged_files` (purged file exists, 0 bytes) and `empty_undetermined_file_runs_clean` (clean run, zero outputs). 69 tests total (27 unit + 29 fastq + 13 pipeline); fmt + clippy clean.
  - **Memory acceptance — MET**: `examples/mem_100m.rs` (release): ExactSet 100 M unique name0s ⇒ **RSS 1.13 GB**; Bloom64(1e8, 1e-4) ⇒ m = 2,300,414,011 bits, k = 16 ⇒ **RSS 0.27 GB** — both far under 2.5 GB. Same-day Python `set[str]` 100 M ⇒ ≈ 13.2 GB peak RSS (~12× the Rust exact figure; matches the §1 wall estimate).
  - **Performance verdict** (same day, 5 M set, 16 cores): Python exact 20.5 s; Rust exact 16t 18.0–20.7 s (99 % CPU — the single-stream inflate floor, §8/M5 wall analysis); Rust approx 16t 19.7 s; Python approx 35.8 s (old direction — it purges the assigned files). Counts identical on all sides. The §8 ≤ 3 s / ≤ 5 s targets are **NOT met on single-lane single-file workloads** (they assume multi-lane workloads); the ≥ 2×-Python criterion is **NOT met at 5 M single-lane** (exact ≈ 1.0×, approx 1.8×) but **is met multi-lane** (4-lane 2.87×, §8/M5); the "≥ 50 M reads" clause is not evaluable — the largest available set is 5 M.
  - **Rollout**: the changelog call-out for D1–D4 + the fixed approx path is a release task on the Python package (the Rust tree has no changelog convention); ship side-by-side per §8/Rollout.
  - **Rollout re-verification (2026-11-10, real NovaSeqX flowcell, legacy mode):** full-log comparison, master Python (`--undetermined-path`/`--assigned-path`) vs Rust master, 4 lanes × 10 k undetermined reads — logs identical except the documented deviations (M3 1–7) and the live `memory_usage` values; purged outputs identical after decompression (note: pointing `--assigned-path` at the demultiplexing root makes every undetermined R1 file a member of its own assigned list in **both** CLIs — 100 % self-match, `Kept 0` — the sheet mode's `Undetermined` exclude exists to prevent exactly this). Found and fixed here: the exact pass-1 line had been `Building membership set from undetermined reads...` since M3 (the M3/M6 gates compared counts and outputs, not the full log); the M1 checklist string `Building python set from undetermined reads...` (what the Python logs) is now what Rust logs. The same run exposed **D10**: the master Python's per-lane `Number of undetermined reads` / `Found {N} duplicates` lines are cumulative (10000 → 40000 over the four lanes — the `main.py` L37 module-global accumulation, the D9 bug family), while Rust logs 10000 / Found 10000 flat (fresh per-lane sets, M3 semantics); documented as D10 (§3) + M3 deviation (6); the cumulative set changes only the reported counts and cannot change the purged outputs; outputs identical on this dataset (100% self-match). A second, **non-degenerate** run (assigned = the four recipe dirs only, so the undetermined files are not members of their own assigned list): L004 `Found 56 duplicates` / `Kept 9944` on **both** sides, L001–L003 zero-dup; purged outputs byte-identical after decompression (record-level identity on real data); the full-log diff again matches only the documented deviations — and it caught one more: the pass-3 writer lines (`Removed`/`Kept`/`Purged fastq file saved to`), which the master exact does not log (they are the master's approx-path writer lines, §7), now documented as M3 deviation (7).

- [x] **M7.1 — SampleSheet.csv parser (done 2026-11-09, `8c87e56`).** `samplesheet.rs`: parser for the alternate input mode (the `samplesheet_support` branch's `parse_sample_sheet`): FCID anchor-line detection, polars-equivalent transforms (Lane `zfill(3)`/`pad_start(4,'L')`, index1/2 character lengths concatenated into the Recipe, the ≥2-distinct-Recipes-per-lane filter, `(Lane, Sample_Project)` stable ordering with names in sheet-row order), quote-aware CSV. Differential vs the branch's parser on polars 1.30.0 (the version pinned by the branch's `uv.lock`; `group_by().all()` list-aggregates there, so the branch's dict is nested `{lane: {proj: [[names]]}}` and the branch's own selection loop flattens it — the Rust `names_flat` matches the effective semantics): identical `(lane, project, name)` structure for basic / CRLF+blank-line / quoted-comma / ragged-row sheets; the no-FCID error is byte-identical (incl. trailing space) with rc 1 on both sides; missing-column and mis-anchored sheets get a clean Rust error where the branch traceback'd in polars (M1 precedent). 9 unit tests; dev example `sheet_dump` for differentials.

- [x] **M7.2 — SampleSheet input mode: CLI + discovery (done 2026-11-09, `de3ab2f`).** `cli.rs`: `--undetermined-path`/`--assigned-path` become optional; new `--flowcell-path`, `--sample-sheet`, `--keep-original` with the branch's help texts (the three shared flags keep their M1/master help texts — the branch's differ slightly, e.g. "The search is recursive."). **D8 — three input modes (superset):** (a) *flowcell* — the sheet is auto-discovered as the last flowcell entry matching the branch's unescaped-dot `SampleSheet.(csv|txt)$` (a clean error where the branch's `IndexError` would fire); undetermined/assigned default to the flowcell dir; (b) *manual* — `--sample-sheet` without `--flowcell-path`: **D5** fixes the branch's manual-mode crash (`validate_args` does `fc.is_dir()` on `fc = None` → `AttributeError`, so the intended else-branch is unreachable in the branch); the port implements the intended behaviour; (c) *legacy* — no flowcell, no sheet: the unchanged M1 directory-discovery flow, with the rc 2 clap usage errors reproduced byte-for-byte by hand (the flags are optional at the parser since M7). `sheetdiscovery.rs` ports the branch's `process_args` discovery: recursive explode with the `^Undetermined` include and `purged`/`Undetermined` excludes, `_L[0-9]{3}_` lane pairing, output-dir creation, `.purged` naming, assigned dedupe (collapse + a byte-verified port of the branch's leftmost `re.sub(r"(_001)?.f(ast)?q(.gz)?", ...)`), lane grouping, per-sheet-lane selection (skip WARNING, the `Lane '...' has 2 undetermined reads and N assigned reads files` INFO, and the details DEBUG block with the branch's trailing-apostrophe typo). `AppConfig` gains `mode` and `keep_original`. 2 unit tests; the legacy usage errors verified rc 2 byte-identical against the pre-M7 binary.

- [x] **M7.3 — sheet-mode per-lane pipeline (done 2026-11-10).** `pipeline/sheet.rs` (dispatched from `pipeline/mod.rs` on `mode == Sheet`): the per-lane sequential loop in sample-sheet order reproducing the branch's `main()` contract — `Processing lane '...'`, the three `PosixPath`-repr DEBUG blocks (`    Undetermined files: [(...)]`, `    Output files: (...)`, `    Assigned files: [...]`), then per lane the exact contract: `Using exact method`, three `memory_usage` DEBUG pairs (from `/proc/meminfo` — live values, normalized in the differential), `Building python set from undetermined reads...` + `Adding '...' to set` + `Number of undetermined reads: N (MB)`, `Using multiprocessing pool...` / `Purging already assigned reads...` / `Creating {n} child processes...`, the per-file pool lines, `Found {n} duplicates in the undetermined file.` (or `Found no duplicates ... Exiting...`, log-only, the loop continues), and the **D6 writer** (`write_sheet`, a new function): the `Writing purged fastq file(s) ...` lines, no `Removed`/`Kept`/`saved` lines, and — unless `--keep-original` — the `Replacing undetermined fastq file(s) ...` lines (trailing space included) with an in-place rename that destroys the originals (the branch's `Path.replace` semantics; no cross-device fallback, as the branch has none). Sheet-mode approx uses the canonical per-lane M4/D4 contract (the branch's approx path is broken for *all* inputs — the undetermined "files" are tuples and the eager f-strings crash — documented, not ported). `--purge-assigned` is rejected with a clean error in sheet mode (the flag is not part of the branch's CLI). 1 unit test (D6 keep-vs-replace writer).
  - **D9 — the branch's cumulative undetermined set is NOT ported.** `ss_main.py` declares a module global `undetermined_set = set()` and `build_undetermined_set` *adds* to it (the local assignment was lost) and returns the global: the set accumulates across lanes (lane 2's count is `|lane1 ∪ lane2|`, not `|lane2|`) and is pickled to the pool children, which intersect the **cumulative** set. The port keeps the intended fresh per-lane set; the cumulative set cannot change the purged outputs (a record purged from a lane's file carries a name from that file, hence in that lane's own set by construction), so the only observable difference is the per-lane count figure.
  - **Observed branch behaviours (documented, ported as observed):** the pool variant calls `pool.map(pool_task, assigned_list)` without the CLI `log_level`, so the child's named `multiprocessing` logger stays at INFO — the child `Process {pid} started: ... with fork` DEBUG line and the child PID memory pair are **never emitted** by the branch (latent level-drop bug); only the child `Process {pid} ('{name}') finished yielding {n} already assigned reads` line (via the module logger, inherited at the CLI level) is visible — the port omits the suppressed lines. `dnaio.open(r1, r2)` validates that the R1/R2 record names are **identical** (`Reads are improperly paired. Read name '...' in file 1 does not match '...' in file 2.` → `FileFormatError`); the Rust zip does not enforce it (M3 precedent). Clean errors where the branch traceback'd: missing flowcell sheet (`IndexError`), an undetermined pattern without `_L[0-9]{3}_` (`AttributeError`), and the fc_b case below (`None.name` in the details block).

- [x] **M7.4 — sheet-mode differential validation (done 2026-11-10).** Differential vs the real branch Python (worktree of `samplesheet_support` with a full venv: `dnaio`, `polars==1.30.0`, `rbloom`, `psutil`, `rich`, `rich-argparse`, `multiprocess`, `loky`, `biopython`), `COLUMNS=2000` (Rich: no wrapping at that width), three runs on a 3-lane flowcell dataset (lane L001 with a paired undetermined + assigned files incl. a collapse-dedup pair and a no-lane-token file; L002 with an undetermined pair but no selected assigned files; L003 with nothing — exercises the dedupe DEBUG, the lane WARNING, the per-lane INFO/DEBUG blocks, the skip WARNING, the zero-duplicates lane, and the in-place replace) plus a 4th dataset (L002 with a single-read undetermined file):
  - **Runs A (INFO), B (DEBUG) and C (`--keep-original`):** after normalizing the branch's Rich rendering (18-char column prefix — 8-char timestamp column blank after the first record of a second via `omit_repeated_times`, 1-space gap, 9-char level field; the `main.py:NNN` location column; trailing-whitespace strip; the one blank line per lane emitted by the transient `rich.progress.track` bar in `build_undetermined_set`; live memory values; child pids; the count-line MB figure — `norm2.sh` in the M7.4 workspace) the log diff is **exactly one line in every run: the D9 count** (branch `5` = cumulative `|L001 ∪ L002|` vs port `2` = fresh per-lane set). Everything else is byte-identical, including the skip WARNING, the dedupe line, the pool lines, and the `Writing`/`Replacing` lines. Run C: identical directory listings (originals intact + the two `.purged` files) and zero `Replacing` lines on both sides.
  - **fc_b crash case:** the branch dies in `process_args` with `AttributeError: 'NoneType' object has no attribute 'name'` (the details-block f-string evaluates `None.name` on the single-read lane before any lane runs — rc 1, nothing written); the port completes and purges the lane (D5 spirit; documented deviation).
  - **Outputs:** the purged files are **byte-identical after decompression** on both sides for every run (the `.gz` bytes differ only in gzip metadata — the M2 finding, writer vs writer); `--keep-original` leaves both the originals and the `.purged` files in place on both sides.
  - 81 tests total (39 unit incl. 9 samplesheet, 2 sheetdiscovery, 1 sheet writer + 29 fastq + 13 pipeline); clippy clean; fmt clean.

- [x] **M8 — full-scale real-data differential (2026-11-26, legacy mode, exact).** Partial NovaSeqX flowcell `20261106_LH00217_0474_A252C37LT3` (2 of 8 lanes, L001 + L002), `Demultiplexing_0`: 4 undetermined files (≈ 27.6 GB gz; L001 79,278,878 reads, L002 80,029,454) and the `R__Oomen_26_01` assigned tree (984 files, ≈ 181 GB gz; 492 R1 files scanned, R2 ignored as paired). The `Demultiplexing/Undetermined_*` symlinks point to `/srv/...` and are broken on the local copy — not used. Run: `--undetermined-path Demultiplexing_0 --assigned-path Demultiplexing_0/R__Oomen_26_01 --method exact --threads 16`, master Python 0.2.0 vs `fastq-purge-rs` 0.3.0 release build, sequential (no I/O contention), 16 cores / 61 GB RAM, input on an ecryptfs home, `COLUMNS=2000`. Scripts, logs and `time -v` files in `~/Downloads/fastq-purge-realtest/`.
  - **Result:** both rc 0; L001 `Number of undetermined reads: 79278878` on both sides; L002 Rust `80029454`, Python `159308332` (= L001 + L002 — the D10 cumulative count, confirmed at scale); `Found 0 duplicates` in both lanes on both sides → no output files written by either (expected: undetermined and assigned come from the same demultiplexing, so they are disjoint by construction).
  - **Log diff** (Python Rich prefix/location/blank lines stripped): only the documented deviations — M3 (2) the MB figure, M3 (5) pass-2 lines (Rust 492 `Processing target file` + 1 `Scanning ...`; Python 2 × `Using multiprocessing pool...`/`Purging already assigned reads...` and 984 + 984 per-file pool lines — Python re-scans the assigned set once per lane), D10 — plus the output-dir path. The stale-`lane` `Assigned files:` label quirk (M1) is reproduced (all entries logged as `L002`).
  - **Performance:** Rust **17:25 wall, 2.96 GB max RSS** (465 % CPU); Python **39:15 wall, 19.38 GB max RSS** (355 % CPU) ⇒ **2.25× faster, 6.6× less memory** at 159 M undetermined reads.
  - **Rust approx** (same inputs, `--method approx --threads 16`, default `--max-items 1e9`/`--fpr 1e-4` ⇒ 2742.31 MB bloom per lane): rc 0, same read counts as exact (79278878 / 80029454), `Found 0 duplicates` in both lanes — matches exact (no structural FP surfaced in the `Found` count at ≈ 8 % bloom load); **14:28 wall, 5.66 GB max RSS** (511 % CPU) — the two pre-sized 1e9-item blooms dominate memory. Python approx not run (see below).
  - **Not covered:** pass 3 (no duplicates ⇒ nothing written; purged-output identity still rests on M3–M7 and the 2026-11-10 subsampled run); Python approx (master writes a purged copy of every assigned file, ≈ 90 GB, more than the free disk — skipped); sheet mode (every lane of `SampleSheet_0.csv` has a single Recipe `151-151`, so the ≥ 2-Recipes filter skips all lanes).

**Differential testing (the gate)**

1. Build a reference dataset where name0s are **globally unique** (fixes the synthetic-data caveat from profiling: the current generator's coordinate space is seed-invariant and produces birthday-collision name0s, which inflates exact-mode duplicate counts).
2. Run Python exact (threads=1) and the Rust binary on the same set; compare purged outputs record-for-record (names + full records), byte-identical when compression settings match (else record-identical).
3. Property tests: for approx, empirical FPR ≤ `--fpr` on a 10 M membership set + 10 M probes; duplicate count monotonic in `--max-items` sanity.

**Edge cases (unit/integration)**

- [x] empty undetermined file (M6 `empty_undetermined_file_runs_clean`); empty assigned dir (M1 validation differential)
- [x] zero duplicates → no output files, exact log string (M3 `no_duplicates_writes_nothing`)
- [x] all reads duplicated → empty (but existing) purged files (M6 `all_duplicated_writes_empty_purged_files`)
- [x] undetermined R2 missing (`(file, None)` path) (M3 `single_read_lane`)
- [x] ungzipped `.fastq` input (output naming without `.gz`) (M1 validation + M2 plain byte-exact differential)
- [x] multi-lane (4 lanes × R1/R2) with cross-lane duplicates (M5 `two_lanes_cross_lane_dups` + the M5 4-lane dataset, byte-identical outputs)
- [x] assigned dir containing I1/I2 files (collapse + warning) (M1 validation differential)
- [x] filenames without `_L00[0-9]{2}_` (warning + skip) (M1 validation differential)
- [x] re-run with outputs present (purged files excluded from input discovery) (M1/M3 differential re-runs)

**Performance acceptance (5 M synthetic set, 16 cores)**

- [ ] exact end-to-end ≤ 3 s (vs 9.1 s Python) — **not met on the 5 M single-lane set** (18.0–20.7 s @ 16t, ≈ same-day Python 20.5 s): one file per pass ⇒ the single-stream inflate floor (§8/M5 wall analysis); the target assumes multi-lane workloads
- [ ] approx end-to-end ≤ 5 s (vs 38.4 s Python) — **not met single-lane** (19.7 s @ 16t; same-day Python approx 35.8 s, old direction); same wall
- [x] 100 M-name membership ≤ 2.5 GB RSS — **MET**: exact 1.13 GB, approx bloom 0.27 GB (`examples/mem_100m.rs`); Python `set[str]` ≈ 13.2 GB
- [x] ≥ 2× the Python wall time at ≥ 50 M reads — **MET on real data (M8, exact, 2 lanes, 159 M undetermined reads): 2.25×** (17:25 vs 39:15 @ 16t; RSS 2.96 vs 19.38 GB). Earlier: not met at 5 M single-lane (exact ≈ 1.0×, approx 1.8×); met multi-lane (4-lane 2.87× vs its own sequential run, §8/M5)

**Rollout**

1. Milestones: (M1) CLI + discovery/validation port; (M2) FASTQ streaming + name0 hashing; (M3) exact pipeline single-threaded; (M4) bloom pipeline; (M5) rayon + progress; (M6) differential + perf validation. ~5–7 working days total.
2. Ship `fastq-purge-rs` side-by-side; keep the Python package for one release cycle; compare on a real run before declaring the Python version deprecated.
3. Changelog must call out D1–D4 and the fixed approx path (a release task on the Python package — the Rust tree has no changelog).

## 9. Risks & open questions

- **R1 — semantic divergence (D1/D2):** existing users may depend on approx's old direction/key. Mitigation: differential tests, `--purge-assigned`, explicit changelog. *Owner: decide before M1.*
- **R2 — name0 hash collisions (exact):** quantified in §6 (2.7×10⁻⁴ expected pairs at 100 M). Accept; optional `--strict` later.
- **R3 — gzip stream assumptions:** multi-member gzip (Illumina sometimes concatenates members) — `flate2` handles it; fixture already committed as `tests/fastq.rs::multi_member_gz_is_supported` (two concatenated gzip members, both records read back). CRLF line endings are not expected; fail loudly if a line lacks `\n` (corrupt input) rather than guessing.
- **R4 — output determinism:** keep writes strictly sequential per file so outputs are reproducible; only file-parallelism, never record-reordering.
- **R5 — `--max-items` default (1e9) allocates a 2.4 GB bloom for a 10 M read run.** Keep for parity (the approx DEBUG line now reads `sized for {n} items`; the 1e9 default stays — change deferred to the next major).
- **Open:** should pass 2 count duplicates per lane (for the `Found {n} duplicates` log) or return the set only? Counts are needed for parity logging — keep counts (cheap: one u64 per lane).