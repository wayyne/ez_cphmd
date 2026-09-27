#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <arex.history.tsv> '<pH values>'" >&2
    echo "Example: $0 arex.history.tsv '4.5 5.0 5.5 6.0 6.5 7.0 7.5 8.0'" >&2
    exit 1
fi

HIST=$1
PHS=$2

if [[ ! -f "$HIST" ]]; then
    echo "Error: history file not found: $HIST" >&2
    exit 1
fi

python3 - "$HIST" "$PHS" <<'PY'
import sys

hist = sys.argv[1]
phs = [float(x) for x in sys.argv[2].split()]
n = len(phs)

if n < 2:
    raise SystemExit("Error: at least two pH values are required.")

rows = []

with open(hist) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        toks = line.split()

        # Skip header:
        # epoch replica_1 replica_2 ...
        if toks[0].lower() == "epoch":
            continue

        try:
            epoch = int(toks[0])
            vals = list(map(int, toks[1:]))
        except ValueError:
            continue

        if len(vals) != n:
            continue

        # Native AREX history is slot_of_replica, 1-based:
        # columns replica_1 ... replica_N tell which slot each replica occupies.
        slot_of_replica = vals

        owner_by_slot = [None] * n
        for rep0, slot in enumerate(slot_of_replica):
            slot0 = slot - 1
            if 0 <= slot0 < n:
                owner_by_slot[slot0] = rep0

        if any(x is None for x in owner_by_slot):
            continue

        rows.append((epoch, owner_by_slot))

if len(rows) < 2:
    raise SystemExit("Error: fewer than two valid history rows were found.")

attempts = [0] * (n - 1)
accepts = [0] * (n - 1)

for k in range(1, len(rows)):
    epoch, curr = rows[k]
    prev_epoch, prev = rows[k - 1]

    # Native first parity appears to be 1-based odd pairs:
    # slots (1,2), (3,4), ... when epoch is odd
    # slots (2,3), (4,5), ... when epoch is even
    #
    # Convert to 0-based pair starts:
    # odd epoch  -> 0,2,4,6
    # even epoch -> 1,3,5
    pair_start = 0 if (epoch % 2 == 1) else 1

    for i in range(pair_start, n - 1, 2):
        attempts[i] += 1
        if curr[i] == prev[i + 1] and curr[i + 1] == prev[i]:
            accepts[i] += 1

print("pair          attempts  accepts  accept_rate")
for i in range(n - 1):
    rate = accepts[i] / attempts[i] if attempts[i] else 0.0
    print(f"{phs[i]:4.1f}-{phs[i+1]:4.1f}  {attempts[i]:8d}  {accepts[i]:7d}  {rate:10.4f}")
PY
