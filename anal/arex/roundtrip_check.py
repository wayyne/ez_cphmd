#!/usr/bin/env python3
import sys

if len(sys.argv) != 3:
    print(
        "Usage: python3 roundtrip_check.py arex.history.tsv "
        "'4.5 5.0 5.5 6.0 6.5 7.0 7.5 8.0'"
    )
    sys.exit(1)

hist = sys.argv[1]
phs = [float(x) for x in sys.argv[2].split()]
n = len(phs)

if n < 2:
    print("Error: at least two pH values are required.", file=sys.stderr)
    sys.exit(1)

rows = []

with open(hist) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        toks = line.split()
        if toks[0].lower() == "epoch":
            continue

        try:
            epoch = int(toks[0])
            slot_of_replica = list(map(int, toks[1:]))
        except ValueError:
            continue

        if len(slot_of_replica) != n:
            continue

        # Convert 1-based slot_of_replica to 0-based.
        slot_of_replica = [s - 1 for s in slot_of_replica]

        if any(s < 0 or s >= n for s in slot_of_replica):
            continue

        rows.append((epoch, slot_of_replica))

if not rows:
    print("No valid rows found.")
    sys.exit(1)

# Occupancy: replica -> slot counts
visits = [[0 for _ in range(n)] for _ in range(n)]

# Complete round trips per replica. The first endpoint reached becomes the anchor;
# a round trip is counted only after the replica reaches the opposite endpoint and
# then returns to its anchor endpoint.
roundtrips = [0 for _ in range(n)]
anchor_endpoint = [None for _ in range(n)]  # 0 or n-1
reached_opposite = [False for _ in range(n)]

# Transition distance per replica
moves = [0 for _ in range(n)]
total_abs_slot_jump = [0 for _ in range(n)]

prev_slots = None

for epoch, slot_of_replica in rows:
    for rep in range(n):
        slot = slot_of_replica[rep]
        visits[rep][slot] += 1

        if slot == 0 or slot == n - 1:
            if anchor_endpoint[rep] is None:
                anchor_endpoint[rep] = slot
            elif slot == anchor_endpoint[rep]:
                if reached_opposite[rep]:
                    roundtrips[rep] += 1
                    reached_opposite[rep] = False
            else:
                reached_opposite[rep] = True

    if prev_slots is not None:
        for rep in range(n):
            jump = abs(slot_of_replica[rep] - prev_slots[rep])
            if jump:
                moves[rep] += 1
                total_abs_slot_jump[rep] += jump

    prev_slots = slot_of_replica

print(f"Rows analyzed: {len(rows)}")
print()

print("Replica slot occupancy fractions")
header = "replica " + " ".join([f"pH{p:>4.1f}" for p in phs])
print(header)
for rep in range(n):
    total = sum(visits[rep])
    fracs = [visits[rep][slot] / total if total else 0.0 for slot in range(n)]
    print(f"{rep+1:7d} " + " ".join(f"{x:6.3f}" for x in fracs))

print()
print("Replica round-trip / movement summary")
print("replica  roundtrips  visited_slots  min_pH_seen  max_pH_seen  moves  avg_abs_jump_when_moved")
for rep in range(n):
    visited_slots = sum(1 for c in visits[rep] if c > 0)
    avg_jump = total_abs_slot_jump[rep] / moves[rep] if moves[rep] else 0.0
    min_seen = min(phs[i] for i, c in enumerate(visits[rep]) if c > 0)
    max_seen = max(phs[i] for i, c in enumerate(visits[rep]) if c > 0)
    print(
        f"{rep+1:7d}  {roundtrips[rep]:10d}  {visited_slots:13d}  "
        f"{min_seen:11.1f}  {max_seen:11.1f}  {moves[rep]:5d}  {avg_jump:23.3f}"
    )

print()
print("Slot occupancy totals")
print("slot/pH  total_visits")
for slot, ph in enumerate(phs):
    total = sum(visits[rep][slot] for rep in range(n))
    print(f"{slot+1:4d}/{ph:4.1f}  {total:12d}")
