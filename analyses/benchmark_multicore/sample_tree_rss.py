#!/usr/bin/env python3
"""Sample the peak aggregate RSS of a process tree.

`/usr/bin/time -f "%M"` reports the largest single process in a tree, not the sum
over concurrently living processes, so it cannot support a claim about total job
memory at a given core count. This walks /proc instead, summing RSS across the
root process and every descendant alive at that instant.

Reports both quantities so they can be compared directly:
  tree_peak_kb  peak of the SUM over all live processes
  proc_peak_kb  peak of the LARGEST single process (comparable to %M)

usage: sample_tree_rss.py <root_pid> <interval_seconds> <out_tsv>
"""

import os
import sys
import time

PAGE_KB = os.sysconf("SC_PAGE_SIZE") // 1024


def children_map():
    """pid -> [child pids], built from /proc/<pid>/stat."""
    kids = {}
    for entry in os.listdir("/proc"):
        if not entry.isdigit():
            continue
        try:
            with open(f"/proc/{entry}/stat", "rb") as fh:
                data = fh.read()
        except (OSError, ProcessLookupError):
            continue
        # comm is parenthesised and may contain spaces; ppid is the 2nd field after it
        close = data.rfind(b")")
        if close == -1:
            continue
        fields = data[close + 2 :].split()
        if len(fields) < 2:
            continue
        try:
            ppid = int(fields[1])
        except ValueError:
            continue
        kids.setdefault(ppid, []).append(int(entry))
    return kids


def rss_kb(pid):
    """Resident set size of one pid, in kB."""
    try:
        with open(f"/proc/{pid}/statm", "rb") as fh:
            return int(fh.read().split()[1]) * PAGE_KB
    except (OSError, IndexError, ValueError):
        return 0


def tree(root, kids):
    """root and every descendant currently alive."""
    out, stack = [], [root]
    while stack:
        pid = stack.pop()
        out.append(pid)
        stack.extend(kids.get(pid, ()))
    return out


def alive(pid):
    return os.path.isdir(f"/proc/{pid}")


def main():
    root = int(sys.argv[1])
    interval = float(sys.argv[2])
    out_path = sys.argv[3]

    tree_peak = proc_peak = n_at_peak = samples = 0
    while alive(root):
        kids = children_map()
        pids = tree(root, kids)
        sizes = [rss_kb(p) for p in pids]
        total = sum(sizes)
        samples += 1
        if total > tree_peak:
            tree_peak, n_at_peak = total, len([s for s in sizes if s])
        biggest = max(sizes, default=0)
        if biggest > proc_peak:
            proc_peak = biggest
        time.sleep(interval)

    with open(out_path, "w") as fh:
        fh.write(f"{tree_peak}\t{proc_peak}\t{n_at_peak}\t{samples}\n")


if __name__ == "__main__":
    main()
