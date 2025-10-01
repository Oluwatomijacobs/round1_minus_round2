#!/usr/bin/env python3 
# the line above lets you run the file directly chmod +x …; ./script.py and calls on the python 3 on your system.

from pathlib import Path # clean path handling
from collections import Counter, defaultdict # counter helps with fast frequency count per cluster id
import re, sys # regular expressions helps with extracting numbers /sys - exit with error, write ith stderr

# ----- folders -----
# .expanduser() turns ~ into your real home directory (Python doesn’t do this i.e expland ~ automatically).
R1 = Path('~/data/fasta/prelim/round1/subtrees1/fastafromFiles').expanduser()
R2 = Path('~/data/fasta/prelim/round2/subtrees2/fastafromFiles2').expanduser()

# Regex: precompiled pattern that finds one or more digits (0–9) in a filename. This will be used repeatedly.
DIGITS_RE = re.compile(r'\d+')

# So for my case, files were renamed in round2 but cluster IDs at the end remained the same. 
# So last numeric token links Round-1 and Round-2 files even if prefixes changed.
def id_lastnum(name: str) -> str:
    """Return the LAST run of digits in the filename stem (cluster ID)."""
    stem = Path(name).stem # this strips the .fa extension
    hits = DIGITS_RE.findall(stem) # returns all digit runs in the stem; we take the last one
    return hits[-1] if hits else stem.lower() # this line ensures If there are no digits (edge case), we fall back to the lowercased stem so the function still returns a key

# Safety: fail fast if the directory doesn’t exist (avoids silently returning zero).
def list_fa_names(d: Path) -> list[str]:
    if not d.is_dir():
        sys.exit(f"ERROR: not a directory: {d}")
    return [p.name for p in d.glob('*.fa') if p.is_file()] #this collects only .fa files (non-recursive), and we return filenames (p.name) because we compare names/IDs, not full paths.

def multiset_missing_by_id(r1_names: list[str], r2_names: list[str]) -> list[str]:
    """Return R1 filenames missing in R2 (matching by last-number ID; handles duplicates)."""
    c1 = Counter(id_lastnum(n) for n in r1_names)
    c2 = Counter(id_lastnum(n) for n in r2_names)

    id2r1 = defaultdict(list) # this keeps for each cluster ID, the actual Round-1 filenames with that ID. This lets us output real filenames (not just IDs) for the missing ones.
    for n in r1_names:
        id2r1[id_lastnum(n)].append(n)

# multi-set logic that helps handles duplicates...well if any
    missing = []
    for cid, cnt1 in c1.items():
        miss = cnt1 - c2.get(cid, 0) # computes how many R1 files are “unmatched”
        if miss > 0: # If positive, this takes that many filenames from the R1 pool and marks them as missing in Round-2.
            missing.extend(id2r1[cid][:miss]) 

    # Sort by numeric ID if possible, then by name (stable output)
    # primarily by numeric cluster ID (so 2 < 10), secondarily by filename, to break ties deterministically.
    def keyfn(n: str):
        stem = Path(n).stem
        hits = DIGITS_RE.findall(stem)
        if hits:
            try:
                return (int(hits[-1]), n)
            except ValueError:
                return (hits[-1], n)
        return (stem.lower(), n)

    return sorted(set(missing), key=keyfn) # set(missing) ensures there are no duplicate filenames (filenames are expected unique within a folder).
    # If you ever truly need to preserve duplicate filenames (rare), remove the set()

def main():
    r1 = list_fa_names(R1)
    r2 = list_fa_names(R2)
    expected = len(r1) - len(r2)

    missing = multiset_missing_by_id(r1, r2) # produces the actual list of Round 1 filenames missing in Round 2 using the last-number matching.

    out_csv = R1 / 'only_in_R1.by_lastnum.csv'
    with out_csv.open('w') as f:
        for name in missing:
            f.write(name + '\n')   # <-- single column: filename only

    print(f"# R1(.fa)={len(r1)}  R2(.fa)={len(r2)}  expected_missing={expected}", file=sys.stderr)
    print(f"# wrote {len(missing)} lines -> {out_csv}", file=sys.stderr)

if __name__ == '__main__':
    main()
