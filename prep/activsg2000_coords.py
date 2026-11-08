# activsg2000_coords.py — bus coordinates for ACTIVSg2000 from TAMU's PowerWorld auxiliary file.
# The MATPOWER case has no latitude/longitude; the .aux file does, on its substations.
#
# Download the ACTIVSg2000 dataset from https://electricgrids.engr.tamu.edu (ACTIVSg2000 page,
# "Download Dataset"), unzip, and copy the .aux file to data/raw/. Then, from the batch0 folder:
#     python3 prep/activsg2000_coords.py data/raw/ACTIVSg2000.aux
# writes data/activsg2000_bus_coords.csv (bus_id, lat, lon), which scenario.jl picks up automatically
# and uses to place each data center at its nearest high-voltage bus.
#
# AUX format: blocks like
#   DATA (Substation, [Number, Name, ..., Latitude, Longitude, ...])
#   { 1 "NAME" ... 31.85 -102.37 ... }
#   DATA (Bus, [BusNum, BusName, ..., SubNumber, ...])
# Field order differs between files, so columns are found by name.

import re
import shlex
import sys
import pandas as pd

path = sys.argv[1] if len(sys.argv) > 1 else "data/raw/ACTIVSg2000.aux"
txt = open(path, encoding="latin-1").read()

def block(obj):
    # header may span lines; rows are between the braces that follow it
    m = re.search(r"DATA\s*\(\s*" + obj + r"\s*,\s*\[(.*?)\]\s*\)\s*\{(.*?)\}", txt, re.S | re.I)
    if m is None:
        raise SystemExit(f"no DATA ({obj}, ...) block in {path}")
    fields = [f.strip().split(":")[0] for f in m.group(1).replace("\n", " ").split(",") if f.strip()]
    rows = []
    for ln in m.group(2).splitlines():
        ln = ln.split("//")[0].strip()
        if ln:
            rows.append(shlex.split(ln, posix=True))
    rows = [r[:len(fields)] for r in rows if len(r) >= len(fields)]
    return pd.DataFrame(rows, columns=fields)

def col(df, *names):
    for n in names:
        for c in df.columns:
            if c.lower() == n.lower():
                return c
    raise SystemExit(f"none of {names} in {list(df.columns)[:20]}...")

sub = block("Substation")
bus = block("Bus")
s = pd.DataFrame({"sub": sub[col(sub, "Number", "SubNumber")].astype(int),
                  "lat": sub[col(sub, "Latitude")].astype(float),
                  "lon": sub[col(sub, "Longitude")].astype(float)})
b = pd.DataFrame({"bus_id": bus[col(bus, "BusNum", "Number")].astype(int),
                  "sub": bus[col(bus, "SubNumber", "Sub Num")].astype(int)})
out = b.merge(s, on="sub", how="left")[["bus_id", "lat", "lon"]]
print(f"{len(out)} buses, {out['lat'].notna().sum()} with coordinates "
      f"(lat {out['lat'].min():.2f}..{out['lat'].max():.2f}, lon {out['lon'].min():.2f}..{out['lon'].max():.2f})")
out.to_csv("data/activsg2000_bus_coords.csv", index=False)
