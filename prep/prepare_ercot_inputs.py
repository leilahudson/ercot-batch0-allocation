# prepare_ercot_inputs.py — build the model's hourly profiles and data center list from
# Thomas Lee's Data-centers data_bundle (github.mit.edu/t-lee/Data-centers), 2025 data.
#
# inputs  (data/raw/, copied from data_bundle/):
#   demand_all.csv   hourly demand by BA and EIA-930 subregion; ERCOT weather zones are
#                    COAS EAST FWES NCEN NRTH SCEN SOUT WEST, plus ERCO total  [MW, Local Standard Time]
#   wind_ba_cf.csv   hourly wind capacity factor by BA (ERCO column)
#   solar_ba_cf.csv  hourly solar capacity factor by BA (ERCO column)
#   dc_at_bus.csv    data center MW at Thomas's network buses: existing, baxtel, cleanview, planned_mw_upper_bound
#   bus100.csv       Thomas's buses: lon, lat, region, geo_eia_zone (weather zone for ERCOT), kV
# outputs (data/):
#   ercot_hourly_2025.csv   8760 rows: zone loads, ERCO total, wind and solar CF (kept so the raw files aren't needed)
#   profiles.csv            the representative hours the model uses (see below)
#   dc_sites.csv            ERCOT data center sites: lat, lon, zone, existing_mw, planned_mw, kV
#
# Representative days (proposal Sec. 4, Hours):
#   summer peak  = day of the highest ERCOT hourly load in Jun-Sep,   stands for 122 days
#   winter peak  = day of the highest ERCOT hourly load in Dec-Feb,   stands for  90 days
#   spring low   = Mar-May day with the lowest daily mean net load,  stands for 153 days (Mar-May, Oct-Nov)
# Net load for picking the spring day uses ERCOT installed wind and solar (WIND_GW, SOLAR_GW below,
# approximate 2025 values); it only affects which day is picked.
#
# run from the batch0 folder:  python3 prep/prepare_ercot_inputs.py

import os
import pandas as pd

RAW, OUT = "data/raw", "data"
ZONES = ["COAS", "EAST", "FWES", "NCEN", "NRTH", "SCEN", "SOUT", "WEST"]
WIND_GW, SOLAR_GW = 40.0, 30.0          # approximate ERCOT installed capacity, 2025
WEIGHTS = {"summer": 122, "winter": 90, "spring": 153}
MIN_SITE_MW = 75.0                      # ERCOT large-load size

# ---- hourly ERCOT series ----
if os.path.exists(f"{RAW}/demand_all.csv"):
    dem = pd.read_csv(f"{RAW}/demand_all.csv", usecols=["Local Standard Time", "ERCO"] + ZONES)
    wcf = pd.read_csv(f"{RAW}/wind_ba_cf.csv", usecols=["Local Standard Time", "ERCO"]).rename(columns={"ERCO": "wind_cf"})
    scf = pd.read_csv(f"{RAW}/solar_ba_cf.csv", usecols=["Local Standard Time", "ERCO"]).rename(columns={"ERCO": "solar_cf"})
    h = dem.merge(wcf, on="Local Standard Time").merge(scf, on="Local Standard Time")
    h = h.rename(columns={"Local Standard Time": "time"})
    h.to_csv(f"{OUT}/ercot_hourly_2025.csv", index=False)
else:
    h = pd.read_csv(f"{OUT}/ercot_hourly_2025.csv")
h["time"] = pd.to_datetime(h["time"])
h["date"] = h["time"].dt.date
h["month"] = h["time"].dt.month
h["netload"] = h["ERCO"] - WIND_GW * 1e3 * h["wind_cf"] - SOLAR_GW * 1e3 * h["solar_cf"]
assert len(h) == 8760, f"expected 8760 hours, got {len(h)}"

def peak_day(months):
    sub = h[h["month"].isin(months)]
    return sub.loc[sub["ERCO"].idxmax(), "date"]

days = {
    "summer": peak_day([6, 7, 8, 9]),
    "winter": peak_day([12, 1, 2]),
    "spring": h[h["month"].isin([3, 4, 5])].groupby("date")["netload"].mean().idxmin(),
}

# zone load factors: each zone's load divided by that zone's 2025 peak, so the case's Pd (taken as
# the zone peak) is scaled by the real hourly shape of its own zone
zpeak = h[ZONES].max()
rows = []
for name, d in days.items():
    day = h[h["date"] == d].sort_values("time")
    assert len(day) == 24, f"{name} day {d} has {len(day)} hours"
    for i, r in enumerate(day.itertuples(index=False), start=1):
        row = {"day": name, "date": str(d), "hour": i, "days_weight": WEIGHTS[name]}
        for z in ZONES:
            row[f"lf_{z}"] = getattr(r, z) / zpeak[z]
        row["load_factor"] = r.ERCO / h["ERCO"].max()       # system-wide, kept for the toy and old code paths
        row["wind_cf"], row["solar_cf"] = r.wind_cf, r.solar_cf
        rows.append(row)
prof = pd.DataFrame(rows)
assert prof["days_weight"].groupby(prof["day"]).first().sum() == 365
prof.to_csv(f"{OUT}/profiles.csv", index=False, float_format="%.5f")

# zone peaks are needed to scale data center MW to the synthetic case (scenario.jl)
zpeak.rename("peak_mw").rename_axis("zone").reset_index().to_csv(f"{OUT}/zone_peaks_2025.csv", index=False)

# ---- data center sites ----
if os.path.exists(f"{RAW}/dc_at_bus.csv"):
    dc = pd.read_csv(f"{RAW}/dc_at_bus.csv")
    bus = pd.read_csv(f"{RAW}/bus100.csv", usecols=["bus_id", "lon", "lat", "region", "geo_eia_zone", "kV", "NAME"])
    s = dc.merge(bus, on="bus_id", how="left")
    s = s[s["region"] == "ERCO"].copy()
    s["existing_mw"] = s["existing"].fillna(0.0)
    s["planned_mw"] = s["planned_mw_upper_bound"].fillna(0.0)          # largest of the Baxtel and Cleanview estimates
    s = s.rename(columns={"bus_id": "src_bus", "geo_eia_zone": "zone", "NAME": "county"})
    s = s[["src_bus", "lat", "lon", "zone", "county", "kV", "existing_mw", "planned_mw"]]
    s.to_csv(f"{OUT}/dc_sites.csv", index=False)
    s = pd.read_csv(f"{OUT}/dc_sites.csv")

    print("representative days:", {k: str(v) for k, v in days.items()})
    print(f"ERCOT 2025 peak {h['ERCO'].max():.0f} MW; zone peaks:", zpeak.round(0).to_dict())
    big = s[s["planned_mw"] >= MIN_SITE_MW]
    print(f"data centers: {len(s)} ERCOT sites, {s['existing_mw'].sum()/1e3:.1f} GW existing, "
          f"{len(big)} planned sites >= {MIN_SITE_MW:.0f} MW totalling {big['planned_mw'].sum()/1e3:.1f} GW")
