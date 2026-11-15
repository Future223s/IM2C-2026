from __future__ import annotations

from pathlib import Path

import pandas as pd
import streamlit as st


ROOT = Path(__file__).resolve().parent
TABLES = ROOT / "outputs" / "tables"
FIGURES = ROOT / "outputs" / "figures"

st.set_page_config(page_title="Autumn phenology coverage", layout="wide")
st.title("Multi-region autumn phenology dataset")
st.caption("Release 0.1 — PhenoCam observations + Daymet environment; verified pending modalities remain explicit.")

sites = pd.read_parquet(TABLES / "sites.parquet")
coverage = pd.read_parquet(TABLES / "coverage_audit.parquet")
ranking = pd.read_parquet(TABLES / "site_selection_ranking.parquet")
inventory = pd.read_parquet(TABLES / "dataset_inventory.parquet")
events = pd.read_parquet(TABLES / "phenology_events.parquet")

c1, c2, c3, c4 = st.columns(4)
c1.metric("Ecological sites", len(sites))
c2.metric("Trajectory-ready site-years", int(coverage["phenocam_autumn_sufficient"].sum()))
c3.metric("Camera event estimates", len(events))
c4.metric("Pending inventory entries", int(inventory["availability"].str.contains("pending|verified accessible|identifier verified", case=False, regex=True).sum()))

st.subheader("Coverage")
st.image(str(FIGURES / "phenocam_autumn_coverage.png"), use_container_width=True)

site_id = st.selectbox("Inspect a site", ranking["site_id"], format_func=lambda value: ranking.set_index("site_id").loc[value, "site_name"])
left, right = st.columns(2)
with left:
    st.write("Site-year audit")
    st.dataframe(coverage[coverage["site_id"].eq(site_id)].sort_values("year"), use_container_width=True, hide_index=True)
with right:
    st.write("Dataset inventory")
    st.dataframe(inventory[inventory["site_id"].eq(site_id)], use_container_width=True, hide_index=True)

st.subheader("RCC90 operational events")
st.warning("These are observed camera RCC90 maxima, not validated landscape peak-color dates.")
st.image(str(FIGURES / "rcc90_maxima_variation.png"), use_container_width=True)

st.subheader("Current-release ranking")
st.dataframe(ranking, use_container_width=True, hide_index=True)

