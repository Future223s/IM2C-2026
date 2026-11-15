from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys

import requests


PRODUCT = "DP1.10055.001"
QUERY_URL = "https://data.neonscience.org/api/v0/data/query"


def collect_urls(value: object) -> list[str]:
    urls: list[str] = []
    if isinstance(value, dict):
        for key, child in value.items():
            if key == "url" and isinstance(child, str):
                urls.append(child)
            else:
                urls.extend(collect_urls(child))
    elif isinstance(value, list):
        for child in value:
            urls.extend(collect_urls(child))
    return urls


def main() -> None:
    parser = argparse.ArgumentParser(description="Download official NEON plant phenology files without flattening individual records.")
    parser.add_argument("--sites", default="HARV,BART,GRSM,ORNL")
    parser.add_argument("--start", default="2013-01")
    parser.add_argument("--end", default="2025-12")
    parser.add_argument("--release", default="RELEASE-2026")
    parser.add_argument("--output", type=Path, default=Path("cache/neon/DP1.10055.001"))
    args = parser.parse_args()

    token = os.environ.get("NEON_API_TOKEN")
    if not token:
        sys.exit("Set NEON_API_TOKEN in the environment. Never commit the token.")
    headers = {"X-API-Token": token, "User-Agent": "IM2C-multi-region-phenology/0.1"}
    args.output.mkdir(parents=True, exist_ok=True)
    for site in [value.strip().upper() for value in args.sites.split(",") if value.strip()]:
        params = {
            "productCode": PRODUCT,
            "siteCode": site,
            "startDateMonth": args.start,
            "endDateMonth": args.end,
            "release": args.release,
            "package": "basic",
            "includeProvisional": "false",
        }
        response = requests.get(QUERY_URL, params=params, headers=headers, timeout=(20, 180))
        response.raise_for_status()
        payload = response.json()
        metadata_path = args.output / f"{site}_query.json"
        metadata_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        for url in sorted(set(collect_urls(payload))):
            filename = url.rsplit("/", 1)[-1].split("?", 1)[0]
            if not filename:
                continue
            path = args.output / site / filename
            if path.exists():
                continue
            path.parent.mkdir(parents=True, exist_ok=True)
            file_response = requests.get(url, headers=headers, timeout=(20, 300))
            file_response.raise_for_status()
            path.write_bytes(file_response.content)
            print(path)


if __name__ == "__main__":
    main()

