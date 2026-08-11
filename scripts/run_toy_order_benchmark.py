#!/usr/bin/env python3
"""Run the deterministic order-dependence benchmark."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from hgsoc_corneto.metabolic.toy import toy_order_benchmark


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--semantics", choices=("published", "bounds_safe"), default="bounds_safe")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    result = toy_order_benchmark(args.semantics)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))
    if not result["order_dependence_observed"]:
        raise SystemExit("Toy benchmark unexpectedly failed to show order dependence")


if __name__ == "__main__":
    main()
