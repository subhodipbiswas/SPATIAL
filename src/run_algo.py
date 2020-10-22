#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
from spatial import run_spatial


def make_runs(options):
    """
    Simulate runs for the search methods.
    """
    run_spatial(options)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-q", "--quiet", action="store_false", dest="verbose",
                        help="don't print status messages to stdout")
    parser.add_argument("-r", "--runs", type=int, default=25)  # number of runs to be simulated
    parser.add_argument("-e", "--seed", type=int, default=17)  # integer seed for random number generator
    parser.add_argument("-a", "--algo", type=str, default="SPATIAL")  # algorithms: SPATIAL, HC, OBA, SA, TS
    parser.add_argument("-s", "--school", type=str, default="ES")  # schools: ES, MS, HS
    parser.add_argument("-i", "--initialization", default=1, type=int)  # 1: seeded, 2: infeasiible 3: existing
    parser.add_argument("-y", "--year", type=int, default=2019)  # school year

    options = parser.parse_args()
    make_runs(options)

    return 0


if __name__ == "__main__":
    print(sys.platform)
    sys.exit(main())
