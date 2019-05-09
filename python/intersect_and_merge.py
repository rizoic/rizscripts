#!/usr/bin/env python3
# --------------- intersect and merge --------------------
# intersect and merge a set of bed file to give regions present
# in atleast given number of files

import argparse
import os
from functools import reduce

import pybedtools


def main(args):
    # Merge the individual files to ensure later in count we dont count two intervals from the same file
    indiv_merged = [pybedtools.BedTool(bedfile).sort().merge(d=args.indivmergedist) for bedfile in args.bedfiles]

    # Concatenate all these merged files into one
    all_cat = reduce(lambda x, y: x.cat(y, postmerge=False, force_truncate=True), indiv_merged)

    # Now merge and get count for the number of files an interval is present in
    all_merged = all_cat.sort().merge(c=1, o="count").filter(lambda x: int(x.name) >= args.minfiles).sort()

    # Print just the chrom start and stop coordinates
    for feature in all_merged:
        print("{f.chrom}\t{f.start}\t{f.stop}".format(f=feature))


if __name__ == '__main__':

    epilog = "EXAMPLE: python " + os.path.basename(__file__) + \
             " --bed /path/to/a.bed --bed /path/to/b.bed --bed /path/to/c.bed" \
             " --minfiles 3 --indivmergedist 0"

    parser = argparse.ArgumentParser(description="Intersect and Merge bed files",
                                     epilog=epilog)

    required_args_group = parser.add_argument_group('required arguments')

    required_args_group.add_argument('-b', '--bed',
                                     dest='bedfiles',
                                     required=True,
                                     help="Bed file to intersect and merge",
                                     action='append')

    parser.add_argument('-m', '--minfiles',
                        dest='minfiles',
                        type=int,
                        help="Interval present in atleast x number of files. Default: all")

    parser.add_argument('-d', '--indivmergedist',
                        dest='indivmergedist',
                        default=1,
                        type=int,
                        help="Distance upon which to merge individual intervals. Default: 1")

    args = parser.parse_args()

    if not args.minfiles:
        args.minfiles = len(args.bedfiles)

    main(args)
