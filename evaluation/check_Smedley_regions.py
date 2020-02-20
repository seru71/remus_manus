#!/usr/bin/env python3

import os
import sys, time
import pybedtools
from remus_client import get_matching_tissues_strict, get_all_tissues, query


class RemusResultMatch:
    def __init__(self, matching_interval, closest, closest_dist, count, length):
        self.matching_interval = matching_interval
        self.closest_interval = closest
        self.closest_interval_dist = closest_dist
        self.count = count
        self.length = length


def format_interval(s):
    s = (s.split('\n')[0]).split('\t')
    return f"{s[0]}:{s[1]}-{s[2]}"


def in_bed(bed, position):
    if len(bed) == 0:
        return RemusResultMatch(None, None, -1, 0, 0)

    bed = bed.sort().merge()
    p_bed = pybedtools.BedTool(f"{bed[0].chrom}\t{position - 1}\t{position}", from_string=True)
    matching_interval = bed.intersect(p_bed, wa=True)
    c = str(p_bed.closest(bed, d=True)).split('\t')
    closest = format_interval('\t'.join(c[3:6]))
    closest_dist = int(c[6])
    count = len(bed)
    length = sum([len(i) for i in bed])

    return RemusResultMatch(None if not matching_interval
                            else format_interval(str(matching_interval)), closest, closest_dist, count, length)


def query_args_builder(genes, tissues, **kwargs):
    request_args = {'genes': genes, 'tissues': tissues, 'genome': "hg19",
                    "genes-select-include-gene-transcripts": None}
    # set default paramters for the query
    request_args.update({"promoters-fantom5-used": "yes",
                         "promoters-fantom5-combine-mode": "any",
                         "promoters-fantom5-kbs-upstream": "100",
                         "promoters-fantom5-kbs-downstream": "3000"})
    request_args.update({"promoters-screen-used": "yes",
                         "promoters-screen-combine-mode": "any",
                         "promoters-screen-kbs-upstream": "100",
                         "promoters-screen-kbs-downstream": "3000"})
    request_args.update({"enhancers-encode-used": "yes",
                         "enhancers-encode-combine-mode": "any",
                         "enhancers-encode-kbs-upstream": 500,
                         "enhancers-encode-kbs-downstream": 500})
    request_args.update({"enhancers-screen-used": "yes",
                         "enhancers-screen-combine-mode": "any",
                         "enhancers-screen-kbs-upstream": 500,
                         "enhancers-screen-kbs-downstream": 500})
    request_args.update({"enhancers-fantom5-used": "yes",
                         "enhancers-fantom5-combine-mode": "any",
                         "enhancers-fantom5-kbs-upstream": 500,
                         "enhancers-fantom5-kbs-downstream": 500})
    request_args.update({"accessible-chromatin-encode-used": "yes",
                         "accessible-chromatin-encode-combine-mode": "any",
                         "accessible-chromatin-encode-kbs-upstream": 500,
                         "accessible-chromatin-encode-kbs-downstream": 500})
    request_args.update({"accessible-chromatin-screen-used": "yes",
                         "accessible-chromatin-screen-combine-mode": "any",
                         "accessible-chromatin-screen-kbs-upstream": 500,
                         "accessible-chromatin-screen-kbs-downstream": 500})
    request_args.update(kwargs)
    return request_args


def query_regulator(gene, track_list, max_attempts=3, upstream_dist=None, downstream_dist=None):
    endpoint = "http://remus.btm.umed.pl"
    tissues = None
    attempt = 1

    extra_args={}
    if upstream_dist:
        extra_args.update({"enhancers-encode-kbs-upstream": str(upstream_dist),
                           "enhancers-fantom5-kbs-upstream": str(upstream_dist),
                           "accessible-chromatin-encode-kbs-upstream": str(upstream_dist)})
    if downstream_dist:
        extra_args.update({"enhancers-encode-kbs-downstream": str(downstream_dist),
                           "enhancers-fantom5-kbs-downstream": str(downstream_dist),
                           "accessible-chromatin-encode-kbs-downstream": str(downstream_dist)})

    while attempt < max_attempts:
        try:
            if not track_list or track_list==[]:
                sys.stderr.write(f"Querying {gene} in ALL tissues. Attempt: {attempt}\n")
                if not tissues:
                    tissues = get_all_tissues(endpoint, "hg19")
            else:
                sys.stderr.write(f"Querying {gene} in tissues: {track_list}. Attempt: {attempt}\n")
                if not tissues:
                    tissues = get_matching_tissues_strict(endpoint, "hg19", track_list)

            return pybedtools.BedTool(query(endpoint, query_args_builder(gene, tissues, **extra_args)).decode(),
                                      from_string=True), tissues

        except Exception as e:
            sys.stderr.write(f"Exception: {e}\n")
            attempt += 1
            time.sleep(1)


def print_result(r, gene, tissues, position, lookup_dist):
    if r.matching_interval:
        print(f"FOUND!. Regulator for {gene} in position {position} is in the result " +
              f"({r.count} features; {r.length / 1000}kb; lookup: {lookup_dist}kb) for tissues {tissues}: {r.matching_interval}.")
    else:
        print(f"Regulator for {gene} in position {position} is not in the result " +
              f"({r.count} features; {r.length / 1000}kb; lookup: {lookup_dist}kb) for tissues {tissues}.\n" +
              f"Closest feature is {r.closest_interval} located {r.closest_interval_dist}bp away.")


def print_result_header():
    print('\t'.join(["gene", "tissues", "num_features", "features_size", "regulator_overlapping", "regulator_closest", "regulator_closest_distance",
                     "tissues_WIDE", "num_features_WIDE", "features_size_WIDE", "regulator_overlapping_WIDE", "regulator_closest_WIDE", "regulator_closest_distance_WIDE",
                     "num_features_ALL", "features_size_ALL", "regulator_overlapping_ALL", "regulator_closest_ALL", "regulator_closest_distance_ALL"]))


def print_result_row(chrom, pos, gene, tissues, results):
    ts, ts_wide = tissues
    res, res_wide, res_all = results
    res_cols = [ts, res.count, res.length, res.matching_interval, res.closest_interval, res.closest_interval_dist]
    res_wide_cols = [ts_wide, res_wide.count, res_wide.length, res_wide.matching_interval, res_wide.closest_interval, res_wide.closest_interval_dist]
    res_all_cols = [res_all.count, res_all.length, res_all.matching_interval, res_all.closest_interval, res_all.closest_interval_dist]
    print('\t'.join([str(e) for e in [chrom, pos, gene] + res_cols + res_wide_cols + res_all_cols]))



if __name__ == '__main__':

    with open(sys.argv[1]) as regs:
        prev_gene_tissue = ("", "", "", 500)
        bed, bed_wide, bed_all = None, None, None
        tissues, tissues_wide = None, None
    
        print_result_header()
  
        for reg in regs:
            chrom, position, gene, track1, tracks, dist = reg.split('\t')
            lookup_dist = max(500, round((int(dist)+100)/1000))

            if prev_gene_tissue != (gene, track1, tracks, lookup_dist):
                prev_gene_tissue = (gene, track1, tracks, lookup_dist)
                time.sleep(1)
                us, ds = (None,None) if lookup_dist==500 else (lookup_dist, lookup_dist)
                bed, tissues = query_regulator(gene, [track1], upstream_dist=us, downstream_dist=ds)
                bed.saveas(os.path.join("result_beds", gene + "_" + track1.replace(" ","_") + ".bed"))
                bed_wide, tissues_wide = query_regulator(gene, tracks.split(';'), upstream_dist=us, downstream_dist=ds)
                bed_wide.saveas(os.path.join("result_beds", gene + "_wide.bed"))
                bed_all, _ = query_regulator(gene, [], upstream_dist=us, downstream_dist=ds)
                bed_all.saveas(os.path.join("result_beds", gene + "_alltissues.bed"))

            r1 = in_bed(bed, int(position))
            r2 = in_bed(bed_wide, int(position))
            r3 = in_bed(bed_all, int(position))
            print_result_row(chrom, position, gene, (';'.join(tissues), ';'.join(tissues_wide)), (r1, r2, r3))
