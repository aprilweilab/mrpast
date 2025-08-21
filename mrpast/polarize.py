from dataclasses import dataclass, field
from typing import List, Optional, TextIO
import pyfaidx


@dataclass
class PolarizeStats:
    unpolarized: List[int] = field(default_factory=list)
    num_allele_notfound: int = 0
    num_no_align: int = 0  # Alignment doesn't match (N, ., - in the FASTA)
    num_after: int = 0
    num_sites: int = 0  # Total sites
    flipped: int = 0
    unphased: int = 0

    def print(self, out, prefix=""):
        print(f"{prefix}Alleles didn't match: {self.num_allele_notfound}", file=out)
        print(f"{prefix}Alignment mismatch (N, ., -): {self.num_no_align}", file=out)
        print(f"{prefix}After end of alignemnt: {self.num_after}", file=out)
        print(f"{prefix}Unphased (skipped): {self.unphased}", file=out)
        print(f"{prefix}Flipped: {self.flipped}", file=out)
        print(f"{prefix}Total sites seen: {self.num_sites}", file=out)


def int_or_miss(string: str):
    if string == ".":
        return string
    return int(string)


def _make_new_genotype(old_genotype: List[str], swap_index: int) -> Optional[List[str]]:
    """
    Swap the ALT swap_index (0-based) with REF. The input is a list like ["0|1", "1|1", ...]
    as is the output.
    """
    new_genotype = []
    for indiv in old_genotype:
        if "/" in indiv:
            return None
        new_indices = []
        for index in map(int_or_miss, indiv.split("|")):
            if index == (swap_index + 1):
                new_indices.append(0)
            elif index == 0:
                new_indices.append(swap_index + 1)
            else:
                new_indices.append(index)
        new_genotype.append("|".join(map(str, new_indices)))
    return new_genotype


def polarize_vcf(
    in_file_obj: TextIO, out_file_obj: TextIO, anc_fasta: str, drop_no_anc: bool = True
) -> PolarizeStats:
    """
    Polarize a phased VCF file.
    """
    stats = PolarizeStats()

    fasta_reader = pyfaidx.Fasta(anc_fasta)
    assert (
        len(fasta_reader.values()) == 1
    ), "Your FASTA file has more than one contig; we only support a single ancestral contig."
    ancestral_str = "X" + str(list(fasta_reader.values())[0])

    for line in in_file_obj:
        line = line.strip()
        if line.startswith("##"):
            print(line, file=out_file_obj)
            continue
        elif line.startswith("#"):
            print(line, file=out_file_obj)
        else:
            stats.num_sites += 1

            data = line.split("\t")
            chrom, position, var_id, ref, alt, qual, filt, info, fmt = data[0:9]
            all_alts = alt.split(",")
            site_pos = int(position)
            new_genotype = None
            skip = False
            if site_pos < len(ancestral_str):
                ref = ref.upper()
                aa = ancestral_str[site_pos].upper()
                if aa in "N.-":
                    stats.num_no_align += 1
                    skip = True
                elif aa != ref:
                    swap_index = -1
                    all_alts = list(map(str.upper, all_alts))
                    for i, alt_allele in enumerate(all_alts):
                        if alt_allele == aa:
                            swap_index = i
                            break
                    if swap_index >= 0:
                        new_genotype = _make_new_genotype(data[9:], swap_index)
                        if new_genotype is not None:
                            stats.flipped += 1
                            all_alts[swap_index] = ref
                            ref = aa
                            alt = ",".join(all_alts)
                        else:
                            skip = True
                            stats.unphased += 1
                    else:
                        stats.num_allele_notfound += 1
                        skip = True
            else:
                skip = True
                stats.num_after += 1

            if new_genotype is None:
                if skip:
                    stats.unpolarized.append(site_pos)
                    if not drop_no_anc:
                        print(line, file=out_file_obj)
                else:
                    print(line, file=out_file_obj)
            else:
                print(
                    "\t".join(
                        [chrom, str(site_pos), var_id, ref, alt, qual, filt, info, fmt]
                        + new_genotype
                    ),
                    file=out_file_obj,
                )
    return stats
