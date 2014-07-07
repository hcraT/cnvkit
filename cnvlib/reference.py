"""Supporting functions for the 'reference' command."""
from __future__ import absolute_import, division

import numpy
from Bio._py3k import map, range, zip

from . import core, metrics, ngfrills, params
from .ngfrills import echo
from .cnarray import CopyNumArray as CNA, row2label


def bed2probes(bed_fname):
    """Create neutral-coverage probes from intervals."""
    cn_rows = [(chrom, start, end, name, 0, 0, 0, 0)
               for chrom, start, end, name in ngfrills.parse_regions(bed_fname)]
    return CNA.from_rows(core.fbase(bed_fname), cn_rows,
                         extra_keys=('gc', 'rmask', 'spread'))


def combine_probes(filenames, has_genome, is_male_normal):
    """Calculate the median coverage of each probe across multiple samples.

    Input:
        List of .cnn files, as generated by 'coverage' or 'import-picard'.
        `has_genome`: reserve columns for GC and RepeatMasker genomic values.
    Returns:
        A single CopyNumArray summarizing the coverages of the input samples,
        including each probe's "average" coverage, "spread" of coverages, and
        genomic GC content.
    """
    # Load coverage from target/antitarget files
    echo("Loading target", filenames[0])
    cnarr1 = CNA.read(filenames[0])

    # Make the sex-chromosome coverages of male and female samples compatible
    chr_x = core.guess_chr_x(cnarr1)
    chr_y = ('chrY' if chr_x.startswith('chr') else 'Y')
    if is_male_normal:
        def shift_sex_chroms(pset):
            """Shift sample X and Y chromosomes for a male-normal reference.

            If sample is male, do nothing.
            If sample is female, chrX -=1, set chrY = -1.
            """
            is_xx = core.guess_xx(pset, chr_x=chr_x)
            if is_xx:
                # Two copies of chrX; halve it
                pset['coverage'][pset.chromosome == chr_x] -= 1.0
                # No chrY; it's all noise, so replace with 1 "flat" copy
                pset['coverage'][pset.chromosome == chr_y] = -1.0

    else:
        def shift_sex_chroms(pset):
            """Shift sample X and Y chromosomes for a female-normal reference.

            If sample is male, chrX += 1.
            If sample is female, set chrY = -1.
            """
            is_xx = core.guess_xx(pset, chr_x=chr_x)
            if is_xx:
                # No chrY; it's all noise, so replace with 1 "flat" copy
                pset['coverage'][pset.chromosome == chr_y] = -1.0
            else:
                # One copy of chrX; double it
                pset['coverage'][pset.chromosome == chr_x] += 1.0

    shift_sex_chroms(cnarr1)
    all_coverages = [cnarr1.coverage]
    for fname in filenames[1:]:
        echo("Loading target", fname)
        pset = CNA.read(fname)
        # Probe information should match across all files
        if not (len(cnarr1) == len(pset)
                and (cnarr1.chromosome == pset.chromosome).all()
                and (cnarr1.start == pset.start).all()
                and (cnarr1.end == pset.end).all()
                and (cnarr1.gene == pset.gene).all()):
            raise RuntimeError("%s probes do not match those in %s"
                               % (fname, filenames[0]))
        shift_sex_chroms(pset)
        all_coverages.append(pset.coverage)
    all_coverages = numpy.vstack(all_coverages)

    echo("Calculating average probe coverages")
    cvg_centers = numpy.apply_along_axis(metrics.biweight_location, 0,
                                         all_coverages)
    echo("Calculating probe spreads")
    spreads = numpy.apply_along_axis(metrics.biweight_midvariance, 0,
                                     all_coverages)
    kwargs = {'spread': spreads}
    # ENH: use the genome FASTA here directly, instead of has_genome?
    if 'gc' in cnarr1:
        # Reuse .cnn GC values if they're already stored (via import-picard)
        kwargs['gc'] = cnarr1['gc']
    elif has_genome:
        # Reserve space for GC values (calculated later)
        kwargs['gc'] = numpy.zeros(len(cnarr1), dtype=numpy.float_)
    if has_genome:
        # Reserve space for RepeatMasker values (calculated later)
        kwargs['rmask'] = numpy.zeros(len(cnarr1), dtype=numpy.float_)
    return CNA("reference", cnarr1.chromosome, cnarr1.start, cnarr1.end,
               cnarr1.gene, cvg_centers, **kwargs)


def warn_bad_probes(probes):
    """Warn about target probes where coverage is poor.

    Prints a formatted table to stderr.
    """
    bad_probes = probes[mask_bad_probes(probes)]
    fg_index = (bad_probes['gene'] != 'Background')
    fg_bad_probes = bad_probes[fg_index]
    if len(fg_bad_probes) > 0:
        # ENH: print coverage and spread too
        bad_pct = 100 * len(fg_bad_probes) / sum(probes['gene'] != 'Background')
        echo("*WARNING*", len(fg_bad_probes), "targets",
             "(%.4f)" % bad_pct + '%', "failed filters:")
        gene_cols = max(map(len, fg_bad_probes['gene']))
        labels = list(map(row2label, fg_bad_probes))
        chrom_cols = max(map(len, labels))
        last_gene = None
        for label, probe in zip(labels, fg_bad_probes):
            if probe['gene'] == last_gene:
                gene = '  "'
            else:
                gene = probe['gene']
                last_gene = gene
            if 'rmask' in probes:
                print("  %s  %s  coverage=%.3f  spread=%.3f  rmask=%.3f"
                      % (gene.ljust(gene_cols), label.ljust(chrom_cols),
                         probe['coverage'], probe['spread'], probe['rmask']))
            else:
                print("  %s  %s  coverage=%.3f  spread=%.3f"
                      % (gene.ljust(gene_cols), label.ljust(chrom_cols),
                         probe['coverage'], probe['spread']))

    # Count the number of BG probes dropped, too (names are all "Background")
    bg_bad_probes = bad_probes[True - fg_index]
    if len(bg_bad_probes) > 0:
        bad_pct = 100 * len(bg_bad_probes) / sum(probes['gene'] == 'Background')
        echo("Antitargets:", len(bg_bad_probes), "(%.4f)" % bad_pct + '%',
             "failed filters")


def mask_bad_probes(probes):
    """Flag the probes with excessively low or inconsistent coverage.

    Returns a bool array where True indicates probes that failed the checks.
    """
    mask = ((probes['coverage'] < params.MIN_BIN_COVERAGE) |
            (probes['spread'] > params.MAX_BIN_SPREAD))
    if 'rmask' in probes:
        mask |= (probes['rmask'] > params.MAX_REPEAT_FRACTION)
    return mask


def get_fasta_stats(probes, fa_fname):
    """Calculate GC and RepeatMasker content of each bin in the FASTA genome."""
    ngfrills.ensure_fasta_index(fa_fname)
    fa_coords = zip(probes.chromosome, probes.start, probes.end)
    echo("Calculating GC and RepeatMasker content in", fa_fname, "...")
    gc_rm_vals = [calculate_gc_lo(subseq)
                  for subseq in ngfrills.fasta_extract_regions(fa_fname,
                                                               fa_coords)]
    gc_vals, rm_vals = zip(*gc_rm_vals)
    return numpy.asfarray(gc_vals), numpy.asfarray(rm_vals)


def calculate_gc_lo(subseq):
    """Calculate the GC and lowercase (RepeatMasked) content of a string."""
    cnt_at_lo = subseq.count('a') + subseq.count('t')
    cnt_at_up = subseq.count('A') + subseq.count('T')
    cnt_gc_lo = subseq.count('g') + subseq.count('c')
    cnt_gc_up = subseq.count('G') + subseq.count('C')
    tot = float(cnt_gc_up + cnt_gc_lo + cnt_at_up + cnt_at_lo)
    if not tot:
        return 0.0, 0.0
    frac_gc = (cnt_gc_lo + cnt_gc_up) / tot
    frac_lo = (cnt_at_lo + cnt_gc_lo) / tot
    return frac_gc, frac_lo


def reference2regions(reference, coord_only=False):
    """Extract iterables of target and antitarget regions from a reference CNA.

    Like loading two BED files with ngfrills.parse_regions.
    """
    cna2rows = (_cna2coords if coord_only else _cna2regions)
    return map(cna2rows, _ref_split_targets(reference))


def _cna2coords(cnarr):
    """Extract the coordinate columns from a CopyNumberArray"""
    return zip(cnarr['chromosome'], cnarr['start'], cnarr['end'])


def _cna2regions(cnarr):
    """Extract the region columns (including genes) from a CopyNumberArray"""
    return zip(cnarr['chromosome'], cnarr['start'], cnarr['end'], cnarr['gene'])


def _ref_split_targets(ref_arr):
    """Split reference into 2 sub-arrays of targets/antitargets."""
    is_bg = (ref_arr.gene == 'Background')
    targets = ref_arr[True - is_bg]
    antitargets = ref_arr[is_bg]
    return targets, antitargets