# Import modules and init Hail
import hail as hl
from hail_init import DEFAULT_REF

# Read variant QC passing matrix table
variantqc_pass = hl.read_matrix_table("variantqc_pass.mt")

# Exact or approximate coordinates
intervals = ["chr10:52765380-52772784", "chr1:100M-200M"]

filtered_mt = hl.filter_intervals(
    variantqc_pass,
    [hl.parse_locus_interval(x, reference_genome=DEFAULT_REF) for x in intervals]) 

# Nucleotide window around locus 
locus = hl.parse_locus("chrX:23833353", DEFAULT_REF)
window = locus.window(100000, 100000) # 100,000 nucleotides before and after

filtered_mt = variantqc_pass.filter_rows(window.contains(variantqc_pass.locus))

# Filter by allelic frequency
filtered_mt = filtered_mt.filter_rows(filtered_mt.variant_qc.AF[1] < 0.01)
