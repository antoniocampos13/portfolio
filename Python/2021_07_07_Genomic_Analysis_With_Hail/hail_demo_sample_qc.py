#
import hail as hl
from hail_demo_init import DEFAULT_REF

from gnomad.utils.filtering import filter_to_autosomes
from gnomad.utils.annotations import add_variant_type

from bokeh.embed import file_html
from bokeh.resources import CDN

# Magic numbers
CALL_RATE = 0.90  # Hail team value = 0.97. GNOMAD value = 0.99
RELATEDNESS = 0.088  # Hail team value. GNOMAD value = 0.08838835
READ_DEPTH = 20  # Hail team value.
FREEMIX_CONTAMINATION = 0.05  # GNOMAD value.
CHIMERIC_READS = 0.05  # GNOMAD value.
MEDIAN_LENGTH = 250  # GNOMAD value.

# Read matrix table
mt = hl.read_matrix_table("recalibrated.mt")

# Show first few samples
mt.s.show()

# Check matrix table fields
mt.describe()

# Mixture of non-empty with empty PL fields causes problems with sample QC for some reason; setting field to all empty
mt = mt.annotate_entries(PL=hl.missing(mt.PL.dtype))

# Add variant-level annotations necessary for variant QC later
## Annotate variants in one of the categories: SNV, multi-SNV, indel, multi-indel, mixed
mt = mt.annotate_rows(**add_variant_type(mt.alleles))

## Number of alleles at the site
mt = mt.annotate_rows(n_alleles = hl.len(mt.alleles))

## Mixed sites (SNVs and indels present at the site)
mt = mt.annotate_rows(mixed_site = hl.if_else(mt.variant_type == "mixed", True, False))

## Spanning deletions
mt = mt.annotate_rows(spanning_deletion=hl.any(lambda a: a == "*", mt.alleles))

# Number of Rows, Columns
mt.count()

# Number of Columns
mt.count_cols()

# Variants breakdown
hl.summarize_variants(mt)

# Split variants with multiple alleles into biallelic configuration. Notice that hl.count() and hl.summarize_variants() will give different numbers after multi-allele sites splitting than before
mt = hl.split_multi_hts(mt)

# Remove monomorphic sites
mt = mt.filter_rows(mt.n_alleles > 1)

# Import sample annotations/metadata
sa = hl.import_table("sample_info.txt", impute=True, key="s")

mt = mt.annotate_cols(sample_info=sa[mt.s])

# Calculate descriptive statistics for sample QC
mt = hl.sample_qc(mt)

# Plot call rate versus mean read depth (coverage)
p = hl.plot.scatter(
    x=mt.sample_qc.dp_stats.mean,
    y=mt.sample_qc.call_rate,
    xlabel="Mean DP",
    ylabel="Call Rate",
    hover_fields={"ID": mt.s},
    size=8,
)

html = file_html(p, CDN, "Chart")

with open("1 Mean Call Rate by Mean DP.html", "w") as f:
    f.write(html)

# Filter by call rate
mt = mt.filter_cols(mt.sample_qc.call_rate >= CALL_RATE)

# Filter by read depth (DP)
mt = mt.filter_cols(mt.sample_qc.dp_stats.mean >= READ_DEPTH)

# Preparing for PCA
for_pca = filter_to_autosomes(mt)
for_pca = for_pca.filter_rows(for_pca.n_alleles == 2)

# Performing the PCA
sample_num = for_pca.cols().count()

_, scores, _ = hl.hwe_normalized_pca(
    for_pca.GT, k=max(1, min(sample_num // 3, 10)), compute_loadings=False
)

relatedness_ht = hl.pc_relate(
    for_pca.GT,
    min_individual_maf=0.01,
    scores_expr=scores[for_pca.col_key].scores,
    block_size=4096,
    min_kinship=0.05,
    statistics="kin",
)

pairs = relatedness_ht.filter(relatedness_ht["kin"] > RELATEDNESS)

related_samples_to_remove = hl.maximal_independent_set(pairs.i, pairs.j, False)

mt = mt.filter_cols(hl.is_defined(related_samples_to_remove[mt.col_key]), keep=False)

# Wrapping up: saving relatednsess Table and dataset MatrixTable to disk
relatedness_ht.write("relatedness.ht", overwrite=True)

mt.write("sampleqc_pass.mt", overwrite=True)