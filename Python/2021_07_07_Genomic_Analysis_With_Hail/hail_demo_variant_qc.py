# Import modules and init Hail
import hail as hl

from gnomad.utils.annotations import bi_allelic_site_inbreeding_expr
from gnomad.variant_qc.random_forest import apply_rf_model, median_impute_features
from gnomad.variant_qc.pipeline import train_rf_model
from hail_init import DEFAULT_REF

# Variant Quality hard filters
INBR_COEFF = -0.3
AB_LOWER_LIM = 0.2
AB_UPPER_LIM = 1 - AB_LOWER_LIM

# Read MatrixTable with sample QC-passing dataset
mt = hl.read_matrix_table("sampleqc_pass.mt")

# Calculate variant statistics
mt = hl.variant_qc(mt)

# Calculate inbreeding coefficient
mt = mt.annotate_rows(inbr_coeff=bi_allelic_site_inbreeding_expr(mt.GT))

# Determine the maximum p-value for sampling the observed allele balance under a binomial model
mt = mt.annotate_rows(
    pab_max=hl.agg.max(
        hl.binom_test(mt.AD[1], mt.DP, 0.5, "two-sided")
    )
)

# Removing variants with excess of heterozygotes
mt = mt.filter_rows(mt.inbr_coeff > INBR_COEFF)

# Removing variants for which no sample had high quality genotypes
mt = mt.filter_rows(hl.agg.any(mt.GQ >= 20))
mt = mt.filter_rows(hl.agg.any(mt.DP >= 10))

mt = mt.annotate_entries(AB=(mt.AD[1] / hl.sum(mt.AD)))

mt = mt.filter_rows(
    hl.agg.any(
        (mt.GT.is_hom_ref() & (mt.AB < AB_LOWER_LIM))
        | (mt.GT.is_het() & (mt.AB >= AB_LOWER_LIM) & (mt.AB <= AB_UPPER_LIM))
        | (mt.GT.is_hom_var() & (mt.AB > AB_UPPER_LIM))
    )
)

# For random forest model training - True positives: PASSing VQSR filter
mt = mt.annotate_rows(tp=hl.if_else(hl.len(mt.filters) == 0, True, False))

# For random forest model training - False positives: not PASSing VQSR
mt = mt.annotate_rows(fp=hl.if_else(hl.len(mt.filters) != 0, True, False))

rf_ht = mt.select_rows(
    mt.inbr_coeff,
    mt.info.SOR,
    mt.info.ReadPosRankSum,
    mt.info.MQRankSum,
    mt.info.QD,
    mt.pab_max,
    mt.variant_type,
    mt.n_alleles,
    mt.mixed_site,
    mt.spanning_deletion,
    mt.tp,
    mt.fp,
).rows()

# List of features for the random forest model training
features = [
    "inbr_coeff",
    "SOR",
    "ReadPosRankSum",
    "MQRankSum",
    "QD",
    "pab_max",
    "variant_type",
    "n_alleles",
    "mixed_site",
    "spanning_deletion",
]

# Impute median values of annotation in empty entries
rf_ht = median_impute_features(rf_ht)

# Chr20 variants will be used for testing the model only
test_intervals = ["chr20"]

test_intervals = [
    hl.parse_locus_interval(x, reference_genome="GRCh38") for x in test_intervals
]

# Training the model
rf_trained_ht, rf_model = train_rf_model(
    rf_ht,
    rf_features=features,
    tp_expr=rf_ht.tp,
    fp_expr=rf_ht.fp,
    test_expr=hl.literal(test_intervals).any(
        lambda interval: interval.contains(rf_ht.locus)
    ),
)

# Joining results
ht = rf_ht.join(rf_trained_ht, how="left")

rf_results = apply_rf_model(
    ht=ht,
    rf_model=rf_model,
    features=hl.eval(rf_trained_ht.features),
    label="rf_label",
    prediction_col_name="rf_prediction",
)

rf_summary_ht = rf_results.group_by(
    "tp", "fp", "rf_train", "rf_label", "rf_prediction"
).aggregate(n=hl.agg.count())

# Join with the sample_qc passing matrix table
variantqc_pass = mt.annotate_rows(**rf_results[mt.locus, mt.alleles])

variantqc_pass.write("variantqc_pass.mt", overwrite=True)