#
import hail as hl
from hail_demo_init import DEFAULT_REF

hl.import_vcf("snp.recalibrated.vcf.gz",
force_bgz=True,
reference_genome=DEFAULT_REF,
array_elements_required=False).write("recalibrated.mt",overwrite=True)