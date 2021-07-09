# %% Load modules
import hail as hl

DEFAULT_REF = "GRCh38"

hl.init(
    idempotent=True,
    quiet=True,
    skip_logging_configuration=True,
    default_reference=DEFAULT_REF,
    spark_conf={
        "spark.executor.cores": "4",
        "spark.driver.memory": "16g",
    },
)