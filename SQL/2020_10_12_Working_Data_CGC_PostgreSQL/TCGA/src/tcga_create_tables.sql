CREATE TABLE allcases (
    case_id TEXT,
    case_primarysite TEXT,
    diagnosis TEXT,
    diagnosis_ageatdiagnosis_1 INT,
    diagnosis_clinicalt_1 TEXT
);

CREATE TABLE demographic (
    case_id TEXT,
    case_primarysite TEXT,
    demographic TEXT,
    demographic_ethnicity_1 TEXT,
    demographic_race_1 TEXT
);

CREATE TABLE follow_up (
    case_id TEXT,
    case_primarysite TEXT,
    followup TEXT,
    followup_primarytherapyoutcomesuccess_1 TEXT
);

CREATE TABLE allfiles (
    case_id TEXT,
    case_primarysite TEXT,
    file_uuid TEXT,
    file_accesslevel_1 TEXT,
    file_datatype_1 TEXT,
    filename TEXT
);
