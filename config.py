INPUT_FASTA = "./test_data/fasta_folder"
INPUT_GFF3 = "./test_data/gff3_folder"
INPUT_ORF = "./test_data/orf_folder"
OUTPUT_FOLDER = "./lumbridge_output"

# features to extract from gff3
GFF3_FEATURES = ["miRNA", "rRNA", "gene"]

# Homer2
HOMER2_BIN_PATH = "/usr/src/app/bin"
HOMER2_OUTPUT_FOLDER = "./lumbridge_output/homer2"
HOMER2_UPSTREAM_GEN_SEQ_LENGTH = 800
HOMER2_DOWNSTREAM_GEN_SEQ_LENGTH = 200
HOMER2_P_THRESHOLD = 0.05
HOMER2_CPU_CORES = 4


# model_data
MAX_FEATURE_OVERLAP = 1  # for given feature, fow now it's tested only for 1


LOG_CONFIG = {
    "version": 1,
    "formatters": {
        "default": {
            "format": "LUMBRIDGE - %(asctime)s - %(levelname)s - %(message)s",
            "datefmt": "%Y-%m-%d %H:%M:%S",
        }
    },
    "handlers": {
        "console": {
            "level": "DEBUG",
            "class": "logging.StreamHandler",
            "formatter": "default",
        },
        "file": {
            "level": "DEBUG",
            "class": "logging.FileHandler",
            "formatter": "default",
            "filename": "./logs/default.txt",  # Specify the file path
        },
    },
    "loggers": {
        "default": {
            "level": "DEBUG",
            "handlers": ["console", "file"],  # Updated to use console handler
            "propagate": False,
        }
    },
    "disable_existing_loggers": False,
}
