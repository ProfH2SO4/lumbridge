INPUT_FASTA = "./test_data/fasta_folder"
INPUT_GFF3 = "./test_data/gff3_folder"
INPUT_ORF = "./test_data/orf_folder"
OUTPUT_FOLDER = "./lumbridge_output"

# Homer2
HOMER2_OUTPUT_FOLDER = "./lumbridge_output/homer2"
HOMER2_UPSTREAM_GEN_SEQ_LENGTH = 800
HOMER2_DOWNSTREAM_GEN_SEQ_LENGTH = 200
HOMER2_P_THRESHOLD = 0.05
HOMER2_CPU_CORES = 4


LOG_CONFIG = {
    "version": 1,
    "formatters": {
        "default": {
            "format": "LUMBRIDGE - %(asctime)s - %(levelname)s - %(message)s",
            "datefmt": "%Y-%m-%d %H:%M:%S",
        }
    },
    "handlers": {
        "sys_logger6": {
            "level": "DEBUG",
            "class": "logging.handlers.SysLogHandler",
            "formatter": "default",
            "address": "/dev/log",
            "facility": "local6",
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
            "handlers": ["sys_logger6", "file"],
            "propagate": False,
        }
    },
    "disable_existing_loggers": False,
}
