from pathlib import Path

from rare_event_detection.api import get_REI_from_testing_scan
import os
from pydantic import BaseModel
import logging
from datetime import datetime
import re

DATA_DIR = Path("/anvil/projects/x-cis220161/datasets/aps/")
EXPERIMENT_DIR = DATA_DIR / "experiment-apr4"
DECOMPRESSED_DIR = DATA_DIR / "apr15-logscale-largefiles" / "decompressed_files"

baseline_scan_path = DATA_DIR / "base.edf.ge5"
test_scan_path = DATA_DIR / "test.edf.ge5"
baseline_dark_path = DATA_DIR / "dark_4_base.edf.ge5"
test_dark_path = DATA_DIR / "dark_4_test.edf.ge5"

KMeansModelPath = EXPERIMENT_DIR / "kmeans_model.pkl"
EmbeddingModelIterDir = EXPERIMENT_DIR / "embedding_model_iter"
trained_encoder_path = EmbeddingModelIterDir / "script-ep00100.pth"


class DetectResult(BaseModel):
    filename: str
    REI_score: float
    time_consumed: float


formatter = logging.Formatter("%(asctime)s %(levelname)s %(message)s")


def setup_logger(name, log_file, level=logging.INFO):
    """To setup as many loggers as you want"""

    handler = logging.FileHandler(log_file)
    handler.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)

    return logger


current_time = datetime.now()
log_prefix = (
    EXPERIMENT_DIR
    / "sz3_detect_logs"
    / f"{current_time.month}-{current_time.day}-{current_time.year}_{current_time.hour}-{current_time.minute}-{current_time.second}"
)
log_prefix.mkdir(parents=True, exist_ok=True)

logger = setup_logger("app_logger", log_prefix / "sz3_detect_app.log")
logger.info("The benchmark for rare event detection has started!")

print("Starting the detection benchmark for decompressed files")


def extract_middle_number(filename):
    match = re.search(r"ge5-([\d.]+)-sz3", str(filename))
    return (
        float(match.group(1)) if match else float("inf")
    )  # fallback in case of no match


dp_files = list(DECOMPRESSED_DIR.iterdir())
sorted_dpfiles = sorted(dp_files, key=extract_middle_number)

results = []
for i, file in enumerate(sorted_dpfiles):
    print(f"DP File {i}: {file}")
    REI_score, time_consumed = get_REI_from_testing_scan(
        trained_encoder_path=trained_encoder_path,
        testing_scan_path=file,
        testing_scan_dark_path=test_dark_path,
        kmeans_model_path=KMeansModelPath,
    )
    results.append(
        DetectResult(
            filename=file.name, REI_score=REI_score, time_consumed=time_consumed
        )
    )
    print(
        f"DP File {i} - {file.name} has REI score {REI_score}, and the time consumed is {time_consumed}"
    )
    logger.info(
        f"DP File {i} - {file.name} has REI score {REI_score}, and the time consumed is {time_consumed}"
    )

REI_score, time_consumed = get_REI_from_testing_scan(
    trained_encoder_path=trained_encoder_path,
    testing_scan_path=test_scan_path,
    testing_scan_dark_path=test_dark_path,
    kmeans_model_path=KMeansModelPath,
)

print("The standard REI score is ", REI_score)
logger.info(f"the standard REI score is {REI_score}")

REIs = []
consumed_times = []
for result in results:
    REIs.append(result.REI_score)
    consumed_times.append(result.time_consumed)
    print(
        f"{result.filename} has REI Score {result.REI_score}, and the time consumed is {result.time_consumed}"
    )
    logger.info(
        f"{result.filename} has REI Score {result.REI_score}, and the time consumed is {result.time_consumed}"
    )

print(f"REI Scores: {REIs}")
print(f"Consumed Times: {consumed_times}")
