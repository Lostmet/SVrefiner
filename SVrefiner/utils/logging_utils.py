import sys
import logging
import datetime
from pathlib import Path

captured_logs = []

class WarningErrorCaptureHandler(logging.Handler):
    def emit(self, record):
        if record.levelno >= logging.WARNING: 
            msg = self.format(record)
            captured_logs.append(msg)

def log_all_warnings_and_errors():
    logger = logging.getLogger()
    if captured_logs:
        logger.info("")
        logger.info("Summary of Warnings and Errors:")
        for msg in captured_logs:
            logger.info(f"  - {msg}")

def log_tqdm_summary(pbar, logger):
    d = pbar.format_dict
    desc = pbar.desc or "Task"
    minutes = int(d["elapsed"] // 60)
    seconds = int(d["elapsed"] % 60)
    logger.info(f"{desc} completed: {d['n']} {d['unit']} processed in {minutes}m {seconds}s.")

class LazyErrorFileHandler(logging.Handler):
    """
    Lazy error file handler:
    Creates the pipeline_error.log file only when the first WARNING or higher record is emitted.
    """
    def __init__(self, log_file: Path):
        super().__init__(level=logging.WARNING)
        self.log_file = Path(log_file)
        self.formatter = logging.Formatter('%(message)s')
        self.file_handler = None

    def emit(self, record):
        if self.file_handler is None:
            self.log_file.parent.mkdir(parents=True, exist_ok=True)
            self.file_handler = logging.FileHandler(self.log_file)
            self.file_handler.setLevel(logging.WARNING)
            self.file_handler.setFormatter(self.formatter)
        self.file_handler.emit(record)

def setup_logging(log_file: Path):
    """
    Setup logging with:
    - console output
    - main log file (INFO level)
    - lazy error log file (WARNING+), only created if needed
    - internal capture handler for warnings/errors
    """
    log_file = Path(log_file)
    log_file.parent.mkdir(parents=True, exist_ok=True)

    logger = logging.getLogger()
    if logger.hasHandlers():
        logger.handlers.clear()
    logger.setLevel(logging.INFO)

    formatter = logging.Formatter('%(message)s')

    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    # Main log file handler
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    # Lazy error file handler
    error_log_file = log_file.with_name("pipeline_error.log")
    lazy_error_handler = LazyErrorFileHandler(error_log_file)
    logger.addHandler(lazy_error_handler)

    # Internal warning/error capture handler
    capture_handler = WarningErrorCaptureHandler()
    capture_handler.setLevel(logging.WARNING)
    capture_handler.setFormatter(formatter)
    logger.addHandler(capture_handler)

def get_logger(name: str) -> logging.Logger:
    """Get a logger instance"""
    return logging.getLogger(name)

def log_step(title: str, width: int = 125, style: str = "box"):
    logger = logging.getLogger()
    if style == "box":
        text = f"[ {title} ]"
        side = (width - len(text)) // 2
        line = "═" * side + text + "═" * (width - len(text) - side)
        logger.info(line)
    elif style == "flat":
        logger.info("=" * width)
        logger.info(f"[ {title} ]")
        logger.info("-" * width)

def get_clean_command() -> str:
    program = Path(sys.argv[0]).name
    args = sys.argv[1:]

    if not args:
        return f"Command:\n  {program}"

    grouped = []
    i = 0
    while i < len(args):
        if args[i].startswith("-"):
            if i + 1 < len(args) and not args[i + 1].startswith("-"):
                grouped.append((args[i], args[i + 1]))
                i += 2
            else:
                grouped.append((args[i], None))
                i += 1
        else:
            grouped.append((args[i], None))
            i += 1

    lines = [f"  {program} {args[0]}"] 
    max_flag_len = max(len(flag) for flag, val in grouped[1:]) if len(grouped) > 1 else 0

    for flag, val in grouped[1:]:
        if val is not None:
            lines.append(f"  {flag.ljust(max_flag_len)}   {val} ")
        else:
            lines.append(f"  {flag}")

    return "\n".join(lines)


def log_summary_block(cmd: str, start: float, duration: float, stats: dict):

    end_time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    runtime_str = f"{int(duration // 3600)}:{int(duration % 3600 // 60):02d}:{int(duration % 60):02d}"
    start_time_str = datetime.datetime.fromtimestamp(start).strftime('%Y-%m-%d %H:%M:%S')
    logger = get_logger(__name__)
    logger.info(f"{'Command:':<5}{cmd}")
    logger.info(f"{'Start time:':<30}{start_time_str}")
    logger.info(f"{'End time:':<30}{end_time}")
    logger.info(f"{'Total runtime:':<30}{runtime_str}")
    logger.info('')

    for k, v in stats.items():
        logger.info(f"{k + ':':<30}{v}")

import multiprocessing
import logging

def log_thread_info(user_threads: int):
    logger = logging.getLogger()
    try:
        total_cores = multiprocessing.cpu_count()
    except NotImplementedError:
        total_cores = 'unknown'

    logger.info(f"Detected {total_cores} logical cores. Currently using {user_threads} thread(s). You can modify this via the --threads option.")

    if isinstance(total_cores, int) and user_threads > total_cores:
        logger.warning(f"You requested {user_threads} threads, but only {total_cores} logical cores are available. Consider reducing --threads.")
