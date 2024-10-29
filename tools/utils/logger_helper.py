"""
Logger helper module
"""

from typing import Optional, Union
import logging
import inspect
import sys
from pathlib import Path

from omegaconf import OmegaConf

from . import files
from .str_helper import to_title_case


RESET = "\033[0m"
BRIGHT_RED = "\033[91m"
CRITICAL_RED = "\033[41m"
BRIGHT_YELLOW = "\033[93m"
BRIGHT_BLUE = "\033[94m"
BRIGHT_INFO = "\033[96m"
GREY = "\033[90m"

TRACE_LEVEL_NUM = 5
TRACE_LEVEL_NAME = "TRACE"


def add_trace_level():
    """
    Add the TRACE level to the logging module.
    """
    logging.addLevelName(TRACE_LEVEL_NUM, TRACE_LEVEL_NAME)

    def trace(self, message, *args, **kwargs):
        if self.isEnabledFor(TRACE_LEVEL_NUM):
            self._log(TRACE_LEVEL_NUM, message, args, **kwargs)

    logging.Logger.trace = trace


add_trace_level()


class LoggerFormatter(logging.Formatter):
    def __init__(self, fmt=None, datefmt=None):
        super().__init__(fmt, datefmt)

    def format(self, record: logging.LogRecord) -> str:
        """
        Format the log record with the output name in title case and color based on log level.

        Args:
            record: Log record

        Returns:
            Formatted log record
        """
        # Add a custom field for title-cased logger name
        record.output_name = to_title_case(record.name)

        # Set color based on log level
        color = self.get_log_colour(record.levelno)

        formatted_message = super().format(record)
        return f"{color}{formatted_message}{RESET}"

    def get_log_colour(self, level: int) -> str:
        """
        Get the log colour based on the log level.

        Args:
            level: Log level

        Returns:
            Log colour
        """
        # Custom log levels
        if level == TRACE_LEVEL_NUM:
            return GREY

        # Standard log levels
        match level:
            case logging.DEBUG:
                return BRIGHT_BLUE
            case logging.INFO:
                return BRIGHT_INFO
            case logging.WARNING:
                return BRIGHT_YELLOW
            case logging.ERROR:
                return BRIGHT_RED
            case logging.CRITICAL:
                return CRITICAL_RED
            case _:
                return RESET


def get_logger_config() -> OmegaConf:
    """
    Get the logger configuration.

    Returns:
        Logger configuration
    """
    configs_folder = files.get_configs_folder()
    files.create_folder_if_not_exists(configs_folder)
    config_path = configs_folder / "logger.yaml"
    if not config_path.exists():
        raise FileNotFoundError(
            f"Logger configuration file not found: {config_path}")

    logger_config = OmegaConf.load(config_path)
    return logger_config


def get_caller_module_name(caller_frame: inspect.FrameInfo) -> str:
    """
    Get the full name of the caller's module.

    Args:
        caller_frame: Caller's frame

    Returns:
        Caller's module name
    """
    caller_module = inspect.getmodule(caller_frame[0])
    caller_name = caller_module.__name__

    if caller_name.startswith("src"):
        # Running module standalone
        # Determine module from path
        relative_path = Path(caller_module.__file__)
        rel_path_str = relative_path.as_posix()
        caller_relative_path = caller_name.replace(".", "/") + ".py"
        if rel_path_str.endswith(caller_relative_path):
            rel_path_str = rel_path_str[: -len(caller_relative_path)]

            if rel_path_str.endswith("/"):
                rel_path_str = rel_path_str[:-1]

        caller_name = rel_path_str.split("/")[-1]

    return caller_name


def map_log_level(level: str) -> int:
    """
    Map the log level string to the corresponding logging level.

    Args:
        level: Log level string

    Returns:
        Logging level
    """
    # See https://docs.python.org/3/library/logging.html#logging.getLevelNamesMapping

    if sys.version_info >= (3, 11):
        # Introduced in Python 3.11
        level_mapping = logging.getLevelNamesMapping()
    else:
        # Python 3.10 and below
        level_mapping = logging._nameToLevel

    log_level = level_mapping.get(level.upper())
    if log_level is None:
        raise ValueError(f"Invalid logging level: {level}")

    return log_level


def get_log_level(level: Union[int, str]) -> Optional[int]:
    """
    Get the logging level based on the caller's module name.

    Args:
        level: Logging level

    Returns:
        Logging level
    """

    loggers = logger_config.loggers

    caller_frame = inspect.stack()[2]
    caller_name = get_caller_module_name(caller_frame)
    top_level_module = caller_name.split(".")[0]

    log_level = None
    if top_level_module in loggers:
        # Override the logging level with the specified level
        log_level = loggers[top_level_module]
        if log_level == False:
            # Disable the logger
            log_level = logging.CRITICAL + 1
        else:
            # Parse the logging level
            log_level = map_log_level(log_level)

    if log_level is None:
        if type(level) == int:
            log_level = level
        elif type(level) == str:
            log_level = map_log_level(level)

    return log_level


def init_logger(level: Union[int, str] = logging.INFO) -> logging.Logger:
    """
    Initialise a named logger with the specified logging level.

    Args:
        level: Logging level

    Returns:
        Logger instance
    """

    # [1] gives the caller of this function
    caller_frame = inspect.stack()[1]
    caller_name = get_caller_module_name(caller_frame)

    # Use the caller's module name for the logger
    if caller_name == "__main__":
        logger_name = "Global"
    elif caller_name == "__mp_main__":
        logger_name = "Multiprocessing"
    else:
        logger_name = caller_name

    logger = logging.getLogger(logger_name)

    level = get_log_level(level)
    logger.setLevel(level)

    logger.propagate = False

    if not logger.hasHandlers():
        attach_formatter(logger)

    return logger


def attach_formatter(logger: logging.Logger) -> None:
    """
    Attach a formatter to the logger.

    Args:
        logger: Logger instance

    Returns:
        None
    """

    formatter = LoggerFormatter(
        "%(asctime)s  %(output_name)-35s %(levelname)-13s%(message)s")

    if not logger.handlers:
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)
    else:
        for handler in logger.handlers:
            handler.setFormatter(formatter)


def disable_logger(logger_name: str) -> None:
    """
    Disable a logger and all its handlers.

    Args:
        logger_name: Logger name

    Returns:
        None
    """

    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.CRITICAL + 1)


def test_logger():
    """
    Test the logger helper functions
    """

    logger = init_logger(TRACE_LEVEL_NUM)

    logger.trace("This is a trace message")
    logger.debug("This is a debug message")
    logger.info("This is an info message")
    logger.warning("This is a warning message")
    logger.error("This is an error message")
    logger.critical("This is a critical message")


if __name__ == "__main__":
    test_logger()

logger_config = get_logger_config()
