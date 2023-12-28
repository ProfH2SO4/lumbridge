import logging, inspect
from logging.config import dictConfig

LOGGER_NAME: str = "default"


def set_up_logger(log_config: dict) -> None:
    """Set up the logger."""
    dictConfig(log_config)
    logger = logging.getLogger(LOGGER_NAME)
    logger.info("Logger, Status: Ready ")

    return None


def get_logger() -> logging.Logger:
    """Get the logger."""
    return logging.getLogger(LOGGER_NAME)


def info(msg: str) -> None:
    """Log an info message."""
    function_before: str = inspect.stack()[1].function
    msg = f"function: {function_before} | {msg}"
    get_logger().info(msg)


def debug(msg: str) -> None:
    """Log a debug message."""
    function_before: str = inspect.stack()[1].function
    msg = f"function: {function_before} | {msg}"
    get_logger().debug(msg)


def warning(msg: str) -> None:
    """Log a warning message."""
    function_before: str = inspect.stack()[1].function
    msg = f"function: {function_before} | {msg}"
    get_logger().warning(msg)


def error(msg: str) -> None:
    """Log an error message."""
    function_before: str = inspect.stack()[1].function
    msg = f"function: {function_before} | {msg}"
    get_logger().error(msg)
