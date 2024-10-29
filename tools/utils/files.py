"""
File handler
"""

from pathlib import Path
from typing import Optional
import logging

logger = logging.getLogger(__name__)


def get_project_root():
    """Get the project folder."""
    this_file_path = Path(__file__).resolve()
    project_path = this_file_path.parent.parent.parent
    return project_path


def get_data_dir():
    """Get the data folder."""
    data_folder = get_project_root() / "data"
    create_folder_if_not_exists(data_folder)
    return data_folder


def get_models_dir():
    """Get the models folder."""
    models_dir = get_data_dir() / "models"
    create_folder_if_not_exists(models_dir)
    return models_dir


def get_configs_folder():
    """Get the configs folder."""
    configs_folder = get_project_root() / "configs"
    create_folder_if_not_exists(configs_folder)
    return configs_folder


def file_exists(file_path: Path) -> bool:
    """
    Checks if a file exists at the specified path.

    Args:
        file_path (Path): The path to the file to check.

    Returns:
        bool: True if the file exists, False otherwise.
    """

    exists = file_path.exists()
    if exists:
        logger.debug("File exists: %s", file_path)
    else:
        logger.debug("File does not exist: %s", file_path)

    return exists


def create_folder_if_not_exists(folder_path: Path):
    """
    Creates a folder at the specified path if it does not already exist.

    Args:
        folder_path (Path): The path to the folder to create.
    """

    if not folder_path.exists():
        folder_path.mkdir(parents=True)
        logger.info("Folder created: %s", folder_path)
    else:
        logger.debug("Folder already exists: %s", folder_path)


def list_files_in_folder(folder_path: Path, file_types: Optional[list[str]] = None) -> list[Path]:
    """
    Lists all files in the specified folder.

    Args:
        folder_path (Path): The path to the folder to list files from.
        file_type (Optional[list[str]]): A list of file types to filter the files by. Defaults to None.

    Returns:
        list[Path]: A list of paths to the files in the folder.
    """

    if file_types is None:
        file_types = [".*"]

    files = [f for f in folder_path.iterdir() if f.is_file() and any(
        f.name.endswith(file_type) for file_type in file_types)]

    logger.debug("Files in folder: %s - %s", folder_path, files)
    return files


def relative_path(file_path: Path) -> Path:
    """
    Returns the relative path of the file from the project root.

    Args:
        file_path (Path): The path to the file.

    Returns:
        Path: The relative path of the file from the project root.
    """

    project_root = get_project_root()
    relative_path = file_path.relative_to(project_root).as_posix()

    return f"./{relative_path}"
