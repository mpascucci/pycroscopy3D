import os
import logging

def create_folders(path, logger=None):
    # create output folder if it does not exist
    if (path != '') and (not os.path.isdir(path)):
        if isinstance(logger, logging.Logger):
            logger.info(f'directory "{path}" does not exist and will be created.')
        os.makedirs(path, exist_ok=True)