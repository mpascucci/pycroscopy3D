import os

def create_folders(path):
    # create output folder if it does not exist
    if (path != '') and (not os.path.isdir(path)):
        print(f'directory "{path}" does not exist and will be created.')
        os.makedirs(path, exist_ok=True)