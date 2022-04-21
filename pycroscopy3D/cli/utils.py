import os

_verbose = False

def  log(s):
    global _verbose
    if _verbose:
        print(s)

def create_folders(path, verbose=True):
    global _verbose
    _verbose = verbose
    # create output folder if it does not exist
    if not os.path.isdir(path):
        log(f'directory "{path}" does not exist and will be created.')
        os.makedirs(path, exist_ok=True)