'''
The script calculate_turkevs.py outputs one json file per point cloud.
To analyse these futher, we need to group them in folders, according to 
the ids.

This script moves all json files in a given folder to folders with names
'id=XXXX', where XXXX is the id corresponding to each file.
'''


import os
import shutil
import data_handling
import re

source = "data/turkevs_20240520"

print('Getting paths of .json files from ' + source, end='', flush=True)
srcpaths = data_handling.get_paths_of_files_in_a_folder(source, '.json')
print('Done.')


ids = []

print('Moving files...')
for path in srcpaths:

    match = re.search(r'(id=\d+)', path)
    if match is None:
        continue
    
    id = match.group(1)
    if id not in ids:
        ids.append(id)
    
    subfolder = os.path.join(source, id)
    if not os.path.isdir(subfolder):
        os.makedirs(subfolder)
        print('Created folder ' + subfolder)

    shutil.move(path, subfolder)

print('Done.')



