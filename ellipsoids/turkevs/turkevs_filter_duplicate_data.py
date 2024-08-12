'''
Removes duplicate files with the same id.
'''

import os
import re
import shutil

parentfolder = 'data/turkevs_20240520'

id = 'id=0004'

# data_folder = find_subfolder_with_given_id(parentfolder, id)
data_folder = os.path.join(parentfolder, id)
destinationfolder = os.path.join(parentfolder, id + '_filtered')
if not os.path.isdir(destinationfolder):
    os.makedirs(destinationfolder)
    print('Created folder ' + destinationfolder)

filenames = [f for f in os.listdir(data_folder) if os.path.isfile(os.path.join(data_folder, f)) if f.endswith('.json')]

print('Total number of files found: ' + str(len(filenames)))


trnsfs = ["std", "trns", "rot", "stretch", "shear", "gauss", "out"]

index_dict = {}
for trnsf in trnsfs:
    index_dict[trnsf] = list(range(200))



filtered_filenames = []
for file in filenames:
    trnsf_current = re.search(r'Turkevs-([a-z]+)', file).group(1)
    index_current = int(re.search(r'Turkevs-[a-z]+-(\d+)', file).group(1))

    if index_current in index_dict[trnsf_current]:
        index_dict[trnsf_current].remove(index_current)
        filtered_filenames.append(file)

if all(index_dict[trnsf] == [] for trnsf in trnsfs):
    print('Success')

for file in filtered_filenames:
    shutil.copy(os.path.join(data_folder,file), destinationfolder)

