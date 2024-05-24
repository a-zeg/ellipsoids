
import json
import visualisation
import barcodePlotting
import matplotlib.pyplot as plt
import matplotlib as mpl

filenames = ['debug_dict__20240522_105148.json',
             'debug_dict__20240522_105207.json',
             'debug_dict__20240522_105214.json']

for filename in filenames:
    with open(filename, "r") as f:
        jsonVars = json.load(f)

    data = jsonVars['data_test_trnsfs_std']
    labels = jsonVars['labels_test']

    print(filename)
    print(len(data))
    print(len(labels))
    print('\n')


filename_crash = 'debug_dict__20240522_121131.json'
with open(filename, "r") as f:
    jsonVars_crash = json.load(f)

data_crash = jsonVars['data_test_trnsfs_std']
labels_crash = jsonVars['labels_test']

print(filename_crash)
print(len(data_crash))
print(len(labels_crash))
print('\n')

print('normal barcode lengths:')
for barcode in data:
    print(len(barcode))

print('\n')

print('crash barcode lengths:')
for barcode in data_crash:
    print(len(barcode))


# data_plot = []
# for barcode in data:
#     data_plot.append([[1, bar] for bar in barcode])
# print(data_plot[0])

# cmap = plt.get_cmap('viridis')
# num_colors = len(data_plot)
# colors = [cmap(i / num_colors) for i in range(num_colors)]

# ax = plt.axes()
# # fig = plt.figure()
# for i, barcode in enumerate(data_plot):
#     barcodePlotting.plot_persistence_barcode(barcode, inf_delta=0.5, fontsize=12, axes=ax,
#                                             axis_start = -0.1, max_intervals=100, colormap=[None, colors[i]]) #(0.1 + xAxisEnd))
# plt.show()






