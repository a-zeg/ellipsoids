
# from src.ellipsoids_vs_rips import labelling

# def test_labelling():
#     meshlist = \
#     [
#         './cat1.mat_nPts=300_seed=0.mat', 
#         './centaur0.mat_nPts=300_seed=0.mat', 
#         './victoria0.mat_nPts=300_seed=0.mat',
#         './cat2.mat_nPts=300_seed=0.mat', 
#         './cat0.mat_nPts=300_seed=0.mat', 
#         './cat3.mat_nPts=300_seed=0.mat', 
#         './centaur1.mat_nPts=300_seed=0.mat'
#     ]

#     correctDict = {'cat': 0, 'centaur': 1, 'victoria': 2}

#     labelDict = labelling(meshlist)

#     assert correctDict == labelDict


from src.main import generate_filename
from src.data_handling import get_timestamp



def test_generate_filename():
    folder = 'data'
    timestamp = get_timestamp()

    dict = {
        'param1': 1,
        'param2': 3.4,
        'param3': 'string'
    }
    generated_filename \
        = generate_filename(filename_parameters=dict, folder=folder, timestamp=timestamp)
    correct_filename \
        = 'data/ellipsoids_param1=1_param2=3.4_param3=string_' + timestamp    

    assert correct_filename == generated_filename

from src.main import filter_dictionary

def test_filter_dictionary():
    var1 = 1
    var2 = 2.0
    var3 = 'string'

    dict_all_vars = locals()

    vars_to_save = ['var1', 'var3']

    correct_dictionary = {
        'var1': 1,
        'var3': 'string'
    }

    filtered_dictionary = filter_dictionary(vars_to_save,dict_all_vars)

    assert correct_dictionary == filtered_dictionary

# if __name__ == "__main__":

#     # test_generate_filename()
#     test_labelling()



