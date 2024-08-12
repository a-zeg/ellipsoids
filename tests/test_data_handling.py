from ellipsoids.data_handling import parse_turkevs_filename
from ellipsoids.data_handling import read_pd0_and_pd1

def test_parse_turkevs_filename():
    filename = 'data/turkevs/ellipsoids_data_type=Turkevs-trns-079_n_pts=300_nbhd_size=9_axes_ratios=[3 3 1]__20240427_103930.json'

    transformation, index = parse_turkevs_filename(filename)

    assert transformation == 'trns'
    assert index == 79



def equal_ignore_order(a, b):
    """ Use only when elements are neither hashable nor sortable! """
    print(a)
    print(b)
    unmatched = list(b)
    for element in a:
        try:
            unmatched.remove(element)
        except ValueError:
            return False
    print(unmatched)
    return not unmatched

# def test_read_pd0_and_pd1():
#     path = 'tests/test_data/turkevs/testEllipsoids_data_type=Turkevs-std-094_n_pts=301_nbhd_size=9_axes_ratios=[3 1]__20240427_100454.json'

#     pd0,pd1,points = read_pd0_and_pd1(path)

#     points_test = [
#         [
#             0.04719981494445757,
#             0.0
#         ],
#         [
#             0.2952852157810678,
#             0.0
#         ],
#         [
#             0.10289952945739707,
#             0.0
#         ],
#         [
#             0.045391384515248295,
#             0.0
#         ],
#         [
#             0.12107783748711236,
#             0.0
#         ],
#         [
#             0.2845905593250435,
#             0.0
#         ]
#     ]

#     pd0_test = [[
#                 0.0,
#                 0.029964225883488276
#             ],
#             [
#                 0.0,
#                 0.0019173261391611301
#             ]]
#     pd1_test = [
#             [
#                 0.09224727005924485,
#                 0.5857452128573868
#             ],
#             [
#                 0.09993224543035839,
#                 0.5857452128573868
#             ],
#             [
#                 0.11120559667633173,
#                 0.5255687631696404
#             ]
#         ]   
    
#     assert equal_ignore_order(pd0, pd0_test)
#     assert equal_ignore_order(pd1, pd1_test)
#     assert equal_ignore_order(points, points_test)


# from src.data_handling import import_turkevs_transformed


# def test_dictionary():
#     mydict = {
#         'std': [[0,1],[0,2],[0,3]]
#     }

#     assert mydict['std'][1] == [0,2]


# # def test_generate_filename():
# #     folder = 'data'
# #     timestamp = get_timestamp()

# #     dict = {
# #         'param1': 1,
# #         'param2': 3.4,
# #         'param3': 'string'
# #     }
# #     generated_filename \
# #         = generate_filename(filename_parameters=dict, folder=folder, timestamp=timestamp)
# #     correct_filename \
# #         = 'data/ellipsoids_param1=1_param2=3.4_param3=string_' + timestamp    

# #     assert correct_filename == generated_filename