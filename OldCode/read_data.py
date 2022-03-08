################################################################################
# Michael Cooney
# 01/13/22
################################################################################

import io
import numpy as np
import pandas as pd

# Requirements: Average number of 1s present in a given window
#               Track the smallest and largest number of 1s in a given profile (column)
#               Track the avg, smallest, and largest number of 1s in a given window (row)

# how to track average number of windows present in a given NP?:
#   Average number of 1s in a given column, probably want a matrix



def extract_csv(file_name):
    data_indexed = pd.read_csv(file_name)
    if 'Unnamed' in data_indexed.columns[0]:
        data = data_indexed.iloc[:, 4:]
    else:
        data = data_indexed.iloc[:, 3:]
    return data_indexed, data

# extract the NPs of interest to us at the moment
# function is aggregate function to perform,
# result is conditional used to drop columns
# data_frame is the pandas data frame to operate on
def extract_data(function, result, data_frame):
    res = data_frame.agg(function)
    for i in res.index:
        if res[i] == result:
            res = res.drop(i)
    return res


def determine_radial_positions(data_frame):
    # sum up the columns
    temp = data_frame.agg('sum', axis=0)
    np_tuples = []
    for i in range(temp.size):
        np_tuples.append((temp[i], data_frame.columns[i]))

    # create list that tags each np with a label, marking its location in the nucleus
    radial_positions = []
    for i in np_tuples:
        if i[0] >= 326:
            radial_positions.append((i[0], 'strongly equatorial'))
        elif i[0] >= 244:
            radial_positions.append((i[0], 'somewhat equatorial'))
        elif i[0] >= 163:
            radial_positions.append((i[0], 'neither apical nor equatorial'))
        elif i[0] >= 81:
            radial_positions.append((i[0], 'somewhat apical'))
        else:
            radial_positions.append((i[0], 'strongly apical'))

    return radial_positions


def determine_compaction(data_frame, sorted_data):
    mean = data_frame.agg('sum', axis=1).values.mean()
    std_dev = data_frame.agg('sum', axis=1).values.std()

    compaction = []
    for i in sorted_data:
        if i[0] <= (mean-std_dev):
            compaction.append((i[0], i[1], 10))
        elif i[0] <= mean:
            compaction.append((i[0], i[1], 9))
        elif i[0] <= (mean+std_dev):
            compaction.append((i[0], i[1], 8))
        elif i[0] <= (mean+2*std_dev):
            compaction.append((i[0], i[1], 7))
        elif i[0] <= (mean+ 3*std_dev):
            compaction.append((i[0], i[1], 6))
        elif i[0] <= (mean+4*std_dev):
            compaction.append((i[0], i[1], 5))
        elif i[0] <= (mean+ 5*std_dev):
            compaction.append((i[0], i[1], 4))
        elif i[0] <= (mean+ 6*std_dev):
            compaction.append((i[0], i[1], 3))
        elif i[0] <= (mean+ 7*std_dev):
            compaction.append((i[0], i[1], 2))
        elif i[0] <= (mean+ 8*std_dev):
            compaction.append((i[0], i[1], 1))

    return compaction



def main():

    indexed, data = extract_csv('./data_set.csv')

    # data = pd.read_csv('./data_set.csv')
    # #drop tables that I don't want currently
    # data = data.drop(['chrom', ' start', ' stop'], axis=1)
    
    # print('Nuclear Profiles: ', len(data.columns))
    # print('Genomic Windows: ', len(data) - 1)

    # # how many windows are present in an np (ones in a column)
    # avgs = data.agg('sum', axis=0)
    # avg = 0
    # for i in avgs:
    #     avg += float(i)
    # avg = avg / len(data.columns)

    # print('Average windows per nuclear profiles: ', avg)
    # print('Least amount of windows in a nuclear profile: ', avgs.agg('min'))
    # print('Most amount of windows in a nuclear profile: ', avgs.agg('max'))

    # # In how many nuclear profiles is a window detected (ones in a row)
    # avgs_row = data.agg('sum', axis=1)
    # avg_row = 0
    # for i in avgs:
    #     avg_row += float(i)
    # avg_row = avg_row / len(data)

    # print('Average amount of window detections in a nuclear profile: ', avg_row)
    # print('Least amount of window detections in a nuclear profile: ', avgs_row.agg('min'))
    # print('Most amount of window detections in a nuclear profile: ', avgs_row.agg('max'))


    ################################################################################
    ################################################################################

    # Nuclear Profiles:  408
    # Genomic Windows:  90876
    # Average windows per nuclear profiles:  5482.811274509804
    # Least amount of windows in a nuclear profile:  31
    # Most amount of windows in a nuclear profile:  21249
    # Average amount of window detections in a nuclear profile:  24.615546287839607
    # Least amount of window detections in a nuclear profile:  0
    # Most amount of window detections in a nuclear profile:  408


    # Determine high detection frequency windows (outliers) - number of NPs in window
    # Estimate radial position of each NP on scale of 1 to 5
    # Estimate compaction of each window

    # Estimate radial position of each NP on scale of 1 to 5
    #  More often detected => closer to center
    #  Have windows sorted by detection rate
    #  Every 20% increase of detection indicates what level - at most 408
    #  >= 0: strongly apical
    #  >= 81.6: somewhat apical
    #  >= 163.2: neither apical nor equatorial
    #  >= 244.8: somewhat equatorial
    #  >= 326.4: strongly equatorial

    # sort pandas data
    # creates a sum of each row (90877 sums), sorts those values
    print(data)
    sorted_data = np.sort(data.agg('sum', axis=1).values)

    temp = data.agg('sum', axis=1)#.values
    sorted_data = []
    for i in range(temp.size):
        sorted_data.append((temp[i], i))
    sorted_data.sort()
    sorted_data = sorted_data[::-1]

    print(sorted_data)
    
    avg_sorted = 0
    for i in sorted_data:
        avg_sorted += i[0]
    avg_sorted = avg_sorted/len(sorted_data)

    print('Average: ', avg_sorted)

    # Radial positions
    radial_positions = determine_radial_positions(data)

    print(radial_positions)

    # Estimation of compaction: Just take avg_sorted from earlier and them ratings from 1-10
    
    compaction = determine_compaction(data, sorted_data)
    print(compaction)

    # for i in range(20):
    #     print(sorted_data[i])

    return


main()

