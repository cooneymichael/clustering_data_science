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

def main():

    data = pd.read_csv('./data_set.csv')
    #drop tables that I don't want currently
    data = data.drop(['chrom', ' start', ' stop'], axis=1)
    
    print('Nuclear Profiles: ', len(data.columns))
    print('Genomic Windows: ', len(data) - 1)

    # how many windows are present in an np (ones in a column)
    avgs = data.agg('sum', axis=0)
    avg = 0
    for i in avgs:
        avg += float(i)
    avg = avg / len(data.columns)

    print('Average windows per nuclear profiles: ', avg)
    print('Least amount of windows in a nuclear profile: ', avgs.agg('min'))
    print('Most amount of windows in a nuclear profile: ', avgs.agg('max'))

    # In how many nuclear profiles is a window detected (ones in a row)
    avgs_row = data.agg('sum', axis=1)
    avg_row = 0
    for i in avgs:
        avg_row += float(i)
    avg_row = avg_row / len(data)

    print('Average amount of window detections in a nuclear profile: ', avg_row)
    print('Least amount of window detections in a nuclear profile: ', avgs_row.agg('min'))
    print('Most amount of window detections in a nuclear profile: ', avgs_row.agg('max'))

    return


main()

