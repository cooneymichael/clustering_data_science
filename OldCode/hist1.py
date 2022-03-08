import pandas as pd
import numpy as np


def main():

    data_indexed = pd.read_csv('./HIST1.csv')
    # data = pd.read_csv('./HIST1.csv')
    #remove unwanted column that may or may not appear 
    if 'Unnamed' in data_indexed.columns[0]:
        data = data_indexed.iloc[:, 4:]
    else:
        data = data_indexed.iloc[:, 3:]

    # extract NPs which contain at least 1 window in region of interest
    sums = data.agg('sum')
    for i in sums.index:
        if sums[i] == 0:
            sums = sums.drop(i)

    # Basic statistics
    print('BASIC STATISTICS')
    print('NPs with at least 1 window in HIST1: ', sums.index)

    print('Average number of windows present: ', sums.mean())

    print('Least number of windows present in an NP: ', sums.agg('min'))
    print('Most number of windows present in an NP: ', sums.agg('max'))

    sums_of_rows = data.agg('sum', axis=1)
    print('Average number of NPs in which a window is detected: ', sums_of_rows.mean())
    print('Least number of NPs in which a window is detected: ', sums_of_rows.min())
    print('Most number of NPs in which a window is detected: ', sums_of_rows.max())


    # What are the most common radial positions of the HIST1 region based on the
    # NPs that captured the area?

    # do I need to use all data, or can I just perform on reduced data?
    all_data = pd.read_csv('./data_set.csv')
    #remove unwanted column that may or may not appear 
    if 'Unnamed' in all_data.columns[0]:
        all_data = all_data.iloc[:, 4:]
    else:
        all_data = all_data.iloc[:, 3:]


    temp = all_data.agg('sum', axis=0)
    np_tuples = []
    for i in range(temp.size):
        np_tuples.append((temp[i], all_data.columns[i]))

    # create list that tags each np with a label, marking its location in the nucleus
    radial_positions = []
    for i in np_tuples:
        if i[0] >= 326:
            radial_positions.append((i[0], i[1], 'strongly equatorial'))
        elif i[0] >= 244:
            radial_positions.append((i[0], i[1], 'somewhat equatorial'))
        elif i[0] >= 163:
            radial_positions.append((i[0], i[1], 'neither apical nor equatorial'))
        elif i[0] >= 81:
            radial_positions.append((i[0], i[1], 'somewhat apical'))
        else:
            radial_positions.append((i[0], i[1], 'strongly apical'))

    hist_profiles_categorized = []
    for i in radial_positions:
        if i[1] in sums.index:
            hist_profiles_categorized.append(i)

    strong_equ = 0
    some_equ = 0
    neither = 0
    some_ap = 0
    strong_ap = 0
    for i in hist_profiles_categorized:
        if i[2] == 'strongly equatorial':
            strong_equ += 1
        elif i[2] == 'somewhat equatorial':
            some_equ += 1
        elif i[2] == 'neither apical nor equatorial':
            neither += 1
        elif i[2] == 'somewhat apical':
            some_ap += 1
        elif i[2] == 'strongly apical':
            strong_ap += 1
        else:
            print('yikes')
    print()
    print('RADIAL POSITION')
    print('strongly equatorial: ', strong_equ)
    print('somewhat equatorial: ', some_equ)
    print('neither: ', neither)
    print('somewhat apical: ', some_ap)
    print('strongly apical: ', strong_ap)

    # What are the typical compactions of the windows within the HIST1 region?


    temp = all_data.agg('sum', axis=1)#.values
    sorted_all_data = []
    for i in range(temp.size):
        sorted_all_data.append((temp[i], i))
    sorted_all_data.sort()
    sorted_all_data = sorted_all_data[::-1]

    # print(sorted_all_data)
    

    # Estimation of compaction: Just take sorted_all_data from earlier and add ratings from 1-10
    mean = all_data.agg('sum', axis=1).values.mean()
    std_dev = all_data.agg('sum', axis=1).values.std()

    compaction = []
    # don't think this is working but I'm too tired right now
    for i in sorted_all_data:
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


    compaction_hist1 = []
    for i in compaction:
        if i[1] in data_indexed.loc[:, 'Unnamed: 0'].values:
            compaction_hist1.append(i)


    average_compaction = 0
    for i in compaction_hist1:
        average_compaction += i[2]

    average_compaction /= len(compaction_hist1)
    print()
    print("COMPACTION:")
    print('Average compaction: ', average_compaction)
    
    return

main()
