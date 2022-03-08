import pandas as pd
import numpy as np
import matplotlib.pyplot as matpl
import random
from quickselect import floyd_rivest as fr
from collections import Counter
import math

# read the csv file and remove undesired columns
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

# input is a pandas data frame
# modified: denominator of jaccard index is:
# min{A, B}, where A = m_11+m_10 for column A and B is same for column B
def jaccard_similarity(data):
    size = data.columns.size
    jaccard_matrix = np.ndarray(shape=(size, size), dtype=float)

    index_1 = 0

    for _, column_A in data.iteritems():
        index_2 = 0
        for __, column_B in data.iteritems():
            A = 0
            B = 0
            m_11 = 0 

            for i in range(column_A.values.size):

                if column_A.values[i] == 1 and column_B.values[i] == 1:
                    m_11 += 1
                    A += 1
                    B += 1

                if column_A.values[i] == 0 and column_B.values[i] == 1:
                    # m_01 += 1
                    B += 1

                if column_A.values[i] == 1 and column_B.values[i] == 0:
                    A += 1
                    # m_10 += 1

            index = m_11 / min([A, B])
            jaccard_matrix[index_1][index_2] = index
            index_2 += 1
            
        index_1 += 1
    return jaccard_matrix


def jaccard_distance(j_matrix):
    # go through jaccard similarity and subtract each value from 1.
    # save that value in the new matrix.  return the matrix
    size = j_matrix.shape[0]
    jaccard_dis = np.ndarray(shape=(size, size), dtype=float)

    for i in range(size):
        for j in range(size):
            jaccard_dis[i][j] = 1 - j_matrix[i][j]

    return jaccard_dis


def assign_clusters(center_list, dists):
    # print("center list assignment: ", center_list)

    cluster_list = []
    for i in center_list:
        cluster_list.append([])

    for i in range(len(dists)):
        temp_list = [0,0,0]
        for j in range(len(center_list)):
            if i == center_list[j][0]:
                temp_list[j] = -100
                # print("guaranteeing center: ", i, " at ", j)
            else:
                temp_list[j] = dists[i][center_list[j][0]]
        index = np.argmin(temp_list)
        cluster_list[index].append(i)
    return cluster_list


def counter_cosine_similarity(c1, c2):
    terms = set(c1).union(c2)
    dotprod = sum(c1.get(k, 0) * c2.get(k, 0) for k in terms)
    magA = math.sqrt(sum(c1.get(k, 0)**2 for k in terms))
    magB = math.sqrt(sum(c2.get(k, 0)**2 for k in terms))
    return dotprod / (magA * magB)


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


def main():
    indexed, data = extract_csv('./HIST1.csv')
    hist1_data = extract_data('sum', 0, data)
    hist1_data_indexed = extract_data('sum', 0, indexed)

    columns = hist1_data.index.intersection(data.columns)
    columns_indexed = hist1_data_indexed.index.intersection(indexed.columns)
    hist1_data_frame = data.loc[:, columns]

    # jaccard index
    jaccard_sim = jaccard_similarity(hist1_data_frame)
    # print(jaccard_sim)

    jaccard_dis = jaccard_distance(jaccard_sim)
    # print(jaccard_dis)

    # get random numbers as centers for the first iteration
    num_clusters = 3

    # Make a global list of all intra-variance stuff (global_sim_list) ->
    #     whenever sentinel condition is true: add sim_list_2 to this
    # make a global list of all clusters -> when sentinel is true: add clusters to this
    # make a global list of all centers -> when sentinel is true: add centers to this
    # Set a for loop
    #    Do the stuff
    # Go through the global_sim_list, sum each set, choose index with highest sum (argmax)
    # make a heatmap
    global_sim_list = []
    global_cluster_list = []
    global_center_list = []

    for outer_index in range(5):
        
        rand_nums = np.ndarray(shape=(num_clusters,1), dtype=int)
        for i in range(len(rand_nums)):
            rand_nums[i][0] = 0
        for i in range(len(rand_nums)):
            center_i = random.randint(0, 162)
            while center_i in rand_nums[:,0]:
                center_i = random.randint(0, 162)
            rand_nums[i][0] = center_i
            
        # start making the clusters
        cluster_list = assign_clusters(rand_nums, jaccard_dis)
        # print(cluster_list)
        # print("========================================================================")
        
        oldest_iteration = cluster_list
        previous_iteration = cluster_list
        
        while True:
            # finding the new centers to get better clusters
            new_centers_wrt_clusters = np.ndarray(shape=(num_clusters,1), dtype=int)
            for i in range(len(new_centers_wrt_clusters)):
                new_centers_wrt_clusters[i][0] = 0
                
            # finding average dissimilarity within clusters
            distance_list = []
            for i in cluster_list:
                distance_list.append([])
            for k in range(len(cluster_list)):
                for i in range(len(cluster_list[k])):
                    dis_sum = 0
                    for j in range(len(cluster_list[k])):
                        dis_sum += jaccard_dis[i][j]
                    distance_list[k].append(dis_sum)
                                            
            # finding the new centers and assigning them to a matrix
            for i in range(len(distance_list)):
                median = fr.nth_largest(distance_list[i], len(distance_list[i])//2)
                index = [x for x,e in enumerate(distance_list[i]) if e == median][0]
                new_centers_wrt_clusters[i][0] = index
            
            # print("new centers wrt clusters: ", new_centers_wrt_clusters)
        
            new_centers = np.ndarray(shape=(num_clusters, 1), dtype=int)
            for i in range(len(new_centers)):
                new_centers[i][0] = 0
                
            for i in range(len(new_centers)):
                index = new_centers_wrt_clusters[i][0]
                new_centers[i][0] = cluster_list[i][index]
                
            # print("new centers wrt population: ", new_centers)
        
            cluster_list = assign_clusters(new_centers, jaccard_dis)
            # print(cluster_list)
            # print("===================================================================")
            
            # sentinel condition
            sim_list = []
            sim_list_2 = []
            for i in cluster_list:
                sim_list.append([])
                # sim_list_2.append([])
                
            for i in range(len(cluster_list)):
                sim_list[i].append(Counter(cluster_list[i]))
                sim_list[i].append(Counter(previous_iteration[i]))
                
            for i in range(len(sim_list)):
                sim_list_2.append(counter_cosine_similarity(sim_list[i][0], sim_list[i][1]))
                
            sentinel_list = []
            for i in cluster_list:
                sentinel_list.append(True)
                
            for i in range(len(sim_list_2)):
                if sim_list_2[i] <= 0.7:
                    # print("the bastard: ", i, " ", sim_list_2[i])
                    sentinel_list[i] = False
                
            if Counter(sentinel_list)[True] >= 2:
                global_sim_list.append(sim_list_2)
                global_cluster_list.append(cluster_list)
                global_center_list.append(new_centers)
                print("DONE")
                break
            else:
                previous_iteration = cluster_list
                
                
            # if cluster_list == previous_iteration or cluster_list == oldest_iteration:
            #     print("DONE")
            #     break
            # else:
            #     oldest_iteration = previous_iteration
            #     previous_iteration = cluster_list
            
        # return
        # end inner loop
    # outer loop
    # print("WE HAVE FINISHED")

    for i in range(len(global_sim_list)):
        global_sim_list[i] = sum(global_sim_list[i])
    # find argmax, use that global_cluster_list index
    best_index = np.argmax(global_sim_list)

    print(global_cluster_list[best_index])

    for i in range(len(global_cluster_list[best_index])):
        local_list = []
        for j in range(len(global_cluster_list[best_index][i])):
            local_list.append(hist1_data_frame.iloc[:, global_cluster_list[best_index][i][j]])
        print(local_list)
        matpl.imshow(local_list, cmap='hot', interpolation='none')
        matpl.colorbar(matpl.pcolor(local_list))
        matpl.title("Cluster map")
        matpl.show()


    # adding to the script: box plots and hist1 vs lads
    # chr13: 21,690,000-24,120,000
    feature_table = pd.read_csv('./Hist1_region_features.csv', sep=',')

    # calculate the percentage of NPs in each cluster that contain HIST1 genes
    # for each cluster:
    #   determine which window the NP is in
    #      need the start value from indexed data
    #      but each NP contains all windows, so what am I doing?
    #   reference that window in feature table
    #   add value of "Hist1" column to some data frame
    #   sum that data_frame and divide by the length of the cluster

    # calculate the percentage of NPs in each cluster that contain LAD genes
    # for each cluster:
    #   determine which window the NP is in
    #   reference that window in the feature table
    #   add value of "LAD" column to some data frame
    #   sum that data_frame and divide by the length of the cluster


    # lists of data frames, sorry for bad names
    # hist1_feature_table_list = []
    # lad_feature_table_list = []

    # create nested list of lists
    # will follow order of assignment:
    # Hist1, Vmn, LAD, RNAPII-S2P, RNAPII-S5P, RNAPII-S7P, Enhancer, H3K9me3, H3K20me3, h3k27me3, H3K36me3, NANOG, pou5f1, sox2, CTCF-7BWU
    columns_of_interest = \
        ["Hist1", "Vmn", "LAD", "RNAPII-S2P", "RNAPII-S5P", "RNAPII-S7P", "Enhancer", "H3K9me3", "H3K20me3", "h3k27me3", "H3K36me3", "NANOG", "pou5f1", "sox2", "CTCF-7BWU"]

    feature_lists = []
    for i in range(15):
        feature_lists.append([])
    

    indexed_hist1_data_frame = indexed.loc[:, columns_indexed]
    
    
    # cycle thru feature table and find windows with hist1
    # for i in feature_table.loc[:, ['Hist1', 'name']].values:
    #     if i[0] >= 1:
    #         hist1_feature_table_list.append(i[1][6:14])

    # for i in feature_table.loc[:, ['LAD','name']].values:
    #     if i[0] >= 1:
    #         lad_feature_table_list.append(i[1][6:14])


    # cycle thru feature table and find windows that contain each column of interest 
    for i in range(len(columns_of_interest)):
        for j in feature_table.loc[:, [columns_of_interest[i], 'name']].values:
            if j[0] >= 1:
                feature_lists[i].append(j[1][6:14])
                
    # cycle through data frame and get corresponding windows, without all info
    # hist1_list = []
    # lad_list = []
    # for i in indexed_hist1_data_frame.loc[:, ' start'].values:
    #     if str(i) in hist1_feature_table_list:
    #         hist1_list.append(i)
    #     elif str(i) in lad_feature_table_list:
    #         lad_list.append(i)

    # cycle through data frame and get corresponding windows, without all info
    interest_lists = []
    for i in range(15):
        interest_lists.append([])

    for i in indexed_hist1_data_frame.loc[:, ' start'].values:
        if str(i) in feature_lists[0]:
            interest_lists[0].append(i)

        if str(i) in feature_lists[1]:
            interest_lists[1].append(i)

        if str(i) in feature_lists[2]:
            interest_lists[2].append(i)

        if str(i) in feature_lists[3]:
            interest_lists[3].append(i)

        if str(i) in feature_lists[4]:
            interest_lists[4].append(i)

        if str(i) in feature_lists[5]:
            interest_lists[5].append(i)

        if str(i) in feature_lists[6]:
            interest_lists[6].append(i)

        if str(i) in feature_lists[7]:
            interest_lists[7].append(i)

        if str(i) in feature_lists[8]:
            interest_lists[8].append(i)

        if str(i) in set(feature_lists[9]):
            interest_lists[9].append(i)

        if str(i) in feature_lists[10]:
            interest_lists[10].append(i)

        if str(i) in feature_lists[11]:
            interest_lists[11].append(i)

        if str(i) in feature_lists[12]:
            interest_lists[12].append(i)

        if str(i) in feature_lists[13]:
            interest_lists[13].append(i)

        if str(i) in feature_lists[14]:
            interest_lists[14].append(i)

    
    # get all dataframe information related to necessary windows
    # hist1_percentage_frame = indexed_hist1_data_frame[indexed_hist1_data_frame[' start'].isin(hist1_list)]
    # lad_percentage_frame = indexed_hist1_data_frame[indexed_hist1_data_frame[' start'].isin(lad_list)]
    percentage_frames = []
    for i in interest_lists:
        percentage_frames.append(indexed_hist1_data_frame[indexed_hist1_data_frame[' start'].isin(i)])
        

    # make nested list (num_cluster lists)
    # hist1_cluster_percents = []
    # lad_cluster_percents = []
    # for i in global_cluster_list[best_index]:
    #     hist1_cluster_percents.append([])
    #     lad_cluster_percents.append([])

    # for i in range(len(global_cluster_list[best_index])):
    #   for j in range(len(global_cluster_list[best_index][i])):
    #       # add 4 to value because first 4 columns are various indices
    #       hist1_cluster_percents[i].append(
    #           hist1_percentage_frame.iloc[:, global_cluster_list[best_index][i][j] + 4].agg('sum') / len(hist1_percentage_frame.index)
    #       )

    #       lad_cluster_percents[i].append(
    #           lad_percentage_frame.iloc[:, global_cluster_list[best_index][i][j] + 4].agg('sum') / len(lad_percentage_frame.index)
    #       )

    # make nested lists (num_cluster lists * 15)
    cluster_percents = []
    for i in range(15):
        cluster_percents.append([[],[],[]])


    # print()
    # print('CLUSTER PERCENTS: ', cluster_percents)


    for k in range(15):
        for i in range(len(global_cluster_list[best_index])):
            for j in range(len(global_cluster_list[best_index][i])):
                # add 4 to value because first 4 columns are various indices
                cluster_percents[k][i].append(
                    percentage_frames[k].iloc[:, global_cluster_list[best_index][i][j] + 4].astype(bool).sum() / len(percentage_frames[k].index)
                )

    # for i in range(len(cluster_percents)):
    #     print(i)
    #     # for j in range(len(cluster_percents[i])):
    #     #     print(j)
    #     print(cluster_percents[i])
    #     print()

    # print('ORIGINAL:')
    # print(hist1_cluster_percents)
    # print('NEW:')
    # # for i in cluster_percents[0]:
    # #     print(i)
    # print(cluster_percents[0])
    # print()

    mean_percents = [[],[],[]]
    for i in range(len(cluster_percents)):
        mean_percents[0].append(np.mean(cluster_percents[i][0]))
        mean_percents[1].append(np.mean(cluster_percents[i][1]))
        mean_percents[2].append(np.mean(cluster_percents[i][2]))

    print(mean_percents)
    # for i in range(len(percentage_frames)):
    #     print(i)
    #     print(percentage_frames[i])
    #     print()


    # make a radar chart using the 15 sub lists, 1 chart for each cluster
    # columns_of_interest is our labels
    # values are the first list of each nested list
    
    # plotting cluster 1 radar
    mean_percents[0] += mean_percents[0][:1]
    angles = np.linspace(0,2 * np.pi, len(columns_of_interest), endpoint=False).tolist()
    angles += angles[:1]

    fix, ax = matpl.subplots(figsize=(15,15), subplot_kw = dict(polar=True))

    ax.grid(color='#4c4c4c')
    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    ax.set_thetagrids(np.degrees(angles[0:len(angles)-1]), columns_of_interest)
    for label, angle in zip(ax.get_xticklabels(), angles):
        if angle in (0,np.pi):
            label.set_horizontalalignment('center')
        elif 0 < angle < np.pi:
            label.set_horizontalalignment('left')
        else:
            label.set_horizontalalignment('right')

    ax.plot(angles, mean_percents[0], color='red', linewidth=1, label='cluster 1')
    ax.fill(angles, mean_percents[0], color='red', alpha=0.25)

    # plotting cluster 2 radar
    mean_percents[1] += mean_percents[1][:1]
    ax.plot(angles, mean_percents[1], color='blue', linewidth=1, label='cluster 2')
    ax.fill(angles, mean_percents[1], color='blue', alpha=0.25)

    # plotting cluster 3 radar
    mean_percents[2] += mean_percents[2][:1]
    ax.plot(angles, mean_percents[2], color='yellow', linewidth=1, label='cluster 3')
    ax.fill(angles, mean_percents[2], color='yellow', alpha=0.25)

    # other elements of the chart
    ax.set_title('Regions of Interest Across Clusters', y=1.08)
    ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1))

    matpl.show()


    # print('================================================================================')
    # print('PERCENTS:')
    # print(hist1_cluster_percents)

    # show results in boxplots
    # matpl.boxplot(hist1_cluster_percents)
    # matpl.show()

    # matpl.boxplot(lad_cluster_percents)
    # matpl.show()

    entire_genome = pd.read_csv('data_set.csv', sep=',')
    if 'Unnamed' in entire_genome.columns[0]:
        entire_genome = entire_genome.iloc[:, 4:]
    else:
        entire_genome = entire_genome.iloc[:, 3:]

    
    radial_positions = determine_radial_positions(entire_genome)

    radial_position_percents = []
    for i in global_cluster_list[best_index]:
        radial_position_percents.append([0,0,0,0,0])

    for i in range(len(global_cluster_list[best_index])):
        for j in range(len(global_cluster_list[best_index][i])):
            if radial_positions[global_cluster_list[best_index][i][j]][1] == 'strongly apical':
                radial_position_percents[i][0] += 1
            elif radial_positions[global_cluster_list[best_index][i][j]][1] == 'somewhat apical':
                radial_position_percents[i][1] += 1
            elif radial_positions[global_cluster_list[best_index][i][j]][1] == 'neither apical nor equatorial':
                radial_position_percents[i][2] += 1
            elif radial_positions[global_cluster_list[best_index][i][j]][1] == 'somewhat equatorial':
                radial_position_percents[i][3] += 1
            elif radial_positions[global_cluster_list[best_index][i][j]][1] == 'strongly equatorial':
                radial_position_percents[i][4] += 1
        for j in range(len(radial_position_percents[i])):
            radial_position_percents[i][j] /= len(global_cluster_list[best_index][i])

    # print(radial_position_percents)
    x=0
    for i in radial_position_percents:
        x += 1
        x_axis = ['strongly apical', 'somewhat apical', 'neither', 'somewhat equatorial', 'strongly equatorial']
        matpl.bar(x_axis, i, alpha=0.7, width=0.5)
        matpl.grid(color='#95a5a6', linestyle='--', linewidth=2, axis='y', alpha=0.7)
        matpl.title('Cluster ' + str(x))
        matpl.show()
    

    return
    
main()


""" 
Plan for radar chart
use nested lists to represent the percent lists of each feature, i.e. list[0] = hist1, list[2] = lad, etc
same for percentage frames (nested list of frames)
when summing, we have to add 1 if > 0, not just sum
"""
