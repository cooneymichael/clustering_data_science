import cluster as cl
import pandas as pd
import numpy as np
from collections import Counter
import math
import matplotlib.pyplot as plt

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
    indexed, data = extract_csv('./CSVs/HIST1.csv')
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

    global_sim_list = []
    global_cluster_list = []
    global_center_list = []

    for outer_index in range(5):
        num_clusters = 3
        cluster = cl.Clusters(num_clusters, jaccard_dis)

        cluster.set_random_centers()
        cluster_list = cluster.assign_clusters()
        cluster.set_cluster_list(cluster_list)

        previous_iteration = cluster_list
        while True:
            cluster.set_new_centers()
            cluster_list = cluster.assign_clusters()
            cluster.set_cluster_list(cluster_list)

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
                # global_cluster_list.append(cluster_list)
                global_cluster_list.append(cluster)
                # global_center_list.append(new_centers)
                print("DONE ", outer_index)
                break
            else:
                previous_iteration = cluster_list

        # end while loop
    # end for loop
    # print(global_cluster_list)
    # find the best set of clusters, i.e. most similar
    for i in range(len(global_sim_list)):
        global_sim_list[i] = sum(global_sim_list[i])
    best_index = np.argmax(global_sim_list)
    cluster_list = global_cluster_list[best_index].get_cluster_list()
    print(len(cluster_list))

    # make some heat maps
    heat_map_lists = [[] for i in range(num_clusters)]
    for i in range(num_clusters):
        for j in range(len(cluster_list[i])):
            heat_map_lists[i].append(hist1_data_frame.iloc[:, cluster_list[i][j]])

    # create and config our plots
    fig, ((ax0, ax1, ax2),(ax3, ax4, ax5)) = plt.subplots(nrows=2, ncols=3)
    plt.style.use('bmh')
    plt.tight_layout()

    # heat maps
    im0 = ax0.imshow(heat_map_lists[0], cmap='binary_r', interpolation='none')
    cbar0 = ax0.figure.colorbar(im0, ax=ax0)
    ax0.set_aspect('auto')
    ax0.set_title('Cluster 1 Similarity')

    im1 = ax1.imshow(heat_map_lists[1], cmap='binary_r', interpolation='none')
    cbar1 = ax1.figure.colorbar(im1, ax=ax1)
    ax1.set_aspect('auto')
    ax1.set_title('Cluster 2 Similarity')

    im2 = ax2.imshow(heat_map_lists[2], cmap='binary_r', interpolation='none')
    cbar2 = ax2.figure.colorbar(im2, ax=ax2)
    ax2.set_aspect('auto')
    ax2.set_title('Cluster 3 Similarity')


    # Making radar chart
    # adding to the script: box plots and hist1 vs lads
    # chr13: 21,690,000-24,120,000
    feature_table = pd.read_csv('CSVs/Hist1_region_features.csv', sep=',')
    columns_of_interest = \
        ["Hist1", "Vmn", "LAD", "RNAPII-S2P", "RNAPII-S5P", "RNAPII-S7P", "Enhancer", "H3K9me3", "H3K20me3", "h3k27me3", "H3K36me3", "NANOG", "pou5f1", "sox2", "CTCF-7BWU"]
    feature_lists = [[] for i in columns_of_interest]
    indexed_hist1_data_frame = indexed.loc[:, columns_indexed]

    # cycle thru feature table and find windows that contain each column of interest 
    for i in range(len(columns_of_interest)):
        for j in feature_table.loc[:, [columns_of_interest[i], 'name']].values:
            if j[0] >= 1:
                feature_lists[i].append(j[1][6:14])

        # cycle through data frame and get corresponding windows, without all info
    interest_lists = [[] for i in columns_of_interest]
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
    percentage_frames = []
    for i in interest_lists:
        percentage_frames.append(indexed_hist1_data_frame[indexed_hist1_data_frame[' start'].isin(i)])
    
    # make nested lists (num_cluster lists * 15)
    cluster_percents = [ [[],[],[]] for i in columns_of_interest]

    for k in range(15):
        for i in range(len(cluster_list)):
            for j in range(len(cluster_list[i])):
                # add 4 to value because first 4 columns are various indices
                cluster_percents[k][i].append(
                    percentage_frames[k].iloc[:, cluster_list[i][j] + 4].astype(bool).sum() / len(percentage_frames[k].index)
                )


    mean_percents = [[] for i in range(num_clusters)]
    for i in range(len(cluster_percents)):
        mean_percents[0].append(np.mean(cluster_percents[i][0]))
        mean_percents[1].append(np.mean(cluster_percents[i][1]))
        mean_percents[2].append(np.mean(cluster_percents[i][2]))


    # plotting cluster 1 on radar chart
    mean_percents[0] += mean_percents[0][:1]
    angles = np.linspace(0,2 * np.pi, len(columns_of_interest), endpoint=False).tolist()
    angles += angles[:1]

    fix, ax6 = plt.subplots(figsize=(15,15), subplot_kw = dict(polar=True))

    # ax6 = fig.add_axes([0.6,0.1,0.4,0.4], polar=True)
    ax6.grid(color='#4c4c4c')
    ax6.set_theta_offset(np.pi / 2)
    ax6.set_theta_direction(-1)
    ax6.set_thetagrids(np.degrees(angles[0:len(angles)-1]), columns_of_interest)
    for label, angle in zip(ax6.get_xticklabels(), angles):
        if angle in (0,np.pi):
            label.set_horizontalalignment('center')
        elif 0 < angle < np.pi:
            label.set_horizontalalignment('left')
        else:
            label.set_horizontalalignment('right')

    ax6.plot(angles, mean_percents[0], color='red', linewidth=1, label='Cluster 1')
    ax6.fill(angles, mean_percents[0], color='red', alpha=0.25)

    # plotting cluster 2 radar
    mean_percents[1] += mean_percents[1][:1]
    ax6.plot(angles, mean_percents[1], color='blue', linewidth=1, label='Cluster 2')
    ax6.fill(angles, mean_percents[1], color='blue', alpha=0.25)

    # plotting cluster 3 radar
    mean_percents[2] += mean_percents[2][:1]
    ax6.plot(angles, mean_percents[2], color='yellow', linewidth=1, label='Cluster 3')
    ax6.fill(angles, mean_percents[2], color='yellow', alpha=0.25)

    # other elements of the chart
    ax6.set_title('Regions of Interest Across Clusters', y=1.08)
    ax6.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1))


    # plotting radial positions of each cluster
    entire_genome = pd.read_csv('CSVs/data_set.csv', sep=',')
    if 'Unnamed' in entire_genome.columns[0]:
        entire_genome = entire_genome.iloc[:, 4:]
    else:
        entire_genome = entire_genome.iloc[:, 3:]

    
    radial_positions = determine_radial_positions(entire_genome)
    radial_position_percents = [ [0,0,0,0,0] for i in cluster_list ]
    for i in range(len(cluster_list)):
        for j in range(len(cluster_list[i])):
            if radial_positions[cluster_list[i][j]][1] == 'strongly apical':
                radial_position_percents[i][0] += 1
            elif radial_positions[cluster_list[i][j]][1] == 'somewhat apical':
                radial_position_percents[i][1] += 1
            elif radial_positions[cluster_list[i][j]][1] == 'neither apical nor equatorial':
                radial_position_percents[i][2] += 1
            elif radial_positions[cluster_list[i][j]][1] == 'somewhat equatorial':
                radial_position_percents[i][3] += 1
            elif radial_positions[cluster_list[i][j]][1] == 'strongly equatorial':
                radial_position_percents[i][4] += 1
        for j in range(len(radial_position_percents[i])):
            radial_position_percents[i][j] /= len(cluster_list[i])

    x_axis = ['strongly\n apical', 'somewhat\n apical', 'neither', 'somewhat\n equatorial', 'strongly\n equatorial']

    for i in radial_position_percents:
        ax3.bar(x_axis, i, alpha=0.7, width=0.5)
        ax3.grid(color='#95a5a6', linestyle='--', linewidth=2, axis='y', alpha=0.7)
        ax3.set_title('Cluster 1')

        ax4.bar(x_axis, i, alpha=0.7, width=0.5)
        ax4.grid(color='#95a5a6', linestyle='--', linewidth=2, axis='y', alpha=0.7)
        ax4.set_title('Cluster 2')

        ax5.bar(x_axis, i, alpha=0.7, width=0.5)
        ax5.grid(color='#95a5a6', linestyle='--', linewidth=2, axis='y', alpha=0.7)
        ax5.set_title('Cluster 3')


    plt.show()



    
    

main()        
