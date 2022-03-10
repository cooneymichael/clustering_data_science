import pandas as pd
import numpy as np
import matplotlib as plt
import random
from quickselect import floyd_rivest as fr

class Clusters:
    def __init__(self, k, distance_matrix):
        # might not use iterations internally - consider removing

        self.num_clusters = k
        # self.iterations = iterations
        self.distance_matrix = distance_matrix
        self.centers = np.ndarray(shape=(self.num_clusters, 1), dtype=int)
        self.centers = np.append(self.centers, [[] for i in range(self.num_clusters)])
        self.previous_centers = []
        self.cluster_list = []
        self.global_cluster_list = []
        for i in range(self.num_clusters):
            self.global_cluster_list.append([])

    def get_cluster_list(self):
        """ 
        allow the user to see the cluster list
        """
        return self.cluster_list

    def set_cluster_list(self, new_list):
        # why am I allowing this?
        """ 
        allow the user to set the cluster list
        """
        if len(new_list) == self.num_clusters:
            self.cluster_list = new_list

    def assign_clusters(self):
        """
        take a list of centers and use it to determine where data points should be placed
        """
        local_cluster_list = []
        for i in self.centers:
            local_cluster_list.append([])

        for i in range(len(self.distance_matrix)):
            temp_list = [0] * len(local_cluster_list)
            for j in range(len(self.centers)):
                # if i == self.centers[j][0]:
                if i == self.centers[j]:
                    # guarantee that a center gets assigned to its cluster
                    temp_list[j] = -100
                else:
                    # temp_list[j] = self.distance_matrix[i][self.centers[j][0]]
                    temp_list[j] = self.distance_matrix[i][int(self.centers[j])]
            index = np.argmin(temp_list)
            local_cluster_list[index].append(i)
        return local_cluster_list

    def set_random_centers(self):
        """
        randomly choose centers, making sure there are no instances of two clusters having the same center at the same time
        """
        for i in range(len(self.centers)):
            center_i = random.randint(0, 162)
            #while center_i in self.centers[:, 0]:
            while center_i in self.centers:
                center_i = random.randint(0,162)
            # self.centers[i][0] = center_i
            self.centers[i] = center_i

    def set_new_centers(self):
        """
        find new centers using the current iteration of clusters and their dissimilarity
        """
        new_centers_wrt_clusters = np.ndarray(shape=(self.num_clusters,1), dtype=int)
        for i in range(len(new_centers_wrt_clusters)):
            new_centers_wrt_clusters[i][0] = 0

        # find the average dissimilarity within clusters
        distance_list = [[] for i in range(self.num_clusters)]
        for k in range(len(self.cluster_list)):
            for i in range(len(self.cluster_list[k])):
                dis_sum = 0
                for j in range(len(self.cluster_list[k])):
                    dis_sum += self.distance_matrix[i][j]
                distance_list[k].append(dis_sum)
                
        # find new centers and assign them to a matrix
        for i in range(len(distance_list)):
            median = fr.nth_largest(distance_list[i], len(distance_list[i])//2)
            index = [x for x,e in enumerate(distance_list[i]) if e == median][0]
            new_centers_wrt_clusters[i][0] = index

        for i in range(len(self.centers)):
            index = new_centers_wrt_clusters[i][0]
            self.centers[i] = self.cluster_list[i][index]
            


    
        


"""
Plan:
Make code object oriented:
- make a class
- what is my class?:
-   Clusters class
-   internal variables:
-     amount of clusters
-     amount of cluster iterations
-     lists for each clusters
-     cosine similarities for each clusters

-   internal functions/functionality
-     create clusters with random centers randomly
-     iterate thru and find better centers
-     terminate when best clusters found
-     choose best iteration of clusters
-     create plots
-     display/save plots



Add version control system
- make this a git repo
- reorganize files
- make files into packages for easy access/orgaization?



Update program
- at each iteration, don't allow random centers to be repeated from previous iterations
- change similarity measure to be intra-cluster similarity
- try to figure out why we need best 2 out of 3 clusters to terminate loop



Can we write this in R?
"""
