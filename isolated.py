import os
import networkx as nx
import csv
import pandas as pd

# just a helper function to return i + 1 if at index i in the array it has the maximum value
def largestD(arr):
    temp = 0
    for i in range(1,len(arr)):
        if arr[i] > arr[temp]:
            temp = i
    return temp + 1

# return isolated points given an assignment DataFrame
def find_isolated(assign_DF):
    isolated = []
    
    for i in range(1,9):
        district = 'd' + str(i)
        G = nx.Graph()

        node_list = list(assign_DF[assign_DF['dc']==district]['i'])

        G.add_nodes_from(node_list)

        for i in node_list:
            for j in adj_dict[i]:
                if j in node_list:
                    G.add_edge(i,j)

        for i in node_list:
            if not(nx.has_path(G,root[district],i)):
                isolated.append(i)
            
    return isolated

####################################
#Start of the main programming flow#
####################################
# convert assignment gdx to csv to make it easier to communicate with python
os.system("gdxdump initial_assignment.gdx format=csv output=initial_assignment.csv symb=assign")

##########################################################
#create an empty dictionary to store each node's neighbor#
##########################################################
adj_dict = {}

for i in range(1,7079):
    node_name = 'F'+str(i)
    adj_dict[node_name] = []
    
with open('border_new.csv', 'r') as csvfile:
    adj = csv.reader(csvfile)
    temp = 0
    for i in adj:
        adj_dict[i[0]].append(i[1])
        adj_dict[i[1]].append(i[0])
        temp = temp + 1
        
for i in range(1,7079):
    node_name = 'F'+str(i)
    adj_dict[node_name] = set(adj_dict[node_name])


#################################
#Read in the assignment csv file#
#################################
assign_df = pd.read_csv('initial_assignment.csv')
assign_df = assign_df[assign_df['Val']==1]


#################################
#Specify Center (to be modified)#
#################################
root = {}
root['d1'] = 'F1'
root['d2'] = 'F1000'
root['d3'] = 'F2000'
root['d4'] = 'F3000'
root['d5'] = 'F4000'
root['d6'] = 'F5000'
root['d7'] = 'F6000'
root['d8'] = 'F7000'

##########################
#Reassign isolated points#
##########################

isolated = find_isolated(assign_df)

# keep reassigning isolated nodes until no node is isolated
while(len(isolated) != 0):
    for i in isolated: # parse all isolated wards
        assign_vote = [0,0,0,0,0,0,0,0] # keep number of adjacent unisolated wards for each district
        index = assign_df[assign_df['i'] == i].index[0]
        for j in adj_dict[i]:
            if j not in isolated:
                assigned = str(assign_df[assign_df['i'] == j].iloc[0,1][-1:])       
                assign_vote[int(assigned)-1] += 1
        assign_df.loc[index,'dc'] = 'd' + str(largestD(assign_vote))
    isolated = find_isolated(assign_df)


assign_df.to_csv('isolated_reassign.csv',index=False)

# convert assignment csv to gdx
os.system("csv2gdx isolated_reassign.csv output=isolated_reassign.gdx id=isolated_reassign index=1,2 useHeader=y value=3")



