#!/usr/bin/env python
# coding: utf-8

# In[ ]:

#This code requires you to be inside of an All of Us Jupyter Notebook. 
#This version of the code is best for obtaining the CSLVs of whole chromosomes.
#For targeted CSLVs, look for another script in my Github.

#Initializing code

from datetime import datetime
import glob
import os
import pandas as pd
import hail as hl 
import numpy as np
import re
start = datetime.now()


# In[ ]:


#bucket associated to the workbench address
bucket = os.getenv('WORKSPACE_BUCKET')
bucket


# In[ ]:


genomic_location = os.getenv("CDR_STORAGE_PATH")
genomic_location


# In[ ]:


#Importing hail and initializing Spark
hl.init(default_reference="GRCh38", idempotent=True)


# In[ ]:


# we use env varibale MICRO.... to load microarray matrix table
# Directory of individual paths for genetic data, https://support.researchallofus.org/hc/en-us/articles/4616869437204-Controlled-CDR-Directory
mt_array_path = os.getenv("MICROARRAY_HAIL_STORAGE_PATH")
mt_array_path


# In[ ]:


#read hail matrix table and store it mt_array
mt= hl.read_matrix_table(mt_array_path)
# to count the number of rows and columns in MatrixTable.
mt.count()


# In[ ]:


mt.show(n_rows=10,n_cols=5)


# In[ ]:


# drop entry fields that we dont need to process, just keep LRR for each patient which is a colum with person id

mt= mt.drop("BAF","GT","IGC","NORMX","NORMY","R","THETA","X","Y")


# In[ ]:

#The new method I'm using for CSLV Generation requires very precise base pair selection.
#I used the UCSC Genome Browser in order to find the last base pair for each chromosome.

chr1_end = 248956422
chr2_end = 242193529
chr3_end = 198295559
chr4_end = 190214555
chr5_end = 181538259
chr6_end = 170805979
chr7_end = 159345973
chr8_end = 145138636
chr9_end = 138394717
chr10_end = 133797422
chr11_end = 135086622
chr12_end = 133275309
chr13_end = 114364328
chr14_end = 107043718
chr15_end = 101991189
chr16_end = 90338345
chr17_end = 83257441
chr18_end = 80373285
chr19_end = 58617616
chr20_end = 64444167
chr21_end = 46709983
chr22_end = 50818468

chr_end = [
    248956422,
    242193529,
    198295559,
    190214555,
    181538259,
    170805979,
    159345973,
    145138636,
    138394717,
    133797422,
    135086622,
    133275309,
    114364328,
    107043718,
    101991189,
    90338345,
    83257441,
    80373285,
    58617616,
    64444167,
    46709983,
    50818468
]


# In[ ]:

#This is the main part of the CSLV Generation Code. 
#Keep in mind, this code is best for extracting CSLVs from the entire chromosome.
#The general process is as follows:
#1. Creating intervals based on the chromosomes you're interested in and the size of the CSLV
#2. Creating CSLVs based on the intervals created previously
#3. Saving those CSLVs in your Workspace Bucket to be combined later with another script
#It also comes with a time tracker to see how long each CSLV takes to run, along with the total run time.

#Input Array (The current chromosomes are there for an example)
chrArray = ['chr16','chr17', 'chr18', 'chr19', 'chr20', 'chr21','chr22']
cslvSize = 10000000
intervals = []

#Use Regex to extract just the number value of the chromosome
#The compilation step makes it easier and quicker to find the number values

for chr in chrArray:
    regex = re.compile("chr([0-9]{1,2})")
    match = regex.search(chr)
    chrIndex = int(match.group(1))

    #Use the index we found to get max value allowed for said chromosome
    maxVal = chr_end[chrIndex - 1]
    lastIndex = maxVal // cslvSize + 1
    
    intervals = []
    
    # Generate Intervals to loop over
    for i in range(0, lastIndex, 1):
        intervals.append("chr{chrIndex}:{start}-{end}".format(chrIndex=chrIndex, start=i*cslvSize+1, end = ({True: maxVal, False: (i+1)*cslvSize}[i + 1 == lastIndex])))

    # Use the Intervals generated for whatver you want
    print(intervals)
    
    # Create a counter for the next step
    counter = 1
    
    # Attempt to start the collection of CSLVs given the selected chromosome and CSLV size
    for v in intervals:
        print("Beginning of run #" + str(counter) + ' out of ' + str(lastIndex))
        print(v)
        mt1 = hl.filter_intervals(
            mt,
            [hl.parse_locus_interval(v)] , keep=True)
        
        # Make sure everything is working
        print(mt1.count())
        mt1.show(n_rows=10,n_cols=5)
        
        # Create CSLVs by aggregating LRRs
        mt1_with_avg_LRR = mt1.annotate_cols(avg_of_LRR_by_chr=hl.agg.mean(mt1.LRR))
        result = mt1_with_avg_LRR.cols()
        df= result.to_pandas()
        df.rename(columns = {'s': 'person_id'}, inplace = True)
        
        # Save CSLV for later
        df.to_csv(f'{bucket}/data/Average_LRR_{chr}_CSLVSize_{str(cslvSize)}_Split_{str(counter)}.csv', index = False)
        
        # Keep track of how much time this particular section takes
        stop = datetime.now()
        total_time = str(stop - start)
        print("Time of finishing run #" + str(counter) + ":")
        print(total_time)
        print("End of run #" + str(counter) + ' out of ' + str(lastIndex))
        
        # Increase counter for next CSV
        counter += 1


# In[ ]:


import os
import subprocess
import numpy as np
import pandas as pd




# In[ ]:
#The code from here on out is not necessary for CSLV creation.
#I like to make sure that my CSLVs saved properly after the code is finished.

# This snippet assumes that you run setup first
# This code lists objects in your Google Bucket

# Get the bucket name
my_bucket = os.getenv('WORKSPACE_BUCKET')

# List objects in the bucket
print(subprocess.check_output(f"gsutil ls -r {my_bucket}", shell=True).decode('utf-8'))




# In[ ]:


# This snippet assumes you run setup first

# This code copies file in your Google Bucket and loads it into a dataframe

# Replace 'test.csv' with THE NAME of the file you're going to download from the bucket (don't delete the quotation marks)
# Average_LRR_chr22_CSLVSize_25000000_Split_1.csv is the average template
# Of course, this template will change depending on the chromosome used, CSLV size, and how many splits it takes for it to handle the entire chromosome

name_of_file_in_bucket = 'test.csv'

########################################################################
##
################# DON'T CHANGE FROM HERE ###############################
##
########################################################################

# get the bucket name
my_bucket = os.getenv('WORKSPACE_BUCKET')

# copy csv file from the bucket to the current working space
os.system(f"gsutil cp '{my_bucket}/data/{name_of_file_in_bucket}' .")

print(f'[INFO] {name_of_file_in_bucket} is successfully downloaded into your working space')
# save dataframe in a csv file in the same workspace as the notebook
my_dataframe = pd.read_csv(name_of_file_in_bucket)
print(my_dataframe.shape[0])
print(my_dataframe.shape[1])
my_dataframe.head()