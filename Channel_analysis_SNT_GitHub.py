


# Imports 

import csv

import pandas as pd
import seaborn as sns
import numpy as np
import glob 


#Data Path
path = "add the path of acquired data from SNT"

# Two separate files are generated and stored:
#one for pixel intensity calculation and the second file contains data for AIS length, Nav1.6 length etc.

path_calculation='***/Intensity_calculation_Nav1.6_Total_TEST.csv'
Aggregated_data='***/AggegratedData_Nav1.6_Total_TEST.csv'


celltype =''
firsttime = True

    
 # reads all csv files in the same folder
      
for file_name in glob.glob(path +'*.csv'): 
    print('test')
    x = np.genfromtxt(file_name,delimiter=',')[:,2] 
    df = pd.read_csv(file_name, encoding='utf-8', engine='python')
  
 #searches for filenames containing specific names   
    if file_name.find("_AcD")==-1:
        celltype ="nonAcD"
    else:
        celltype="AcD"
   
    
 # Parameters for channel analysis 

    #threshold for local average intensity 
    threshold = 40
    #number of pixels to consider in the local average
    numberpixels = 5
    
 # Variables   
  
    average_valuesCh1 = df.mean()[1]
    average_valuesCh2 = df.mean()[2]
    average_valuesCh3 = df.mean()[3]
    length = len(df.index)
    startindex_Ch1 = 0
    endindex_Ch1 = length-1
    startindex_Ch2 = 0
    endindex_Ch2 = length-1
    # Normalizes the intensity data of each channel  
    normalized_df=(df-df.min())/(df.max()-df.min())
    
    
# Determine start and end point for length measurment for each channel by checking, if the local average intensity is above the defined threshold    
 
 
    for i in df.index:
        LocalSet= normalized_df[(normalized_df.index <= (i + numberpixels)) & (normalized_df.index >= i)]
        LocalAverageCh1 = LocalSet['Ch1'].mean()
        if LocalAverageCh1*100 > threshold:
            print('startindex_Ch1:',i)
            print('LocalAverageCh1:',LocalAverageCh1)
            print('average_valuesCh1:',average_valuesCh1)
            startindex_Ch1 = i
            break
    for i in df.index:
        LocalSet= normalized_df[(normalized_df.index <= (i + numberpixels)) & (normalized_df.index >= i)]
        LocalAverageCh2 = LocalSet['Ch2'].mean()
        if LocalAverageCh2*100 > threshold:
            print('startindex_Ch2:',i)
            print('LocalAverageCh2:',LocalAverageCh2)
            print('average_valuesCh2:',average_valuesCh2)
            startindex_Ch2 = i
            break
            
    for i in df.index:
        LocalSet = normalized_df[(normalized_df.index >= (length-1-i - numberpixels)) & (normalized_df.index <= length-1-i)]
        LocalAverageCh1 = LocalSet['Ch1'].mean()
        if LocalAverageCh1*100 > threshold:
            print('endindex_Ch1:',length-1-i)
            print('LocalAverageCh1:',LocalAverageCh1)
            print('average_valuesCh1:',average_valuesCh1)
            endindex_Ch1 = length-1-i
            break
    for i in df.index:
        LocalSet = normalized_df[(normalized_df.index >= (length-1-i - numberpixels)) & (normalized_df.index <= length-1-i)]
        LocalAverageCh2 = LocalSet['Ch2'].mean()
        if LocalAverageCh2*100 > threshold:
            print('endindex_Ch2:',length-1-i)
            print('LocalAverageCh2:',LocalAverageCh2)
            print('average_valuesCh2:',average_valuesCh2)
            endindex_Ch2 = length-1-i
            break
 # adding the information of distance between manually selected start point (Soma) to start of AIS (or any selected channel)               
    SomaAISDistance = df['Distance'][startindex_Ch1]
    
# calculation of AIS length
    AISLength = df['Distance'][endindex_Ch1] - df['Distance'][startindex_Ch1]
    print('AISLength:', AISLength)

#calculcation of Nav1.6 length
    NaV1_6_Length = df['Distance'][endindex_Ch2] - df['Distance'][startindex_Ch2]
    print('NaV1_6_Length:', NaV1_6_Length)


#  repeat normalization of intensity (work around)  
    normalized_df['Distance'] = df['Distance']
    
    AISSet = normalized_df[(normalized_df.index >= startindex_Ch1) & (normalized_df.index <= endindex_Ch1)]
    
    normalized_AISSet=(AISSet-AISSet.min())/(AISSet.max()-AISSet.min())
    
    
    AISSet['Distance'] = normalized_AISSet['Distance']
 
    
 
#plotting the results of normalized channels

    AIS_plot = normalized_df
    sns.lineplot(data= AIS_plot, x=df["Distance"], y=normalized_df["Ch1"])
    
    normalized_df_rounded = AISSet
    normalized_df_rounded_temp =  round(AISSet, 2)
    normalized_df_rounded["Distance"] = normalized_df_rounded_temp["Distance"] 
    
    Ch2_plot = normalized_df
    sns.lineplot(data= Ch2_plot, x=df["Distance"], y=normalized_df["Ch2"])
    normalized_df_rounded["celltype"] = celltype
    print (normalized_df_rounded)
    
    
    
    normalized_df_rounded.to_csv(path_or_buf=path_calculation, mode='a', header=firsttime)
 
     
    
 # creating csv files containing all acquired information  
 
    with open(Aggregated_data, 'a', encoding="ISO-8859-1", newline='') as myfile:
          wr = csv.writer(myfile)
          if firsttime:
              wr.writerow(("Celltype","startindex_Ch1","endindex_Ch1","startindex_Ch2","endindex_Ch2","AISlength","Channel_length"," SomaAISDistance", "path"))  
          wr.writerow((celltype,
              startindex_Ch1,
            endindex_Ch1,
            startindex_Ch2,
            endindex_Ch2,
            AISLength,
            NaV1_6_Length,
            SomaAISDistance,
            path+file_name)
              )
print(celltype)
