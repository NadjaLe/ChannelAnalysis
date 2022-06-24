


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
    threshold_Ch2 = 30 #in my case Nav1.6  
    #number of pixels to consider in the local average
    numberpixels = 5
    
 # Variables   
  
    average_valuesCh1 = df.mean()[1]
    average_valuesCh2 = df.mean()[2]
    average_valuesCh3 = df.mean()[3]
    length = len(df.index)
    startindex_Ch1 = None
    endindex_Ch1 = None
    startindex_Ch2 = None
    endindex_Ch2 = None
    
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
    
# Exclusion of further analysis + defining boundaries; in my case: AIS length
  if not(startindex_Ch1 is None or endindex_Ch1 is None):
    AISSet = normalized_df[(normalized_df.index >= startindex_Ch1) & (normalized_df.index <= endindex_Ch1)]
    
 # Create index list for length calculation of Nav1.6 within the AIS.   
       for i in AISSet.index:
            LocalSet= normalized_df[(normalized_df.index <= (i + numberpixels)) & (normalized_df.index >= i)]
            LocalAverageCh2 = LocalSet['Ch2'].mean()
            if LocalAverageCh2*100 > threshold_channel2:
                startindex_Ch2 = i 
                break
                    indexlist = AISSet.index.tolist()
        indexlist.reverse()         
        for i in indexlist:
            LocalSet = normalized_df[(normalized_df.index >= (i - numberpixels)) & (normalized_df.index <= i)]
            LocalAverageCh2 = LocalSet['Ch2'].mean()
            if LocalAverageCh2*100 > threshold_channel2:
                endindex_Ch2 = i
                break
                
# Calculation of Soma- AIS distance:
    if startindex_Ch1 is None:
        SomaAISDistance = None
    else:
         SomaAISDistance = df['Distance'][startindex_Ch1]
                
# calculation of AIS length; Exclusion of further analysis if AIS length is < 10 µm length
    if startindex_Ch1 is None or endindex_Ch1 is None:
        AISLength = None
    else:
     AISLength = df['Distance'][endindex_Ch1] - df['Distance'][startindex_Ch1]   
     if AISLength <= 10:
        AISLength = None

#calculcation of Nav1.6 length
  if startindex_Ch2 is None or endindex_Ch2 is None:
        NaV1_6_Length = None
    else:    
        NaV1_6_Length = df['Distance'][endindex_Ch2] - df['Distance'][startindex_Ch2]
        
 # calculation of distance between AIS start and Nav1.6 start
    if AISLength is None or NaV1_6_Length is None:  
        AIS_Nav1_6_distance = None
    else:   
        AIS_Nav1_6_distance = df ['Distance'][startindex_Ch2] - df['Distance'][startindex_Ch1]
        if startindex_Ch2 == endindex_Ch2:
            AIS_Nav1_6_distance = None
            print ('AIS_Nav1_6_distance', AIS_Nav1_6_distance)


    
    
 
    
 
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
    
    
 # creating csv file which contains all raw information of the channel analysis
    normalized_df_rounded.to_csv(path_or_buf=path_calculation, mode='a', header=firsttime)
 
     
    
 # creating csv files containing all acquired information  
 
    with open(Aggregated_data, 'a', encoding="ISO-8859-1", newline='') as myfile:
          wr = csv.writer(myfile)
          if firsttime:
              wr.writerow(("Celltype","startindex_Ch1","endindex_Ch1","startindex_Ch2","endindex_Ch2","AISlength","Channel_length"," SomaAISDistance","AIS_Nav_distance", "path"))  
          wr.writerow((celltype,
              startindex_Ch1,
            endindex_Ch1,
            startindex_Ch2,
            endindex_Ch2,
            AISLength,
            NaV1_6_Length,
            SomaAISDistance,
            AIS_Nav1_6_distance,
            path+file_name)
              )
print(celltype)
