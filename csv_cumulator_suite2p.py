
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 19:35:50 2022

@author: BioCraze
"""
import csv
import os




files = [i.path for i in os.scandir("E:/BAPTA NE/BAPTA_spont/analysis/data-peaks/") if i.path.endswith('PEAKS.csv')]

def cumulate(files, output):
    cumulative = []
    for file in files:
        f = open(file)
        read = csv.reader(f)
        for row in read:
            row.insert(0, file)
            cumulative.append(row)
    
            
    
    new_csv = open(output,'w',newline='')
    
    write_to = csv.writer(new_csv)
    for row in cumulative:
        write_to.writerows([row])
    
    new_csv.close()

cumulate(files, "")


def count_active(files, output):
    active_count=[]
    for file in files:
        f = open(file)
        read = csv.reader(f)
        count=sum([1 for i in read if len(i)>1]) #subtracting 1 from length of row because each row starts with the roi idx
        animal=file[file.rfind("/")+1:file.rfind("_")]
        if "control" in animal:
            treatment="control"
        else:
            treatment="bapta"
            
        active_count.append([animal,treatment,count])
        
    new_csv = open(output,'w',newline='')
    
    write_to = csv.writer(new_csv)
    for row in active_count:
        write_to.writerows([row])
    
    new_csv.close()
        
count_active(files, "E:/BAPTA NE/BAPTA_spont/summaries/bapta_spont_active_glia_count.csv")
    
    