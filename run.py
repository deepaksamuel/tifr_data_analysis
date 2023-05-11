#use this code to run on all files in a given folder
#python3 run.py
import os
dir_list = os.listdir("/DATA/tifr_data_clsf")
print(dir_list)
for file in dir_list:
    if "ire" in file:
        print(file)
        os.system("./analysis /DATA/tifr_data_clsf/{}".format(file))
