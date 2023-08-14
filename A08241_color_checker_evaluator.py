## Import python libraries
from time import sleep
import csv
import numpy
from pylab import *
from datetime import date
import os.path
import matplotlib.pyplot as plt
import cv2
import sys
from optparse import OptionParser
import re
import glob
import pandas as pd
import skimage
from skimage import feature
#import plotnine
import math
import scipy

pwd = os.getcwd()
infile_query = pwd + "/*d}.tif"
#infile_query = pwd + "/*.jpg"
files = (glob.glob(infile_query))
files.sort()

## Lets make a header for the output
r_names=["r_" + str(s) for s in list(range(1, 25))]
g_names=["g_" + str(s) for s in list(range(1, 25))]
b_names=["b_" + str(s) for s in list(range(1, 25))]

## Here is the header
output = ["img_name"] + r_names + g_names + b_names

for f in files:
    print(f)
    ## source
    img = cv2.imread(f)
    img = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
    #plt.imshow(img)
    #plt.show()
    ## Lets crop the image box
    #img = img[1900:3700, 1650:4400]
    ## Lets crop the image std black background
    #img = img[0:1900,1000:5000]
    ## Scanner doesn't need to be cropped
    #plt.imshow(img)
    #plt.show()
    df, start_coord, spacing = find_color_card(img, 'adaptgauss', 125, False, 'light')
    if spacing == 0:
        continue
    imgC_copy = rotate_color_card(df, spacing, img)
    df, start_coord, spacing = find_color_card(imgC_copy, 'adaptgauss', 125, False, 'light')
    if spacing == 0:
        continue
    ## spacing = 260
    #spacing = (260, 260)
    mask = create_color_card_mask(rgb_img=imgC_copy, radius=5, start_coord=start_coord, spacing=spacing, ncols=4, nrows=6, exclude=[])
    #mask = create_color_card_mask(rgb_img=imgC_copy, radius=5, start_coord=start_coord, spacing=spacing, ncols=6, nrows=4, exclude=[])
    #plt.imshow(mask,cmap='gray', vmin=0, vmax=255)
    #plt.show()
    maskedImg = cv2.bitwise_and(imgC_copy, imgC_copy, mask=mask)
    #plt.imshow(maskedImg)
    #plt.show()
    col_headers, col_matrix = get_color_matrix(rgb_img=imgC_copy, mask=mask)
    col_matrix = orient_color_card(col_matrix, ncols = 4)
    #col_matrix = orient_color_card(col_matrix, ncols=6)
    img_name_list = f.split('/')
    img_name = img_name_list[-1]
    col_matrix = col_matrix[:,1:]
    entry = []
    entry = np.hstack((entry, img_name))
    for col in col_matrix.T:
        entry = np.hstack((entry, col))
        #print(entry)
    output = np.vstack((output, entry))
out_table_path = pwd + "/color_checker_data_box.csv"
#out_table_path = pwd + "/color_checker_data_std.csv"
#out_table_path = pwd + "/color_checker_data_scanner_rotated.csv"
#out_table_path = pwd + "/color_checker_data_scanner.csv"


np.savetxt(out_table_path, output, fmt='%s', delimiter=",")


