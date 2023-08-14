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
import scipy
pwd = os.getcwd()
infile_query = pwd + "/*d}.tif"
files = (glob.glob(infile_query))
files.sort()
r_names=["r_" + str(s) for s in list(range(0, 256))]
g_names=["g_" + str(s) for s in list(range(0, 256))]
b_names=["b_" + str(s) for s in list(range(0, 256))]
x_names = ["x_" + str(x) for x in list(range(1, 101))]
y_names = ["y_" + str(y) for y in list(range(1, 101))]
s_names = x_names + y_names

summary_table = pd.DataFrame(columns=['img_name', 'clone', 'rep', 'side', 'light', 'tuber', 'cmx', 'cmy', 'area', 'perimeter', 'length', 'width', 'ratio', 'eccentricity', 'red_ave', 'green_ave', 'blue_ave', 'red_sd', 'green_sd', 'blue_sd'])

r_table = pd.DataFrame(columns=r_names)
g_table = pd.DataFrame(columns=g_names)
b_table = pd.DataFrame(columns=b_names)
s_table = pd.DataFrame(columns=s_names)

for f in files:
    img = cv2.imread(f)
    image_name = f.split("/")[-1]
    image_name = image_name.replace(".tif", "")
    image_name = image_name.split("}")[0]
    image_name = image_name.replace("{fileName=", "")
    imgC = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
    #plt.imshow(imgC)
    #plt.show()
    ix, iy, iz = np.shape(img)
    #imgL = cv2.line(imgC, (1350,0), (1350, 4000), (0,255,0), 10)
    #imgL = cv2.line(imgL, (4700, 0), (4700, 4000), (0, 255, 0),10)
    #imgL = cv2.line(imgL, (0, 1900), (6000, 1900), (0, 0, 255),10)
    #imgL = cv2.line(imgL, (0, 3700), (6000, 3700), (0, 0, 255), 10)
    #215: 2600, 1500: 5050
    #imgM = imgC[0:1900,1200:4700]
    img_crop = imgC[0:1900, 1200:3500]
    ix, iy, iz = np.shape(img)
    blur = cv2.GaussianBlur(img_crop, (5, 5), 3)
    imgHSV = cv2.cvtColor(blur, cv2.COLOR_BGR2HSV)
    imgLAB = cv2.cvtColor(blur, cv2.COLOR_BGR2Lab)
    h, s, v = cv2.split(imgHSV)
    l, a, b = cv2.split(imgLAB)
    r = blur[:, :, 0]
    g = blur[:, :, 1]
    bl = blur[:, :, 2]
    #s_ret, s_th = cv2.threshold(s, 100, 255, cv2.THRESH_BINARY)
    s_ret, s_th = cv2.threshold(s, 65, 255, cv2.THRESH_BINARY)
    b_ret_inv, b_th_inv = cv2.threshold(b, 125, 255, cv2.THRESH_BINARY_INV)
    chip_ret, chip_th = cv2.threshold(b, 140, 255, cv2.THRESH_BINARY)
    mask = cv2.bitwise_and(s_th, b_th_inv)
    kernel = np.ones((5, 5), np.uint8)
    mask_er = cv2.erode(mask, kernel, iterations=1)
    chip_er = cv2.erode(chip_th, kernel, iterations=1)
    mask_dil = cv2.dilate(mask_er, kernel, iterations=1)
    chip_dil = cv2.dilate(chip_er, kernel, iterations=1)
    contours, hierarchy = cv2.findContours(mask_dil, cv2.RETR_TREE, cv2.cv2.CHAIN_APPROX_NONE)
    chip_contour, chip_hierarchy = cv2.findContours(chip_dil, cv2.RETR_TREE, cv2.cv2.CHAIN_APPROX_NONE)
    #chip_solid = cv2.drawContours(chip_dil, chip_contour, -1, (255,255,255), -1)
    #chip_s_contour, chip_s_hierarchy = cv2.findContours(chip_dil, cv2.RETR_TREE, cv2.cv2.CHAIN_APPROX_NONE)
    chipSorted = sorted(chip_contour, key=lambda x: cv2.contourArea(x), reverse=True)
    chip_contour = chipSorted[0]
    contours2 = []
    it = 0
    cXY = []
    cntsSorted = sorted(contours, key=lambda x: cv2.contourArea(x), reverse=True)
    contours = cntsSorted[0:5]
    for cnt in contours:
        area = cv2.contourArea(cnt)
        if area > 3000:
            contours2.append(cnt)
            M = cv2.moments(cnt)
            cX = int(M["m10"] / M["m00"])
            cY = int(M["m01"] / M["m00"])
            #print(cX)
            #print(cY)
            entry = [it, cX, cY]
            #print(entry)
            cXY.append(entry)
            it += 1
    df = pd.DataFrame(cXY, columns=['cnt', 'cmx', 'cmy'])
    df = df.sort_values(by='cmx', ascending=False)
    right_side = df[:3]
    right_side = right_side.sort_values(by='cmy')
    left_side = df[3:5]
    left_side = left_side.sort_values(by='cmy')
    cnt_order = right_side['cnt'].tolist()
    cnt_order = cnt_order + left_side['cnt'].tolist()
    cnt_order.append(5)
    img_draw = img_crop.copy()
    contours2.append(chip_contour)
    ## This subroutine is where the values are measured
    summary, r_values, g_values, b_values, s_values = get_A08241_measurements(cnt_order, contours2, img_draw, image_name)
    #summary, r_values, g_values, b_values = get_A08241_measurements(cnt_order, contours2, img_draw, image_name)
    summary_table =  summary_table.append(summary, )
    r_table =  r_table.append(r_values, )
    g_table =  g_table.append(g_values, )
    b_table =  b_table.append(b_values, )
    s_table =  s_table.append(s_values, )

summary_table = pd.concat([summary_table, r_table], axis =1)
summary_table = pd.concat([summary_table, g_table], axis =1)
summary_table = pd.concat([summary_table, b_table], axis =1)
summary_table = pd.concat([summary_table, s_table], axis =1)
# out_table_path = outfile_path + "/test_potato_measurements.csv"
out_table_path = pwd + "/A08241_potato_measurements_std_shape_ColorCorrected_2022-02-09.csv"
summary_table.to_csv(out_table_path, mode='a', header=True, encoding='utf-8')
