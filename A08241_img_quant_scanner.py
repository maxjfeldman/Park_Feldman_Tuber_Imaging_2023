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
infile_query = pwd + "/*.jpg"
files = (glob.glob(infile_query))
files.sort()
r_names=["r_" + str(s) for s in list(range(0, 256))]
g_names=["g_" + str(s) for s in list(range(0, 256))]
b_names=["b_" + str(s) for s in list(range(0, 256))]

summary_table = pd.DataFrame(columns=['img_name', 'clone', 'rep', 'tuber', 'cmx', 'cmy', 'area', 'perimeter', 'length', 'width', 'ratio', 'eccentricity', 'red_ave', 'green_ave', 'blue_ave', 'red_sd', 'green_sd', 'blue_sd'])

r_table = pd.DataFrame(columns=r_names)
g_table = pd.DataFrame(columns=g_names)
b_table = pd.DataFrame(columns=b_names)

files2=files


for f in files[260:361]:
    img = cv2.imread(f)
    image_name = f.split("/")[-1]
    image_name = image_name.replace(".jpg", "")
    imgC = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
    ## Color correction routine goes here....
    #plt.imshow(imgC)
    #plt.show()
    ix, iy, iz = np.shape(img)
    #imgL = cv2.line(imgC, (1350,0), (1350, 4000), (0,255,0), 10)
    #imgL = cv2.line(imgL, (4700, 0), (4700, 4000), (0, 255, 0),10)
    #imgL = cv2.line(imgL, (0, 1900), (6000, 1900), (0, 0, 255),10)
    #imgL = cv2.line(imgL, (0, 3700), (6000, 3700), (0, 0, 255), 10)
    #215: 2600, 1500: 5050
    #imgM = imgC[0:1900,1200:4700]
    ix, iy, iz = np.shape(imgC)
    blur = cv2.GaussianBlur(imgC, (5, 5), 3)
    imgHSV = cv2.cvtColor(blur, cv2.COLOR_BGR2HSV)
    imgLAB = cv2.cvtColor(blur, cv2.COLOR_BGR2Lab)
    h, s, v = cv2.split(imgHSV)
    l, a, b = cv2.split(imgLAB)
    r = blur[:, :, 0]
    g = blur[:, :, 1]
    bl = blur[:, :, 2]
    #plt.hist(v.ravel(),256,[0,256]); plt.show()
    v_ret, v_th = cv2.threshold(v, 50, 255, cv2.THRESH_BINARY)
    #plt.imshow(v_th)
    #plt.show()
    b_ret_inv, b_th_inv = cv2.threshold(b, 125, 255, cv2.THRESH_BINARY_INV)
    #plt.hist(l.ravel(),256,[0,256]); plt.show()
    l_ret, l_th = cv2.threshold(l, 50, 255, cv2.THRESH_BINARY)
    #plt.imshow(l_th)
    #plt.show()
    r_ret, r_th = cv2.threshold(r, 50, 255, cv2.THRESH_BINARY)
    #plt.imshow(r_th)
    #plt.show()
    ## Maybe a is a good channel for the chip
    chip_ret, chip_th = cv2.threshold(a, 140, 255, cv2.THRESH_BINARY)
    mask = cv2.bitwise_and(v_th, b_th_inv)
    kernel = np.ones((5, 5), np.uint8)
    mask_er = cv2.erode(mask, kernel, iterations=5)
    chip_er = cv2.erode(chip_th, kernel, iterations=3)
    mask_dil = cv2.dilate(mask_er, kernel, iterations=3)
    chip_dil = cv2.dilate(chip_er, kernel, iterations=5)
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
    df = df.sort_values(by='cmy', ascending=True)
    top_side = df[:3]
    top_side = top_side.sort_values(by='cmx', ascending=False)
    bottom_side = df[3:5]
    bottom_side = bottom_side.sort_values(by='cmx', ascending=False)
    cnt_order = top_side['cnt'].tolist()
    cnt_order = cnt_order + bottom_side['cnt'].tolist()
    cnt_order.append(5)
    img_draw = imgC.copy()
    contours2.append(chip_contour)
    ## This subroutine is where the values are measured
    summary, r_values, g_values, b_values = get_A08241_measurements_scanner(cnt_order, contours2, img_draw, image_name)
    summary_table =  summary_table.append(summary, )
    r_table =  r_table.append(r_values, )
    g_table =  g_table.append(g_values, )
    b_table =  b_table.append(b_values, )

summary_table = pd.concat([summary_table, r_table], axis =1)
summary_table = pd.concat([summary_table, g_table], axis =1)
summary_table = pd.concat([summary_table, b_table], axis =1)
# out_table_path = outfile_path + "/test_potato_measurements.csv"
out_table_path = pwd + "/A08241_potato_measurements_scanner_2022-06-09.csv"
summary_table.to_csv(out_table_path, mode='a', header=True, encoding='utf-8')

