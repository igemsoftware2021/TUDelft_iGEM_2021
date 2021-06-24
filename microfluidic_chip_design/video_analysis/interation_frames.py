import cv2 
import numpy as np

#define video file name
filename = './microfluidic_chip_design/video_analysis/test_wells.mp4'
cap = cv2.VideoCapture(filename)     #load the video

total_frames = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
sample_rate = 
for i in range(0, total_frames, sample_rate)