import cv2 
import numpy as np

#define video file name
filename = './microfluidic_chip_design/video_analysis/test_wells.mp4'
cap = cv2.VideoCapture(filename)     #load the video

total_frames = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
sample_rate = 5
for i in range(0, total_frames, sample_rate):
    #get one specific frame number i
    cap.set(cv2.CAP_PROP_POS_FRAMES, i)
    _, image = cap.read()
    frametest = cv2.imshow('Frame', image)
    print(frametest)

