import cv2 
import numpy as np

#define video file name
filename = './microfluidic_chip_design/video_analysis/test_wells.mp4'
cap = cv2.VideoCapture(filename)     #load the video

total_frames = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
sample_rate = 1
for i in range(0, total_frames, sample_rate):
    #get one specific frame number i
    cap.set(cv2.CAP_PROP_POS_FRAMES, i)
    _, image = cap.read()

    # turn RGB(?) image grey
    gray = cv2.cvtColor(image, cv2.COLOR_RGB2GRAY)

    blur = cv2.GaussianBlur(gray, (3,3), 0) 

    # create binary image
    _, bi = cv2.threshold(gray, 80, 255, cv2.THRESH_BINARY)

     # show gray image
    cv2.imshow('Frame', bi)
    if cv2.waitKey(0) & 0xFF == ord('q'): #press q to quit
        #waitKey is a method which displays the frame for specified milliseconds
        #‘0xFF == ord(‘q’)’ inside the ‘if’ statement is a special syntax to provide the ‘while’ loop break, by a keyboard key pressing event
        
        break
    else:
        break

    