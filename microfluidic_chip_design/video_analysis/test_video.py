import cv2 
import numpy as np
#import matplotlib

#define video file name
filename = './microfluidic_chip_design/video_analysis/test_wells.mp4'
cap = cv2.VideoCapture(filename)     #load the video

# count = 0
# print(type(cap))

# Check if camera opened successfully
if (cap.isOpened()== False):
    print("Error opening video stream or file")

#create a non-trivial cycle: read until video is completed
while cap.isOpened():               #play the video by reading frame by frame
    ret, frame = cap.read()
    if ret == True:
        
        #optional: image processing here
        
        cv2.imshow('Frame', frame)       #shows the video
        # cv2.imwrite("frame%d.jpg" % count, frame)
        # count = count + 1

    if cv2.waitKey(25) & 0xFF == ord('q'): #press q to quit
        #waitKey is a method which displays the frame for specified milliseconds
        #‘0xFF == ord(‘q’)’ inside the ‘if’ statement is a special syntax to provide the ‘while’ loop break, by a keyboard key pressing event
        
        break
    else:
        break

# print(frame)

# cap.release()
# cv2.destroyAllWindows()
# print(count)
# print('loulotte')
"""
8
# Check if camera opened successfully
9
if (cap.isOpened()== False):
10
  print("Error opening video stream or file")
11
 
12
# Read until video is completed
13
while(cap.isOpened()):
14
  # Capture frame-by-frame
15
  ret, frame = cap.read()
16
  if ret == True:
17
 
18
    # Display the resulting frame
19
    cv2.imshow('Frame',frame)
20
 
21
    # Press Q on keyboard to  exit
22
    if cv2.waitKey(25) & 0xFF == ord('q'):
23
      break
24
 
25
  # Break the loop
26
  else:
27
    break
28
 
29
# When everything done, release the video capture object
30
cap.release()
31
 
32
# Closes all the frames
33
cv2.destroyAllWindows()
"""
# """
# count = 0
# while cap.isOpened():
#     ret,frame = cap.read()
#     cv2.imshow('window-name', frame)
#     cv2.imwrite("frame%d.jpg" % count, frame)
#     count = count + 1
#     if cv2.waitKey(10) & 0xFF == ord('q'):
#         break
# """

frame_width = 10
frame_height = 8
# Define the codec and create VideoWriter object.The output is stored in 'outpy.avi' file.
# Define the fps to be equal to 10. Also frame size is passed.
out = cv2.VideoWriter('output_test.avi', cv2.VideoWriter_fourcc('M','J','P','G'), 10, (frame_width,frame_height))
