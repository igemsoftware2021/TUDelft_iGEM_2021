import cv2 
import numpy as np

#define video file name
filename = './microfluidic_chip_design/video_analysis/wells_geenwhatsapp.MOV'
cap = cv2.VideoCapture(filename)     #load the video

total_frames = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
sample_rate = 1

# Detection of wells and squares 
# Detection at frame 2 
cap.set(cv2.CAP_PROP_POS_FRAMES, 2)
_, image = cap.read()

# turn RGB(?) image grey
gray = cv2.cvtColor(image, cv2.COLOR_RGB2GRAY)

blur = cv2.GaussianBlur(gray, (3,3), 0) 

# create binary image
_, bi = cv2.threshold(blur, 120, 255, cv2.THRESH_BINARY)

# find rectangles in image
# first find the edges to find the lines
# edges = cv2.Canny(bi, threshold1=30, threshold2=100)
contours = cv2.findContours(bi, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)[0]
# cv2.drawContours(image, contours,-1,(128,255,0),3)
# print(contours)
# Loop over contours to find the rectangles
cntrrect = []
cntrcircle = []
for c in contours:
    peri = cv2.arcLength(c, True)
    approx = cv2.approxPolyDP(c, 0.02*peri, True)
    # approx2 = cv2.approxPolyDP(c, 0.02*peri, False)
    (x, y, w, h) = cv2.boundingRect(approx)
    ar = w / float(h)
    area = cv2.contourArea(c)

    # DETECT SQUARES
    if len(approx) == 4 and ar>=0.95 and ar<=1.05 and area > 500:
        cv2.drawContours(image, c, -1, (0,255,0), 2)
        cntrrect.append(approx)

# DETECT CIRCLES (WELLS)
circles = cv2.HoughCircles(bi, cv2.HOUGH_GRADIENT, 1.2, 100, param1=50, param2=30, minRadius=10, maxRadius=150)  
circles = np.around(circles).astype("int")
for i in circles[0,:]:
    # draw the outer circle
    cv2.circle(image, (i[0],i[1]),i[2],(0,0,255),2)
    # draw the center of the circle
    cv2.circle(image, (i[0],i[1]),2,(0,0,255),3)

# show image
cv2.imshow('Frame', image)
if cv2.waitKey(0) & 0xFF == ord('q'): #press q to quit
    cv2.destroyAllWindows()



# Iterate over frames and calculate distance
# for i in range(0, total_frames, sample_rate):
#     #get one specific frame number i
#     cap.set(cv2.CAP_PROP_POS_FRAMES, i)
#     _, image2 = cap.read()

#     # turn RGB(?) image grey
#     gray2 = cv2.cvtColor(image2, cv2.COLOR_RGB2GRAY)

#     blur2 = cv2.GaussianBlur(gray2, (3,3), 0) 

#     # create binary image
#     _, bi2 = cv2.threshold(blur2, 120, 255, cv2.THRESH_BINARY)




    

    