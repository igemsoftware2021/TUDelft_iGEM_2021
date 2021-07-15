import cv2 
import numpy as np

#define video file name
filename = './microfluidic_chip_design/video_analysis/wells_geenwhatsapp.MOV'
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
    _, bi = cv2.threshold(blur, 120, 255, cv2.THRESH_BINARY)


    # find rectangles in image
    # first find the edges to find the lines
    # edges = cv2.Canny(bi, threshold1=30, threshold2=100)
    contours = cv2.findContours(bi, cv2.RETR_TREE, cv2.CHAIN_APPROX_NONE)[0]
    # cv2.drawContours(image, contours,-1,(128,255,0),3)
    # print(contours)
    # Loop over contours to find the rectangles
    cntrrect = []
    cntrwell = []
    for c in contours:
        peri = cv2.arcLength(c, True)
        approx = cv2.approxPolyDP(c, 0.02*peri, True)
        # approx2 = cv2.approxPolyDP(c, 0.02*peri, False)
        (x, y, w, h) = cv2.boundingRect(approx)
        ar = w / float(h)
        area = cv2.contourArea(c)

        # DETECT SQUARES
        if len(approx) == 4 and ar>=0.95 and ar<=1.05 and area > 500:
            cv2.drawContours(image, c, -1, (128,255,0), 2)
            cntrrect.append(c)

        # DETECT CIRCLES (WELLS)
        # elif len(approx) > 5 and area > 250:
        #     cv2.drawContours(image, c, -1, (128,255,0), 2)
        #     cntrwell.append(c)
    # print(cntrrect)        

    # lines = cv2.HoughLinesP(edges, 1, np.pi/180, 60, np.array([]), 50, 5)

     # show image
    cv2.imshow('Frame', image)
    if cv2.waitKey() & 0xFF == ord('q'): #press q to quit
        #waitKey is a method which displays the frame for specified milliseconds
        #‘0xFF == ord(‘q’)’ inside the ‘if’ statement is a special syntax to provide the ‘while’ loop break, by a keyboard key pressing event
        
        break
    else:
        break

    