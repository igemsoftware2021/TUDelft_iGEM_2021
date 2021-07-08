import numpy as np

#detection of shape
def shapedetection(cntr, sp):
    # cntr is list of cntrs
    # sp is the wanted shapes 3 = triangle, 4 = rectangle, 5 = pentagon, 0 = circle
    for i in cntr:
        peri = cv2.arclength(c,True)
