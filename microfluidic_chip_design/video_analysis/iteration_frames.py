import cv2
import numpy as np

# define video file name
filename = './microfluidic_chip_design/video_analysis/rechte_wells.MOV'
cap = cv2.VideoCapture(filename)  # load the video

total_frames = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
sample_rate = 1

# Detection of wells and squares
# Detection at frame 2
cap.set(cv2.CAP_PROP_POS_FRAMES, 2)
_, image = cap.read()

cimage = np.copy(image)

# turn RGB(?) image grey
gray = cv2.cvtColor(image, cv2.COLOR_RGB2GRAY)

blur = cv2.GaussianBlur(gray, (3, 3), 0)

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
cntrcircle = []
cntrchannel = []
cntrwell = []
for c in contours:
    peri = cv2.arcLength(c, True)
    approx = cv2.approxPolyDP(c, 0.02*peri, True)
    # approx2 = cv2.approxPolyDP(c, 0.02*peri, False)
    (x, y, w, h) = cv2.boundingRect(approx)
    ar = w / float(h)
    area = cv2.contourArea(c)

    # DETECT SQUARES
    if len(approx) == 4 and ar >= 0.95 and ar <= 1.05 and area > 500:
        cv2.drawContours(image, c, -1, (0, 255, 0), 2)
        cntrrect.append(approx)

# DETECT CIRCLES (WELLS)
circles = cv2.HoughCircles(bi, cv2.HOUGH_GRADIENT_ALT, 1.2,
                           5, param1=150, param2=0.9, minRadius=10, maxRadius=0)
circles = np.around(circles).astype("int")
for i in circles[0, :]:
    # draw the outer circle
    cv2.circle(image, (i[0], i[1]), i[2], (0, 0, 255), 2)
    # draw the center of the circle
    cv2.circle(image, (i[0], i[1]), 2, (0, 0, 255), 3)

output = cv2.connectedComponentsWithStats(bi, connectivity=8, ltype=cv2.CV_32S)
(num_labels, labels, stats, centroids) = output
# output = image.copy()
# find height and width of all connected components
cntrmidleft = []
cntrmidright = []
label_store = []
for j in range(2, num_labels):  # START FROM 2, CAUSE INDEX 0 AND 1 ARE BACKGROUND
    X = stats[j, cv2.CC_STAT_LEFT]
    Y = stats[j, cv2.CC_STAT_TOP]
    W = stats[j, cv2.CC_STAT_WIDTH]
    H = stats[j, cv2.CC_STAT_HEIGHT]
    area2 = stats[j, cv2.CC_STAT_AREA]
    (cX, cY) = centroids[j]
    if W > 500:
        cv2.circle(image, (int(cX), int(cY)), 4, (0, 0, 255), -1)
        cv2.rectangle(image, (X, Y), (X+W, Y+H), (0, 0, 255), 3)
        #rect = cv2.minAreaRect()
        cntrchannel.append((int(cX), int(cY)))
        cntrwell.append((X, Y))
        Ymid = Y + H/2
        cntrmidleft.append((X, int(Ymid)))
        cntrmidright.append(((X+W), int(Ymid)))
        label_store.append(j)

print(label_store)

label_contour_dict = {}

# https://stackoverflow.com/questions/37745274/opencv-find-perimeter-of-a-connected-component
# Loop over the label values of the wells with channels
for label in label_store:
    # turn RGB(?) image grey
    cgray = cv2.cvtColor(cimage, cv2.COLOR_RGB2GRAY)

    cblur = cv2.GaussianBlur(cgray, (3, 3), 0)

    # create binary image
    _, cbi = cv2.threshold(cblur, 120, 255, cv2.THRESH_BINARY)
    # Create a mask for the image
    mask_label = (labels == label)

    mask_label_invert = np.invert(mask_label)

    cbi[mask_label_invert] = 0

    contours, hierarchy = cv2.findContours(
        cbi, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)

    temp = contours[0]

    label_contour_dict[label] = contours[0]

# https://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python/13849249#13849249


def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)


def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


label_switch_point_dict = {}
# Vector range to create vectors
vec_range = 11

# Now check where the circles are
for label in label_store:
    contour = label_contour_dict[label]
    contour_adapted = np.concatenate(
        (contour[-vec_range:, :, :], contour, contour[:vec_range, :, :]), axis=0)

    # Create an array and fill it with the value False
    switch_point = np.full(contour.shape[0], False, dtype=bool)

    # Start with vec_range, because you want to start at the first value from contour
    for i in range(vec_range, contour.shape[0]+vec_range):
        # Point is contour[i][0][:]

        # Middle point
        m_x, m_y = contour_adapted[i][0]
        # Before point
        b_x, b_y = contour_adapted[i-vec_range][0]
        # After point
        a_x, a_y = contour_adapted[i+vec_range][0]

        # Vector from middle point to before point
        vector_1 = np.array([b_x - m_x, b_y - m_y], dtype=np.float32)

        # Vector from middle point to after point
        vector_2 = np.array([a_x - m_x, a_y - m_y], dtype=np.float32)

        angle = angle_between(vector_1, vector_2) * (180.0/np.pi)

        if angle < 160:
            switch_point[i-vec_range] = True

    label_switch_point_dict[label] = switch_point


for label in label_store:
    contour = label_contour_dict[label]
    switch_point = label_switch_point_dict[label]

    for i in range(contour.shape[0]):
        # switch_point contains only True and False, so if True
        if switch_point[i]:
            x, y = contour[i][0]
            cimage[y, x] = [0, 0, 255]

# show image
cv2.imshow('Switch points', cimage)
if cv2.waitKey(0) & 0xFF == ord('q'):  # press q to quit
    cv2.destroyAllWindows()


# for k in range(0, len(cntrchannel)):
#     left = cv2.line(image, cntrmidleft[k], cntrchannel[k], (255, 0, 0), 2)
#     right = cv2.line(image, cntrmidright[k], cntrchannel[k], (255, 0, 0), 2)

    # Draw lines in binary image
    # bileft = cv2.line(bi, cntrmidleft[k], cntrchannel[k], (255, 0, 0), 2)
    # biright = cv2.line(bi, cntrmidright[k], cntrchannel[k], (255, 0, 0), 2)
colorcodes = []
# for v in left.shape:

#     if bi[left[0:1][v],0] == 255:
#         colorcodes.append[(left[0:1][v])]


# show image
cv2.imshow('Frame', image)
if cv2.waitKey(0) & 0xFF == ord('q'):  # press q to quit
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

#     bi2[left, 0]
