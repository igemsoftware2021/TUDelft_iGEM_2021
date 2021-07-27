import cv2
import numpy as np
import helper_functions

# This file is to create an info image, so you know what every channel
# number corresponds to
# The file should be changed to a function


# Define video file name
# filename = './microfluidic_chip_design/video_analysis/rechte_wells.MOV'
filename = './microfluidic_chip_design/video_analysis/IMG_1715.MOV'
cap = cv2.VideoCapture(filename)  # load the video

# Retrieve the total number of videos
total_frames = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))

# Set the sample rate
sample_rate = 1

# Detection of wells and squares
# Detection at frame 2
cap.set(cv2.CAP_PROP_POS_FRAMES, 2)
_, image = cap.read()

label_contour_dict, structure_dict = helper_functions.determine_fluidic_part_structure(
    image)

for label in label_contour_dict.keys():
    cnt = label_contour_dict[label]
    M = cv2.moments(cnt)
    cX = int((M["m10"] / M["m00"]))
    cY = int((M["m01"] / M["m00"]))
    cv2.drawContours(image, cnt, -1, (50, 255, 50), 1)
    cv2.putText(image, str(label), (cX+5, cY+5), cv2.FONT_HERSHEY_SIMPLEX,
                0.5, (0, 255, 0), 2)

for label in structure_dict.keys():
    for num in structure_dict[label].keys():
        dim = structure_dict[label][num][0].shape[0]
        x_middle = structure_dict[label][num][1][0]
        y_middle = structure_dict[label][num][0][0]

        color = np.random.randint(0, 256, size=3, dtype=np.int32)
        # Change the int32 values to Python int values, otherwise
        # cv2.putText function does not work
        color = (int(color[0]), int(color[1]), int(color[2]))

        image[structure_dict[label][num]] = color
        cv2.putText(image, str(num), (x_middle+5, y_middle+5), cv2.FONT_HERSHEY_SIMPLEX,
                    0.5, color, 1)


# show the output image
cv2.imshow("Image", image)
cv2.waitKey(0)

# Write the image to a file
# cv2.imwrite(filename, image)
