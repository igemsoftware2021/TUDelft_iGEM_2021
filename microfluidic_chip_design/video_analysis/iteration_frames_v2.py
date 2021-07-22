import cv2
import numpy as np
import helper_functions

# Define video file name
filename = './microfluidic_chip_design/video_analysis/rechte_wells.MOV'
cap = cv2.VideoCapture(filename)  # load the video

# Retrieve the total number of videos
total_frames = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))

# Set the sample rate
sample_rate = 1

# Detection of wells and squares
# Detection at frame 2
cap.set(cv2.CAP_PROP_POS_FRAMES, 2)
_, image = cap.read()

circles = helper_functions.circle_finder(image)
circles = np.around(circles).astype("int")

for i in circles[0, :]:
    # draw the outer circle
    cv2.circle(image, (i[0], i[1]), i[2], (0, 0, 255), 1)
    # draw the center of the circle
    cv2.circle(image, (i[0], i[1]), 2, (0, 0, 255), 1)


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

# fig, ax = plt.subplots()
# for label in label_store:
#     angles = label_angles_dict[label]
#     ax.scatter(np.full(angles.shape[0], label, dtype=np.int32), angles)

# fig.show()
# plt.show()
