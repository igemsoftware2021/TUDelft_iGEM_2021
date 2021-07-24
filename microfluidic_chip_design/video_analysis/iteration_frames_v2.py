import cv2
import numpy as np
import helper_functions

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

dimage = np.copy(image)

# show image
cv2.imshow('Frame', image)
if cv2.waitKey(0) & 0xFF == ord('q'):  # press q to quit
    cv2.destroyAllWindows()

# image_gray = cv2.cvtColor(image, cv2.COLOR_RGB2GRAY)
# image_thres = cv2.adaptiveThreshold(image_gray, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C,
#                                     cv2.THRESH_BINARY, 5, 6)

# # show image
# cv2.imshow('Frame', image_thres)
# if cv2.waitKey(0) & 0xFF == ord('q'):  # press q to quit
#     cv2.destroyAllWindows()

cimage = np.copy(image)
conversion_factor = helper_functions.find_conversion_factor(
    image, length_squares=3.5)  # mm/pixel

labels, label_store, label_contour_dict = helper_functions.find_fluidic_components_and_contours(
    cimage)

for label in label_contour_dict.keys():
    cv2.drawContours(cimage, label_contour_dict[label], -1, (0, 255, 0), 2)

# show image
cv2.imshow('Cimage', cimage)
if cv2.waitKey(0) & 0xFF == ord('q'):  # press q to quit
    cv2.destroyAllWindows()

# show image
cv2.imshow('Image_iteration', dimage)
if cv2.waitKey(0) & 0xFF == ord('q'):  # press q to quit
    cv2.destroyAllWindows()
circles = helper_functions.circle_finder(dimage)
circles = np.around(circles).astype("int")

for i in circles[0, :]:
    # draw the outer circle
    cv2.circle(image, (i[0], i[1]), i[2], (0, 0, 255), 2)
    # draw the center of the circle
    cv2.circle(image, (i[0], i[1]), 1, (0, 0, 255), 1)

temp = helper_functions.determine_fluidic_part_structure(dimage)
print(temp)
# show image
cv2.imshow('Cimage', cimage)
if cv2.waitKey(0) & 0xFF == ord('q'):  # press q to quit
    cv2.destroyAllWindows()

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
