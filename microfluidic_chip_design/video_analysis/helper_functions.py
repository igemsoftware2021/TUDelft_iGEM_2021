import cv2
import numpy as np


def find_conversion_factor(image, length_squares, kernel_size=3, threshold_min=120, threshold_max=255, min_area=500):
    """
    Paramaters
    ----------
    image: ....
    length_squares (float): real size in mm of the squares in paper chip
    kernel_size (int): size of the Gaussian kernel size


    Returns
    -------
    conversion_factor: (float) in mm / pixel which can be used to calculate real distance in
        milimeters from pixel distance
    """
    # turn RGB(?) image grey
    gray = cv2.cvtColor(image, cv2.COLOR_RGB2GRAY)

    # blur the gray image
    blur = cv2.GaussianBlur(gray, (kernel_size, kernel_size), 0)

    # create binary image
    _, bi = cv2.threshold(blur, threshold_min,
                          threshold_max, cv2.THRESH_BINARY)

    # find contours in image
    contours = cv2.findContours(bi, cv2.RETR_TREE, cv2.CHAIN_APPROX_NONE)[0]

    # create empty list for storing the contours
    reccontour = []

    # loop over contours to find the squares
    for c in contours:
        peri = cv2.arcLength(c, True)
        approx = cv2.approxPolyDP(c, 0.02*peri, True)
        (x, y, w, h) = cv2.boundingRect(approx)
        ar = w / float(h)
        area = cv2.contourArea(c)

        # DETECT SQUARES
        if len(approx) == 4 and ar >= 0.95 and ar <= 1.05 and area > min_area:
            cv2.drawContours(image, c, -1, (0, 255, 0), 2)
            reccontour.append(c)

    # calculate the perimeter of the 4 detected squares and divide by 4 to obtain the width
    width_square = np.zeros(len(reccontour))
    for i in range(len(reccontour)):
        square = reccontour[i]
        peri_i = cv2.arcLength(square, True)
        width_square[i] = peri_i / 4

    # Calculate the average width and the conversion factor, multiply pixel distance by conversion factor
    # to obtain the length in mm
    width_av = np.mean(width_square)  # pixel
    conversion_factor = length_squares / width_av  # mm/pixel

    return conversion_factor


def update_labels(labels, label_store, label_contour_dict):
    """
    Function that updates the labels array and the label store function,
    labels could now be [30,31,32,33,60,61,62,63].
    After the function it will be [1,2,3,4,5,6,7,8]

    """

    labels_updated = np.zeros(labels.shape, dtype=np.int16)
    label_store_updated = []
    label_contour_dict_updated = {}

    new_label = 1
    for label in label_store:
        # Create a mask for the image
        mask_label = (labels == label)
        # Update the labels array
        labels_updated[mask_label] = new_label

        label_store_updated.append(new_label)

        contour = label_contour_dict.get(label)
        label_contour_dict_updated[new_label] = contour

        new_label += 1

    return labels_updated, label_store_updated, label_contour_dict_updated


def find_fluidic_components_and_contours(image, kernel_size=3, threshold_min=120, threshold_max=255, min_area=500):
    """
    Paramaters
    ----------
    image: ....
    kernel_size (int): size of the Gaussian kernel size


    Returns
    -------
    labels: (nd.array) array where every pixel has a label value
    label_store: labels that correspond to a fluidic part
    label_contour_dict: (dict) key is a label value and value is a array of the contour points
    """
    # turn RGB into gray image
    gray = cv2.cvtColor(image, cv2.COLOR_RGB2GRAY)

    # Blur the gray image
    gray = cv2.GaussianBlur(gray, (kernel_size, kernel_size), 0)

    # create binary image of the blurred image
    _, binary = cv2.threshold(
        gray, threshold_min, threshold_max, cv2.THRESH_BINARY)

    # Find the connected components
    num_labels, labels, stats, centroids = cv2.connectedComponentsWithStats(
        binary, connectivity=8, ltype=cv2.CV_32S)

    label_store = []
    for i in range(2, num_labels):  # START FROM 2, CAUSE INDEX 0 AND 1 ARE BACKGROUND
        W = stats[i, cv2.CC_STAT_WIDTH]
        H = stats[i, cv2.CC_STAT_HEIGHT]
        area = stats[i, cv2.CC_STAT_AREA]  # Number of pixels of component

        # component_ratio is the ratio of all the pixels within the bounding box
        # of a connected component that actually belong to a connected component
        component_ratio = area / (W*H)

        # H/W is to filter out the squares around the microfluidic parts
        if area > min_area and component_ratio < 0.9:
            label_store.append(i)

    # Dictionary where the key is the label number and the value is the contour array
    label_contour_dict = {}
    # Dictionary where the key is the label number and the value is the array filled
    # with True and False

    # https://stackoverflow.com/questions/37745274/opencv-find-perimeter-of-a-connected-component
    # Loop over the label values of the wells with channels
    for label in label_store:
        # Copy the binary image
        cbinary = np.copy(binary)

        # Create a mask for the image
        mask_label = (labels == label)

        mask_label_invert = np.invert(mask_label)

        cbinary[mask_label_invert] = 0

        contours, hierarchy = cv2.findContours(
            cbinary, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)

        label_contour_dict[label] = contours[0]

    labels_updated, label_store_updated, label_contour_dict_updated = update_labels(
        labels, label_store, label_contour_dict)

    return labels_updated, label_store_updated, label_contour_dict_updated


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


def circle_finder(image, kernel_size=3, threshold_min=120, threshold_max=255, min_area=500, vector_range=9, max_angle=173):
    cv2.imshow('Image_function', image)
    if cv2.waitKey(0) & 0xFF == ord('q'):  # press q to quit
        cv2.destroyAllWindows()

    labels, label_store, label_contour_dict = find_fluidic_components_and_contours(
        image, kernel_size=kernel_size, threshold_min=threshold_min, threshold_max=threshold_max, min_area=min_area)

    # Dictionary where the key is the label number and the value is the array filled
    # with True and False
    label_switch_point_dict = {}

    # Now check where the circles are
    for label in label_contour_dict.keys():
        contour = label_contour_dict[label]
        contour_adapted = np.concatenate(
            (contour[-vector_range:, :, :], contour, contour[:vector_range, :, :]), axis=0)

        # Create an array and fill it with the value False
        switch_point = np.full(contour.shape[0], False, dtype=bool)

        # Start with vector_range, because you want to start at the first value from contour
        for i in range(vector_range, contour.shape[0]+vector_range):
            # Point is contour[i][0][:]

            # Middle point
            m_x, m_y = contour_adapted[i][0]
            # Before point
            b_x, b_y = contour_adapted[i-vector_range][0]
            # After point
            a_x, a_y = contour_adapted[i+vector_range][0]

            # Vector from middle point to before point
            vector_1 = np.array([b_x - m_x, b_y - m_y], dtype=np.float32)

            # Vector from middle point to after point
            vector_2 = np.array([a_x - m_x, a_y - m_y], dtype=np.float32)

            angle = angle_between(vector_1, vector_2) * (180.0/np.pi)

            if angle < max_angle:
                switch_point[i-vector_range] = True

        label_switch_point_dict[label] = switch_point

    # Clean up the switch point array
    # You look backward and forward from a certain switch_point that is True
    # and see how many switch points behind and infront are also True.
    # If a certain threshold value is reached than you keep storing it as
    # a switch point, otherwise set the value on False
    thres_value = 0.75*vector_range
    for label in label_switch_point_dict.keys():
        switch_point = label_switch_point_dict[label]
        switch_point_adapted = np.concatenate(
            (switch_point[-vector_range:], switch_point, switch_point[:vector_range]), axis=0)

        clean_switch_point = np.full(switch_point.shape[0], False, dtype=bool)
        for i in range(vector_range, switch_point.shape[0]+vector_range):
            # If the switch point already has a True value
            if switch_point_adapted[i]:
                before_part = switch_point_adapted[(i-vector_range):i]
                after_part = switch_point_adapted[i:(i+vector_range)]
                # Booleans are treated as 0 and 1 in arithmic operations, hence we can use .sum()
                if (before_part.sum() < thres_value and after_part.sum() > thres_value) or (before_part.sum() > thres_value and after_part.sum() < thres_value):
                    clean_switch_point[i-vector_range] = False
                else:
                    clean_switch_point[i-vector_range] = True

        label_switch_point_dict[label] = clean_switch_point

    # An image where all the points are black except for
    # the switch points, which will be white
    switch_point_image = np.zeros(
        (image.shape[0], image.shape[1]), dtype=np.uint8)

    for label in label_contour_dict.keys():
        contour = label_contour_dict[label]
        switch_point = label_switch_point_dict[label]

        for i in range(contour.shape[0]):
            # switch_point contains only True and False, so if True
            if switch_point[i]:
                x, y = contour[i][0]
                switch_point_image[y, x] = 255

    # # show image
    # cv2.imshow('switch_point_image', switch_point_image)
    # if cv2.waitKey(0) & 0xFF == ord('q'):  # press q to quit
    #     cv2.destroyAllWindows()

    # Please DO NOT mess with these settings
    circles = cv2.HoughCircles(switch_point_image, cv2.HOUGH_GRADIENT, 1, 15,
                               param1=100, param2=12.5, minRadius=5, maxRadius=75)

    circles = np.around(circles).astype("int")

    return circles


def determine_fluidic_part_structure(image, kernel_size=3, threshold_min=120, threshold_max=255):

    # Structure:
    # {label #num: {channel #num: indices of the line, 2: indices of the line, 3: indices of the line, ...}, label #num: {...}}
    #

    # turn RGB into gray image
    gray = cv2.cvtColor(image, cv2.COLOR_RGB2GRAY)

    # Blur the gray image
    gray = cv2.GaussianBlur(gray, (kernel_size, kernel_size), 0)

    # create binary image of the blurred image
    _, binary = cv2.threshold(
        gray, threshold_min, threshold_max, cv2.THRESH_BINARY)

    labels, label_store, label_contour_dict = find_fluidic_components_and_contours(
        image)

    circles = circle_finder(image)

    black_image = np.zeros(binary.shape, dtype=np.uint8)

    fluidic_part_structure = {}

    # Store the (x, y) position of every circle first in a list
    circles_pos = [(x, y) for x, y, r in circles[0, :]]
    # Change to list to a set, this is for easy comparison later
    circles_pos_set = set(circles_pos)

    # Retrieve all the unique label values from the array labels
    different_labels = np.unique(labels)
    # Remove the label 0 from the list, because that is just
    # the label for the background
    different_labels = different_labels[1:]

    for label in different_labels:

        fluidic_part_structure[label] = dict()

        # This part is to determine which circles are in the
        # component of a certain label

        # Retrieve the indices of labels where labels == label
        label_indices = (labels == label).nonzero()

        # Store all the (x, y) coordinates for a certain label in a set
        label_pos_set = set([(x, y)
                             for x, y in zip(label_indices[1], label_indices[0])])

        # Retrieve the circle positions that are in a certain component
        label_circles = list(label_pos_set.intersection(circles_pos_set))

        count = 1
        for i in range(len(label_circles)):
            for j in range(i, len(label_circles)):

                if i != j:
                    # Create a copy of the completely black image
                    cblack_image = np.copy(black_image)
                    # Create a copy of the binary image
                    cbinary = np.copy(binary)

                    circle_pos_1 = label_circles[i]
                    circle_pos_2 = label_circles[j]
                    cv2.line(cblack_image, circle_pos_1,
                             circle_pos_2, 255, thickness=1)

                    cv2.line(cbinary, circle_pos_1,
                             circle_pos_2, 0, thickness=1)

                    line_mask = (cblack_image == 255)

                    if np.all(binary[line_mask] == 255):
                        fluidic_part_structure[label][count] = (
                            cblack_image == 255).nonzero()
                        count += 1

                    # # show image
                    # cv2.imshow('Frame', cbinary)
                    # if cv2.waitKey(0) & 0xFF == ord('q'):  # press q to quit
                    #     cv2.destroyAllWindows()
    return fluidic_part_structure
