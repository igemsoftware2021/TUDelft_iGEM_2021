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
    _, bi = cv2.threshold(blur, threshold_min, threshold_max, cv2.THRESH_BINARY)
    
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
    
    #calculate the perimeter of the 4 detected squares and divide by 4 to obtain the width
    width_square = np.zeros(len(reccontour))
    for i in range(len(reccontour)):
        square = reccontour[i]
        peri_i = cv2.arcLength(square, True)
        width_square[i] = peri_i / 4

    # Calculate the average width and the conversion factor, multiply pixel distance by conversion facctor
    # to obtain the length in mm
    width_av = np.mean(width_square) #pixel
    conversion_factor = length_squares / width_av #mm/pixel

    return conversion_factor


def find_connected_components_contours(image, kernel_size=3, threshold_min=120, threshold_max=255, min_area=300):
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
        area = stats[i, cv2.CC_STAT_AREA]
        if area > min_area:
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

    return labels, label_store, label_contour_dict


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


def circle_finder(image, kernel_size=3, threshold_min=120, threshold_max=255, min_area=300, vector_range=10, max_angle=170):
    # turn RGB into gray image
    gray = cv2.cvtColor(image, cv2.COLOR_RGB2GRAY)

    # Blur the gray image
    gray = cv2.GaussianBlur(gray, (kernel_size, kernel_size), 0)

    # create binary image of the blurred image
    _, binary = cv2.threshold(
        gray, threshold_min, threshold_max, cv2.THRESH_BINARY)

    labels, label_store, label_contour_dict = find_connected_components_contours(
        image, kernel_size=kernel_size, threshold_min=threshold_min, threshold_max=threshold_max, min_area=min_area)

    # Dictionary where the key is the label number and the value is the array filled
    # with True and False
    label_switch_point_dict = {}

    # Now check where the circles are
    for label in label_store:
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

    # An image where all the points are black except for
    # the switch points, which will be white
    switch_point_image = np.zeros(
        (image.shape[0], image.shape[1]), dtype=np.uint8)

    for label in label_store:
        contour = label_contour_dict[label]
        switch_point = label_switch_point_dict[label]

        for i in range(contour.shape[0]):
            # switch_point contains only True and False, so if True
            if switch_point[i]:
                x, y = contour[i][0]
                switch_point_image[y, x] = 255

    circles = cv2.HoughCircles(switch_point_image, cv2.HOUGH_GRADIENT, 1, 20,
                               param1=100, param2=30, minRadius=15, maxRadius=500)

    return circles
