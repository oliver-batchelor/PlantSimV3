import cv2

def pyrDownSmple(file_path):

    ann_image = cv2.imread(file_path)
    ann_greyscale = cv2.cvtColor(ann_image, cv2.COLOR_BGR2GRAY)

    ann_dwnsmpl = cv2.pyrDown(ann_greyscale, (300, 300)) * 10
    ann_dwnsmpl = cv2.pyrDown(ann_dwnsmpl, (150, 150)) * 10
    ann_dwnsmpl = cv2.pyrDown(ann_dwnsmpl, (75, 75))

    ann_blured = cv2.blur(ann_dwnsmpl, (3, 3))
    cv2.imwrite(file_path, ann_blured)