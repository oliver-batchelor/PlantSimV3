import cv2
import pathlib
import random
import numpy as np


def CropPadCenter(np_img, img_shape=(600, 600)):
    y,x,c = np_img.shape
    if y < img_shape[1]:
        np_img = np.pad(np_img, [(0, img_shape[1]-y), (0, 0), (0, 0)], mode='constant')
    if x < img_shape[0]:
        np_img = np.pad(np_img, [(0, 0), (0, img_shape[0]-x), (0, 0)], mode='constant')
    y, x, c = np_img.shape
    startx = x//2-(img_shape[0]//2)
    starty = y//2-(img_shape[1]//2)
    return np_img[starty:starty+img_shape[1], startx:startx+img_shape[0]]


def pyrDownSmple(file_path):

    ann_image = cv2.imread(file_path)
    ann_greyscale = cv2.cvtColor(ann_image, cv2.COLOR_BGR2GRAY)

    ann_dwnsmpl = cv2.pyrDown(ann_greyscale, (300, 300)) * 10
    ann_dwnsmpl = cv2.pyrDown(ann_dwnsmpl, (150, 150)) * 10
    ann_dwnsmpl = cv2.pyrDown(ann_dwnsmpl, (75, 75))

    ann_blured = cv2.blur(ann_dwnsmpl, (3, 3))
    cv2.imwrite(file_path, ann_blured)


def AddRandomBackground(frgnd_img_path, bckgnd_img_dir):
    """Adds a random image from the image directory to the black background"""

    background_img_paths = []
    [background_img_paths.append(path) for path in pathlib.Path(bckgnd_img_dir).iterdir()]

    frgnd_image = cv2.imread(frgnd_img_path)
    ret, frgnd_mask = cv2.threshold(cv2.cvtColor(cv2.GaussianBlur(frgnd_image, (3, 3), 1),
                                                 cv2.COLOR_BGR2GRAY),10,255,cv2.THRESH_BINARY_INV)

    bckgnd_img = cv2.imread(str(random.choice(background_img_paths)))
    bckgnd_img = CropPadCenter(bckgnd_img)
    bckgnd_masked = cv2.bitwise_and(bckgnd_img, bckgnd_img, mask=frgnd_mask)

    frgnd_bckgnd = frgnd_image + bckgnd_masked
    cv2.imwrite(frgnd_img_path, cv2.GaussianBlur(frgnd_bckgnd, (3, 3), 1))