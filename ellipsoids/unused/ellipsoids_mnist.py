# mnist stuff
import numpy as np


def imageToCoordinates(image: np.array, scaling = 10):
    [height,width] = image.shape
    positiveValues = sum(sum(image>0))
    coordinates = np.zeros([positiveValues,2], dtype=float)

    k = 0
    for i in np.arange(width):
        for j in np.arange(height):
            if image[i,j] > 0:
                coordinates[k,:] = [scaling*j,scaling*i]
                k = k+1

    return np.flip(coordinates,1)


def import_mnist():
    from keras.datasets import mnist
    (train_X, train_y), (test_X, test_y) = mnist.load_data()
    coordinates = imageToCoordinates(train_X[1])
    return coordinates