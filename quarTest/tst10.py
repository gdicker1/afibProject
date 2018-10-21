import os
import numpy as np
from matplotlib import pyplot as plt
from CMAES_evolStrat import getQuarantined

if __name__ == "__main__":
    ablnFile = os.path.join(os.getcwd(), "quarTest", "abln_10.txt")
    print(ablnFile)
    if(os.path.exists(ablnFile)):
        cells = np.loadtxt(ablnFile, dtype=int, delimiter=",").reshape(10,10)

        plt.imshow(1-cells)
        plt.show()

        quarantined = getQuarantined(cells, 10)
        print("There are ", quarantined, "cells separated from the rest")

        plt.imshow(cells)
        plt.show()
