import os
import numpy as np
from matplotlib import cm, pyplot as plt
from CMAES_evolStrat import getQuarantined

if __name__ == "__main__":
    quarDir = os.path.join(os.getcwd(), "quarTest")
    ablnFiles = []
    for x in os.listdir(quarDir):
        if "abln" in x:
            ablnFiles.append(os.path.join(quarDir, x))

    current_cmap = cm.get_cmap('Greys')
    current_cmap = current_cmap.set_bad("red")

    for x in ablnFiles:
        if "10" in x:
            width = 10
        else:
            width = 80
        cells = np.loadtxt(x, dtype=float).reshape(width, width)
        cpcells = np.copy(cells)
        for i in range(width):
            for j in range(width):
                if cpcells[i, j] == 1:
                    cpcells[i, j] = np.nan
        plt.imshow(cpcells, cmap="Greys")
        plt.show()

        getQuarantined(cells, width)
        cpcells = np.copy(cells)
        for i in range(width):
            for j in range(width):
                if cpcells[i, j] == 1:
                    cpcells[i, j] = np.nan
        plt.imshow(cpcells, cmap="Greys")
        plt.show()
