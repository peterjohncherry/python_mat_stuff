import matplotlib.pyplot as plt
import numpy as np

class Matrix_Images:

    name_ = "unnamed"
    mat_imgs_ = []

    def __init__(self, name):
        self.name_ = name

    def add_plot(self, mat_in):
        self.mat_imgs_.append(self.plot_matrix(mat_in))

    def save_imgs(self):
        ii = 1
        for mat_img in self.mat_imgs_:
            fig = mat_img.get_figure()
            fig.savefig(self.name_+'_'+str(ii)+'.png')
            ii+=1

    @staticmethod
    def plot_matrix(mat_in):
        return plt.matshow(mat_in, cmap='gray')

