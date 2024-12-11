import matplotlib.pyplot as plt
import sys
import os

class Point():

    def __init__(self, x : float, y : float, z:float):
        self.x = x
        self.y = y
        self.z = z

pointsArr = []

def readPoints(sourcePath : str):
    f = open(sourcePath, 'r')
    lines = f.readlines()
    for line in lines:
        data = line.split(" ")
        pointsArr.append(Point(float(data[0]), float(data[1]), float(data[2])))
        
def drawScatter():
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    for point in pointsArr:
        ax.scatter(point.x, point.y, point.z, color="black")
        ax.set_xlabel('X Label')
        ax.set_ylabel('Y Label')
        ax.set_zlabel('Z Label')
    plt.show()

if __name__ == "__main__":
    args = sys.argv
    #print (type(os.path.abspath(args[1])))
    readPoints(os.path.abspath(args[1]))
    drawScatter()