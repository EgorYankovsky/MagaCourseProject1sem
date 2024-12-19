import matplotlib.pyplot as plt
import sys
import os

class Point():
    def __init__(self, x : float, y : float, z:float):
        self.x = x
        self.y = y
        self.z = z

class Rib():
    def __init__(self, r1 : int, r2 : int):
        self.r1 = r1
        self.r2 = r2

pointsArr = []
ribsArr = []

def readPoints(sourcePath : str):
    f = open(sourcePath, 'r')
    lines = f.readlines()
    lines.pop(0)
    for line in lines:
        data = line.split(" ")
        pointsArr.append(Point(float(data[0]), float(data[1]), float(data[2])))
        
def readRibs(sourcePath : str):
    f = open(sourcePath, 'r')
    lines = f.readlines()
    lines.pop(0)
    for line in lines:
        data = line.split(" ")
        ribsArr.append(Rib(int(data[0]), int(data[1])))
    
def drawScatter():
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    for point in pointsArr:
        ax.scatter(point.x, point.y, point.z, color="black")
    for rib in ribsArr:
        x = [pointsArr[rib.r1].x, pointsArr[rib.r2].x]
        y = [pointsArr[rib.r1].y, pointsArr[rib.r2].y] 
        z = [pointsArr[rib.r1].z, pointsArr[rib.r2].z]
        ax.plot(x, y, z, color='black')
    plt.show()

if __name__ == "__main__":
    args = sys.argv
    readPoints(os.path.abspath(args[1]))
    readRibs(os.path.abspath(args[2]))
    drawScatter()