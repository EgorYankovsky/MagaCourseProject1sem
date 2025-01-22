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

_bmax = -1.0 * float(sys.float_info.max)
_bmin = float(sys.float_info.max)

def readPoints(sourcePath : str):
    global _bmax, _bmin
    f = open(sourcePath, 'r')
    lines = f.readlines()
    lines.pop(0)
    for line in lines:
        data = line.split(" ")
        pointsArr.append(Point(float(data[0]), float(data[1]), float(data[2])))
        if float(data[0]) > _bmax: _bmax = float(data[0])
        if float(data[0]) < _bmin: _bmin = float(data[0])
        
        if float(data[1]) > _bmax: _bmax = float(data[1])
        if float(data[1]) < _bmin: _bmin = float(data[1])
        
        if float(data[2]) > _bmax: _bmax = float(data[2])
        if float(data[2]) < _bmin: _bmin = float(data[2])
        
        
def readRibs(sourcePath : str):
    f = open(sourcePath, 'r')
    lines = f.readlines()
    lines.pop(0)
    for line in lines:
        data = line.split(" ")
        ribsArr.append(Rib(int(data[0]), int(data[1])))
    
def drawScatter():
    global _bmax, _bmin
    fig = plt.figure(figsize=(19.80, 10.80))
    ax = fig.add_subplot(projection="3d")
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    ax.set_xlim(_bmin, _bmax)
    ax.set_ylim(_bmin, _bmax)
    ax.set_zlim(_bmin, _bmax)
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