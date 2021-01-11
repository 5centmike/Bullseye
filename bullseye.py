#!python3
#lower = int(args.lower)

import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import sys
import argparse
import signal
import math
from scipy.stats import zscore
from bisect import bisect_left
parser = argparse.ArgumentParser(description='Creates normal heatmaps and creates bullseye plots which account for the manhattan distance between each site.')
parser.add_argument('-i', action = 'store', dest = 'input_file', required = True, type=str, help='input file')
parser.add_argument('-o', action = 'store', dest = 'output_file', required = True, type=str, help='prefix for output plots')
parser.add_argument('-c', action = 'store', dest = 'colorscheme', required = False, help='matplotlib_colors', default = "Red")
parser.add_argument('-z', dest='znorm', action='store_true', help='znorm each ring in bullseye plot', default = False)
parser.add_argument('-s', dest='squareBullseye', action='store_true', help='Trim edges to make a square bullseye plot', default = False)
parser.add_argument('-l', dest='lower', action='store', required = False, help='minvalue for color scale', type=float)
parser.add_argument('-u', dest='upper', action='store', required = False, help='maxvalue for color scale', type=float)
args = parser.parse_args()


spot = str(args.input_file)
outpre = str(args.output_file)
heatcol = str(args.colorscheme)

#lower = int(args.lower)
#upper = int(args.upper)
znorm = args.znorm
square = args.squareBullseye


def takeClosest(myList, myNumber):
    """
    Assumes myList is sorted. Returns closest value to myNumber.

    If two numbers are equally close, return the smallest number.
    """
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return 0
    if pos == len(myList):
        return -1
    before = myList[pos - 1]
    after = myList[pos]
    if after - myNumber < myNumber - before:
       return pos
    else:
       return pos - 1
mat = np.loadtxt(spot)

#this is the size of the starting matrix
matsize = len(mat)
if (matsize % 2) == 0:
    print("Bullseye requires a central bin, therefore the number of rows and columns in the input matrix should be odd")
    sys.exit(0)

#this is the number of rings in the plot which corresponds to the Manhattan distance
rings = math.floor(math.sqrt(((matsize//2)**2)/2))

#This is the maximum number of rings that will be needed if square is set to True.
maxrings = math.ceil(math.sqrt(rings**2 + rings**2))

#this is the number of radial slices the circle is divided in to. Too few and the segments will be distorted.
#slices = 2000
slices = (matsize-1)*100
#this toggles on z-score normalization of each ring.
#znorm = False

#this toggles between "closest-angle" mapping and "uniform segment length" mapping
uniform = True

#this makes a square plot by trimming the edges of the plot. 
#square = False

#First we compute the Manhattan distance and arctan of every point in the starting matrix
arctan = np.zeros((matsize,matsize))
manhattandistance = np.zeros((matsize,matsize))
for i in range(0,len(manhattandistance)):
    for j in range(0,len(manhattandistance)):
        manhattandistance[i,j] = abs(j-int((matsize-1)/2)) + abs(i-((matsize-1)/2))
        arctan[i,j] = np.arctan2(((j-int((matsize-1)/2))),((i-((matsize-1)/2))))

      
#We rearrange the Manhattan distances and arctans into lists based on the Manhattan distances.
#This means each ring is now represented by a list of scores and arctans.
manhat_arctans = []
manhat_scores = []
for a in range(matsize):
    manhat_arctans.append([])
    manhat_scores.append([])
    for b in range(matsize):
        for c in range(matsize):
            if (manhattandistance[b,c] == a):
                manhat_scores[a].append(mat[b,c])
                manhat_arctans[a].append(arctan[b,c])
                
#Next we construct the meshgrid for the circular plot.
rad = np.linspace(0, maxrings, maxrings+1)
azm = np.linspace(0, 2 * np.pi, slices)
r, th = np.meshgrid(rad, azm)

#now we trim the edges of the meshgrid to make a square
if square:
    for i in range(slices):
        for j in range(maxrings+1):
            #compute the cartesian coordinates of each quad to see if they need trimming
            nrad = r[i,j]
            theta = th[i,j]
            cosine = math.cos(theta)
            x = nrad * cosine
            sine = math.sin(theta)
            y = nrad * sine
            if abs(x) > abs(y):
                if abs(x) > rings:
                    nrad = abs(rings/cosine)
            else:
                if abs(y) > rings:
                    nrad = abs(rings/sine)
            r[i,j] = nrad

#This will be the new matrix that will be plotted. We need to fill it with the appropriate values.

C = np.zeros([slices,int(maxrings)+1])


dscoretmp1 = 0
dscoretmp2 = 0
dscore = 0
#For each ring we convert the arctans into thetas by making them all positive
#and adjusting for the slices.
for ring in range(int(maxrings)):
    arctans = manhat_arctans[ring]
    thetas = [(x + np.pi) * ((slices/2)/np.pi) for x in arctans]
    scores = manhat_scores[ring]
    #Z-score normalization toggled by boolean value above.
    if znorm:
        if ring:
            scores = zscore(scores).tolist()
        else:
            scores = [0]

    #In order for the values to appropriately wrap around to 0 we need to add the last score to the beginning.
    thetas.append(0)
    scores.append(scores[0])
    #Here we sort the scores and thetas by the theta in order to help finding the closest one.
    sorted_scores = [x for _,x in sorted(zip(thetas,scores))]
    if len(sorted_scores) > 1:
        quatrang = int(float(len(sorted_scores)-1)/4)
        quatstart = quatrang+1
        quatend = quatrang+quatrang
        for hloc in sorted_scores:
            
            dscoretmp1+=1

        for hloc in sorted_scores[quatstart:quatend]:
            if hloc >= 1:
                dscoretmp2+=1 
    sorted_thetas = sorted(thetas)
    #Here to produce a uniform map we map the actual angles to idealized angles.
    if uniform:
        if ring:
            slices_per_segment = slices/(len(scores)-1)
            for t in range(len(sorted_thetas)):
                sorted_thetas[t] = slices_per_segment * t
                
    
    #Now for each slice in the ring we find the score corresponding to the nearest theta
    #and assign it to the appropriate coordinate in our matrix to be plotted: C.
    for i in range(slices):
        idx = takeClosest(sorted_thetas,i)
        C[i,ring] = sorted_scores[idx]
        
#create a juicebox-like red to white colormap
if heatcol == "Red":
    br = LinearSegmentedColormap.from_list("bright_red",[(1,1,1),(1,0,0)])
else:
    br = heatcol
        
ax = plt.subplot(projection="polar")

if znorm:
    plt.pcolormesh(th,r,C,cmap="seismic",vmin=-4,vmax=4)
    if (dscoretmp1 > 0):
        dscore=float(dscoretmp2/dscoretmp1)*100
    print("ADA score is %s\n" % (str(round(dscore,2))))
else: 
    if args.lower is not None and args.upper is not None:
        plt.pcolormesh(th,r,C,cmap=br,vmin=args.lower,vmax=args.upper)
    elif args.lower is not None:
        plt.pcolormesh(th,r,C,cmap=br,vmin=args.lower)
    elif  args.upper is not None:
        plt.pcolormesh(th,r,C,cmap=br,vmax=args.upper)
    else:
        plt.pcolormesh(th,r,C,cmap=br)
if square:
    ax.axis('off')
else:
    xlabels=["","","Left Anchor","","Right Anchor","","",""]
    ax.set_xticklabels(xlabels, rotation=40, ha="right")
    ax.set_yticklabels([])
    plt.colorbar()
    
#Rotate the plot so that the anchors go back to their original positions.
ax.set_theta_offset(np.pi/2)

plt.savefig(outpre + '_bullseye.png')

plt.clf()
center = (matsize-1)//2

left = center-rings

right = center+rings+1

newmat = (mat[int(left):int(right),int(left):int(right)])
if args.lower is not None and args.upper is not None:
    plt.imshow(newmat,cmap=br,vmin=args.lower,vmax=args.upper)
elif args.lower is not None:
    plt.imshow(newmat,cmap=br,vmin=args.lower)
elif  args.upper is not None:
    plt.imshow(newmat,cmap=br,vmax=args.upper)
else:
    plt.imshow(newmat,cmap=br)

plt.colorbar()
plt.savefig(outpre + '_normal.png')
