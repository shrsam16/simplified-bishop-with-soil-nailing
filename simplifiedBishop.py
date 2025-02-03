import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Arc
def directionAngle(x1: float,y1:float,x2:float,y2:float)->float:
    delx = x2-x1
    dely = y2-y1
    #edge cases
    if delx == 0:
        if dely==0:
            raise ValueError
        elif dely>0:
            return 90.0
        else:
            return 270.0
    if dely == 0:
        if delx>0:
            return 0.0
        else:
            return 180.0
    if delx>0:
        if dely>0:
            return np.rad2deg(np.arctan(abs(dely/delx)))
        else:
            return 360-np.rad2deg(np.arctan(abs(dely/delx)))
    if delx<0:
        if dely>0:
            return 180-np.rad2deg(np.arctan(abs(dely/delx)))
        else:
            return 180+np.rad2deg(np.arctan(abs(dely/delx)))

xAxisLimits = {'min': -50, 'max': 50}
yAxisLimits = {'min': -50, 'max': 50}

##soil parameters
UnitWeight = 18 #kN/m^3
Cohesion = 40 #kPa
FrictionAngle = 0 #degrees

slopeHeight = 7.0
slopeAnlge = 80
slopeBottom = (0,0)
# slopeTop = (1.234,7)
slopeTop = (slopeBottom[0]+np.cos(np.deg2rad(slopeAnlge))*slopeHeight, slopeBottom[1]+np.sin(np.deg2rad(slopeAnlge))*slopeHeight)
center = (-1.96, 11.09)
resolution = 0.001
noOfSlices = 90
sliceWidth = 1.3 #[1.28,1.22,1.35,1.41,1.56,2.88]

radius = np.sqrt((center[0]-slopeBottom[0])**2+(center[1]-slopeBottom[1])**2)

#to draw bottom horizontal
bottomX = np.arange(xAxisLimits['min'], slopeBottom[0], resolution)
bottomY = [slopeBottom[1] for each in bottomX]

#to draw top horizontal
topX = np.arange(slopeTop[0], xAxisLimits['max'], resolution)
topY = [slopeTop[1] for each in topX]

#to draw the slope
slopeEqn = np.poly1d(np.polyfit((slopeBottom[0], slopeTop[0]), (slopeBottom[1], slopeTop[1]), 1))
# ^ creates a function to get y along the slope for given x(point or list)
SlopeX = np.arange(slopeBottom[0], slopeTop[0], resolution)
SlopeY = slopeEqn(SlopeX)

#to draws the arc of slip surface

startAngle = directionAngle(center[0], center[1], slopeBottom[0], slopeBottom[1]) #360+np.rad2deg(np.arctan((center[1]-slopeBottom[1])/(center[0]-slopeBottom[0])))
endAngle = 360-np.rad2deg(np.arcsin((center[1]-slopeTop[1])/radius))
SlipSurface = Arc(center, radius*2, radius*2, theta1 = startAngle, theta2 = endAngle)

#points of interest for calculation
arcEndX = center[0] + np.cos(np.deg2rad(endAngle))*radius
if noOfSlices:
    sliceWidth = (arcEndX-slopeBottom[0])/noOfSlices
noOfSlicesBeforeTop = int((slopeTop[0]-slopeBottom[0])/sliceWidth)
sliceWidthBeforeTop = (slopeTop[0]-slopeBottom[0])/noOfSlicesBeforeTop
sliceWidthAfterTop = (arcEndX-slopeTop[0])/(noOfSlices-noOfSlicesBeforeTop)
sliceX = np.append(np.arange(slopeBottom[0], slopeTop[0], sliceWidthBeforeTop), np.arange(slopeTop[0], arcEndX, sliceWidth))
sliceX = np.append(sliceX, arcEndX)

## points along the slope and top 
sliceYTop = slopeEqn(sliceX)

for n, j in enumerate(sliceYTop):
    if j>slopeTop[1]:
        sliceYTop[n] = slopeTop[1]


## to find the points along the curve
midAngle = np.rad2deg(np.arcsin((sliceX-center[0])/radius))
sliceYBottom = center[1] - (radius * np.cos(np.deg2rad(midAngle)))

#plot each point of interest
for i, ju, jd in zip(sliceX, sliceYTop, sliceYBottom):
    plt.plot(i,ju,'go')
    plt.plot(i,jd,'rx')

    y = np.arange(jd,ju,resolution)
    x = [i for each in y]

    plt.plot(x,y,'g')




#plot all other stuff
plt.plot(bottomX, bottomY, color= 'black')
plt.plot(topX, topY, color= 'black')
plt.plot(SlopeX, SlopeY, 'black')
plt.plot(center[0], center[1], 'go')
plt.gca().add_patch(SlipSurface)


#calculation
## alpha
alphaX = np.diff(sliceX)
alphaY = np.diff(sliceYBottom)
alpha = alphaY/alphaX
alpha = np.rad2deg(np.arctan(alpha))



gamma = UnitWeight
c= Cohesion
phi = FrictionAngle

#base widths
b = alphaX

#mid-height
height = sliceYTop-sliceYBottom
heightcopy = np.append(np.copy(height[1:]), 0)
h = ((height + heightcopy)/2)[:-1]

#weight of slices and other calculations
w = gamma * b * h #col 4
sinAlpha = np.sin(np.deg2rad(alpha)) #col 5
wsinAlpha = w*sinAlpha #col 6

sumWsinAlpha = np.sum(wsinAlpha)
l = np.sqrt(alphaX**2 + alphaY**2)
cl =  c * l#col 7

F = 1
Fnew = 5
tolerance = 0.005

#calculation loop
while abs(F-Fnew)>=tolerance:
    F = Fnew
    col8 = w - (sinAlpha * cl)/F #col 8
    cosAlpha = np.cos(np.deg2rad(alpha)) #col 9
    tanPhiTanAlpha = np.tan(np.deg2rad(phi))*np.tan(np.deg2rad(alpha)) #col 10
    col11 = cosAlpha/(1+tanPhiTanAlpha/F) #col 11
    N = col8/col11 #col 12
    NtanPhi = N * np.tan(np.deg2rad(phi))  #col 13
    col14 = cl + NtanPhi  #col 14
    sumcol14 = np.sum(col14)

    Fnew = sumcol14/sumWsinAlpha
    print(f'{Fnew}\n')
plt.show()