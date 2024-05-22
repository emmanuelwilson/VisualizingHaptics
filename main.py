# <Visualizing Haptics, a tool to help visualize haptic environments and interactions>
# Copyright (C) <2024>  <Emmanuel Wilson>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

from HaplyHAPI import Board, Device, Mechanisms, Pantograph
import time
import serial.tools.list_ports
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d import axes3d
from mayavi import mlab
import numpy as np
import skimage.measure
import math
import time
#from sklearn.preprocessing import normalize


CW = 0
CCW = 1
pixelPerMeter = 4000
width = 600                                                                     #Width of the workspace
height = 600                                                                    #Height of the workspace
quiverW = 15                                                                    #Density of the quiver plot width
quiverH = 15                                                                    #Density of the quiver plot height
haptics = True
ind = np.indices((height,width))
indQ = np.indices((quiverW,quiverH))
indLinear = np.zeros((quiverW * quiverH,2))
indLinear[:,0] = indQ[0].reshape((quiverW * quiverH))
indLinear[:,1] = indQ[1].reshape((quiverW * quiverH))
springForces = np.zeros( (height, width, 2) )
textureForces = np.zeros( (height, width,2) )
frictionForces = np.zeros( (height, width,2) )
forceAmp = np.zeros( (height, width) )
forceDir = np.zeros( (height*width,2) )
lastPos = [0,0]

hardware_version = 3  # 2 for the metallic plate device, 3 for the newer plastic device

haplyBoard = Board
device = Device
SimpleActuatorMech = Mechanisms
pantograph = Mechanisms

com_ports = list(serial.tools.list_ports.comports())
# if only one port is available, use that one, otherwise ask the user to select
# the correct one
if len(com_ports) == 1:
    port = com_ports[0].device
else:
    # try:
    #     print("Select the COM port for the Haply board, if no Haply board press ENTER:")
    #     for i, port in enumerate(com_ports):
    #         print(str(i) + ": " + port.device)
    #     port = com_ports[int(input())].device
    # except:
        print("Invalid port selection, entering NO DEVICE authoring mode")
        haptics = False
        port = 0


# circle spring
# INPUT: x centroid, y centroid, r radius, y radius, spring constant, tilt angle
# OUTPUT: width by height by 2 matrix containing the X and Y force direction vector of the circle
def circleSpring(Xc, Yc, rx, ry, k, tilt=0):
    calculateR = ((((ind[1]-Xc)*np.cos(tilt)+(ind[0]-Yc)*np.sin(tilt))**2)/rx**2 +
              ((-(ind[1]-Xc)*np.sin(tilt)+(ind[0]-Yc)*np.cos(tilt))**2)/ry**2)
    # Calculate the distance in X from the centroid to any point of the shape
    circleForceX = np.where(calculateR<=1,
             ind[1]-Xc,0).reshape((width*height,1))

    # Calculate the distance in Y from the centroid to any point of the shape
    circleForceY = np.where(calculateR<=1,
             ind[0]-Yc,0).reshape((width*height,1))

    # Calculate the amplitude 0-1 * abs(K)
    circleForceAmp = np.where(calculateR <=1,
             (1-calculateR)*np.abs(k),0).reshape((height,width))

    # Calculate the amplitude 0-1 * K USED FOR 3D VISUALIZATION
    circleForceAmpSign = np.where(calculateR <=1,
             (1-calculateR)*k,0).reshape((height,width))

    # circleForceAmp[Yc,Xc] = 1
    if (k < 0):
        circleForceX = -circleForceX
        circleForceY = -circleForceY

    # Normalize direction vectors
    circleForceDir = (np.concatenate((circleForceX,circleForceY),axis=1))
    div =  np.sqrt(np.square(circleForceDir).sum(axis=1))
    div = np.where(div==0,1,div).reshape(width*height,1)

    circleForceDir = circleForceDir/div

    return circleForceAmp, circleForceDir, circleForceAmpSign

# Set the friction of a circle to a value between 0 and 1, it will apply a dampening force proportional to the velocity of the end effector
# INPUT: X centroid, Y centroid, width radius, height radius, friction coefficient [0-1], ellipse tilt angle
def circleFriction(Xc,Yc,rx,ry, friction=0.0,tilt=0):
    calculateR = ((((ind[1]-Xc)*np.cos(tilt)+(ind[0]-Yc)*np.sin(tilt))**2)/rx**2 +
              ((-(ind[1]-Xc)*np.sin(tilt)+(ind[0]-Yc)*np.cos(tilt))**2)/ry**2)
    # Set all points within the circle to friction
    frictionMap = np.where(calculateR<=1,
             friction,0).reshape((width*height,1))

    return frictionMap

# Set texture as a series of small springs within the designated area
# INPUT: X centroid, Y centroid, width radius, height radius, texture max amplitude, texture density = approximate number of springs present in the space,
#        radius size of texture springs, friction coefficient [0-1], ellipse tilt angle
# OUTPUT: Amplitude and direction matrices of size width and height
def circleTexture(Xc,Yc,rx,ry, textK, textureDensity, texutreRadius,tilt=0):
    calculateR = ((((ind[1]-Xc)*np.cos(tilt)+(ind[0]-Yc)*np.sin(tilt))**2)/rx**2 +
              ((-(ind[1]-Xc)*np.sin(tilt)+(ind[0]-Yc)*np.cos(tilt))**2)/ry**2)
    xrand = np.random.randint(Xc-rx,Xc+rx,textureDensity)
    yrand = np.random.randint(Yc-ry,Yc+ry,textureDensity)

    textAmp = np.zeros((width,height,1))
    textDir = np.zeros((width*height,2))
    for i in range(textureDensity):
      if (calculateR[yrand[i],xrand[i]] <= 1):
        texAmptemp, textDirtemp = circleSpring(xrand[i],yrand[i],texutreRadius,texutreRadius,textK)
        textAmp, textDir = addForces(textAmp,textDir,texAmptemp,textDirtemp)


    return textAmp, textDir

# Set texture as a series of small springs within the designated area
# INPUT: X centroid, Y centroid, width radius, height radius, texture max amplitude, texture density = approximate number of springs present in the space,
#        radius size of texture springs, friction coefficient [0-1], ellipse tilt angle
# OUTPUT: Amplitude and direction matrices of size width and height
def circleTextureLattice(Xc,Yc,rx,ry, textK, textureDensity, texutreRadius,tilt=0):
    calculateR = ((((ind[1]-Xc)*np.cos(tilt)+(ind[0]-Yc)*np.sin(tilt))**2)/rx**2 +
              ((-(ind[1]-Xc)*np.sin(tilt)+(ind[0]-Yc)*np.cos(tilt))**2)/ry**2)
    xlat = np.arange(Xc-rx,Xc+rx,textureDensity)
    ylat = np.arange(Yc-ry,Yc+ry,textureDensity)

    textAmp = np.zeros((width,height,1))
    textDir = np.zeros((width*height,2))
    for xl in range(len(xlat)):
      for yl in range(len(ylat)):
        if (calculateR[ylat[yl],xlat[xl]] <= 1):
          texAmptemp, textDirtemp, blah = circleSpring(xlat[xl],ylat[yl],texutreRadius,texutreRadius,textK)
          textAmp, textDir = addForces(textAmp,textDir,texAmptemp,textDirtemp)


    return textAmp, textDir

# Set texture as a series of small springs within the designated area
# INPUT: X centroid, Y centroid, width radius, height radius, texture max amplitude [0-1], texture density = approximate number of events that occur in a second, ellipse tilt angle
# OUTPUT: Centroid location of each scatter point and the alpha value based on the max output
def circleTextureScatter(Xc,Yc,rx,ry, textMax, textureDensity,tilt=0):
    calculateR = ((((ind[1]-Xc)*np.cos(tilt)+(ind[0]-Yc)*np.sin(tilt))**2)/rx**2 +
              ((-(ind[1]-Xc)*np.sin(tilt)+(ind[0]-Yc)*np.cos(tilt))**2)/ry**2)
    xlat = np.arange(Xc-rx,Xc+rx,int(120/textureDensity))
    ylat = np.arange(Yc-ry,Yc+ry,int(120/textureDensity))
    latgood = np.zeros((len(xlat)*len(ylat),2))
    latgoodind = np.array([])
    countind = 0
    count = 0
    for xl in range(len(xlat)):
      for yl in range(len(ylat)):
        if (calculateR[ylat[yl],xlat[xl]] <= 1):
          latgoodind = np.append(latgoodind,countind).astype(int)
          latgood[countind,:] = [xlat[xl],ylat[yl]]
          count += 1
        countind += 1
    print(np.shape(latgoodind))
    latgood = latgood[latgoodind,:]

    return latgood, textMax

# This function will take circle parameters and output the location of its circumference along with the friction value at every point
# INPUT: friction map, X centroid, Y centroid, width radius, height radius, tilt of ellipse
# OUTPUT: Circumference location of the circle and the friction value at each point
def findCircleCircumferenceFriction(friction, Xc,Yc,rx,ry, tilt = 0):
    calculateR = ((((ind[1]-Xc)*np.cos(tilt)+(ind[0]-Yc)*np.sin(tilt))**2)/rx**2 +
              ((-(ind[1]-Xc)*np.sin(tilt)+(ind[0]-Yc)*np.cos(tilt))**2)/ry**2)

    circleInd = np.where((calculateR<=1) & (calculateR>=0.99))
    frictionMap = friction[circleInd]

    return circleInd, frictionMap

# This function will take circle parameters and output the location of its circumference along with the friction value at every point in circular order
# INPUT: friction map, X centroid, Y centroid, width radius, height radius, tilt of ellipse
# OUTPUT: Circumference location of the circle and the friction value at each point
def findCircleCircumferenceFrictionOrdered(friction, Xc,Yc,rx,ry, tilt = 0):
    calculateR = ((((ind[1]-Xc)*np.cos(tilt)+(ind[0]-Yc)*np.sin(tilt))**2)/rx**2 +
              ((-(ind[1]-Xc)*np.sin(tilt)+(ind[0]-Yc)*np.cos(tilt))**2)/ry**2)
    xInd = np.array([])
    yInd = np.array([])

    longestdist = int(math.hypot(height,width))
    circleInd = np.where((calculateR<=1) & (calculateR>=0.99))
    if rx > ry:
      rdist = np.arange(0,rx+5,1)
    else:
      rdist = np.arange(0,ry+5,1)

    for i in np.arange(0,360,1):
      x = np.round(rdist*np.cos(np.deg2rad(i)) + Xc).astype(int).reshape(len(rdist),1)
      y = np.round(rdist*np.sin(np.deg2rad(i)) + Yc).astype(int).reshape(len(rdist),1)
      for j in range(len(rdist)):
        indLineX = np.where((circleInd[1] == x[j]))
        indLineY = np.where((circleInd[0] == y[j]))
        indLine = np.unique(np.intersect1d(indLineX,indLineY))
        if len(indLine) > 0:
          xInd = np.append(xInd,circleInd[1][indLine]).astype(int)
          yInd = np.append(yInd,circleInd[0][indLine]).astype(int)
          break

    circleIndOrdered = np.concatenate((yInd.reshape(len(yInd),1),xInd.reshape(len(xInd),1)),axis=1)

    frictionMap = friction[circleIndOrdered[:,0],circleIndOrdered[:,1]]
    return circleIndOrdered, frictionMap

# Calculate the force direction and amplitude based on the end effector position and heading
# INPUT: haplyBoard, device, force amplitude, force direction
# OUTPUT: Will apply forces onto the 2DIY 
#-----------------NOT IN WORKING ORDER YET----------------
def applyForces(haplyBoard,device, forceAmp,forceDir, frictionForces):
    if (haplyBoard.data_available()):
            device.device_read_data()
            motorAngle = device.get_device_angles()
            device_position = device.get_device_position(motorAngle)*pixelPerMeter
            # Apply Spring Forces
            forces = forceDir[device_position]*forceAmp[device_position]

            # Apply friction forces
            # forces += frictionForces[device_position]

            device.set_device_torques(forces)
            device.device_write_torques()
            print(forces)

# Calculate the resultant overlapping forces using force amplitude and direction
# INPUT: force amplitude object1 , force direction object1, force amplitude object1, force direction object 2 MUST BE OF SAME SIZE
# OUTPUT: Force amplitude and force direction of the same size as inputs
def addForces(forceAmp1,forceDir1,forceAmp2,forceDir2):
    forceAmp1 = forceAmp1.reshape(width*height)
    forceAmp2 = forceAmp2.reshape(width*height)

    forceAmpOut = forceAmp1 + forceAmp2
    forceDirOut = forceDir1 + forceDir2

    overlapIndex = np.where( (forceAmp1 != 0) & (forceAmp2 != 0))
    overlapIndex = np.transpose(np.array(overlapIndex))

    if (np.size(overlapIndex) > 0):

        forceOverlapX = forceDir1[overlapIndex,1]*forceAmp1[overlapIndex] + forceDir2[overlapIndex,1]*forceAmp2[overlapIndex]
        forceOverlapY = forceDir1[overlapIndex,0]*forceAmp1[overlapIndex] + forceDir2[overlapIndex,0]*forceAmp2[overlapIndex]

        overlapDir = (np.concatenate((forceOverlapY,forceOverlapX),axis=1)).reshape(len(overlapIndex),2)
        mag =  np.sqrt(np.square(overlapDir).sum(axis=1)).reshape(len(overlapIndex),1)
        forceAmpOut[overlapIndex] = mag
        forceDirOut[overlapIndex.reshape(len(overlapIndex))] = (overlapDir/mag)
        # plt.show()

    return forceAmpOut.reshape(height,width), forceDirOut

# Calculate the resultant overlapping forces using force amplitude and direction
# INPUT: force amplitude object1 , force direction object1, force amplitude object1, force direction object 2 MUST BE OF SAME SIZE
# OUTPUT: Force amplitude, positive values for pushing forces negative values for pulling forces
def addForcesSign(forceAmp1,forceDir1,forceAmp2,forceDir2):
    forceAmp1 = forceAmp1.reshape(width*height)
    forceAmp2 = forceAmp2.reshape(width*height)

    forceAmpOut = forceAmp1 + forceAmp2
    forceDirOut = forceDir1 + forceDir2

    overlapIndex = np.where( (forceAmp1 != 0) & (forceAmp2 != 0))
    overlapIndex = np.transpose(np.array(overlapIndex))

    if (np.size(overlapIndex) > 0):

        forceOverlapX = forceDir1[overlapIndex,1]*np.abs(forceAmp1[overlapIndex]) + forceDir2[overlapIndex,1]*np.abs(forceAmp2[overlapIndex])
        forceOverlapY = forceDir1[overlapIndex,0]*np.abs(forceAmp1[overlapIndex]) + forceDir2[overlapIndex,0]*np.abs(forceAmp2[overlapIndex])

        overlapDir = (np.concatenate((forceOverlapY,forceOverlapX),axis=1)).reshape(len(overlapIndex),2)
        mag =  np.sqrt(np.square(overlapDir).sum(axis=1)).reshape(len(overlapIndex),1)
        anglediff1 = np.arccos(np.transpose((np.transpose(forceDir1[overlapIndex,0]) * overlapDir[:,0] + np.transpose(forceDir1[overlapIndex,1]) * overlapDir[:,1])) / (mag))
        anglediff2 = np.arccos(np.transpose((np.transpose(forceDir2[overlapIndex,0]) * overlapDir[:,0] + np.transpose(forceDir2[overlapIndex,1]) * overlapDir[:,1])) / (mag))

        # overlapConstructive = np.where((anglediff1 < (np.pi/2)) | (anglediff1 > (np.pi - np.pi/2)), 1, 0)
        overlapConstructiveInd1 = np.where(((anglediff1 < (np.pi/2)) | (anglediff1 > (np.pi - np.pi/2))) & ((forceAmp1[overlapIndex]<0)))
        overlapConstructiveInd2 = np.where(((anglediff1 < (np.pi/2)) | (anglediff1 > (np.pi - np.pi/2))) & ((forceAmp2[overlapIndex]<0)))
        overlapConstructiveInd = np.unique(np.concatenate((overlapIndex[overlapConstructiveInd1],overlapIndex[overlapConstructiveInd2]),axis=0))
        overlapConstructive = np.zeros(np.shape(forceAmp1))
        overlapConstructive[overlapConstructiveInd] = 1


        forceAmpOut[overlapIndex] = mag
        overlapPushPull = (((forceAmp1 < 0) | (forceAmp2 < 0)) & ((forceAmp1 > 0) | (forceAmp2 > 0)) | ((forceAmp1 < 0) & (forceAmp2 < 0)))
        # overlapConstructive = (np.abs(forceAmpOut) > np.abs(forceAmp1)) & (np.abs(forceAmpOut) > np.abs(forceAmp2))
        test = np.where(((overlapPushPull == 1) & (overlapConstructive==1)), 2, 0)
        test = np.where(((overlapPushPull == 1) & (overlapConstructive==0)), 1, test)
        test = np.where(((overlapPushPull == 0) & (overlapConstructive==1)), 1, test)

        indoverlapConstructivePull = (overlapConstructive.astype(int) & overlapPushPull.astype(int))


        forceAmpOutBefore = forceAmpOut.copy()
        forceAmpOut = np.where(indoverlapConstructivePull == 1,-np.abs(forceAmpOut), forceAmpOut)

    return forceAmpOut.reshape(height,width)

# Adds friction maps together and doesn't allow values to go above 1
# INPUT: friction map 1, friction map2 MUST BE SAME SIZE
# OUTPUT: summed friction map with a cap of 1 of same dims as the inputs
def addFriction(friction1,friction2):
    frictionOut = np.where( (friction1+friction2) > 1, 1 , friction1+friction2)

    return frictionOut

# Fetch the proper color value for the friction visualization
# INPUT: width height of the workspace and Friction map at that pixel
# OUTPUT: friction color values of the same size as X
def cMap(ind,friction):
    c = plt.cm.autumn_r((np.clip(friction[ind[1], ind[0]], 0, 1)))

    return c

def main():
    print("Starting the application!")

    # Set your objects here:
    circleAmp1,circleDir1, circleAmp1Sign  = circleSpring(30,30,20,20,2)

    circleAmp2,circleDir2, circleAmp2Sign = circleSpring(300,300,100,100,-2)

    ellipseAmp,ellipseDir, ellipseAmpSign = circleSpring(400,400,100,200,2,np.deg2rad(-20))

    ellipseAmp2,ellipseDir2, ellipseAmpSign2 = circleSpring(200,200,100,200,-2,np.deg2rad(20))

    # Set Texture forces
    textCircleCentroids, textCircleAlpha = circleTextureScatter(100,250,60,100,0.5,5)
    textCircleCentroids2, textCircleAlpha2 = circleTextureScatter(500,150,80,50,0.9,3,np.deg2rad(20))

    forceAmp,forceDir = addForces(circleAmp1,circleDir1,circleAmp2,circleDir2)

    #Add 3D signed forces
    forceAmpSign = addForcesSign(circleAmp1Sign,circleDir1,circleAmp2Sign,circleDir2)
    forceAmpSign = addForcesSign(forceAmpSign,forceDir,ellipseAmpSign,ellipseDir)
    forceAmpSign = addForcesSign(forceAmpSign,forceDir,ellipseAmpSign2,ellipseDir2)

    textAmp = forceAmpSign[textCircleCentroids[:, 0].astype(int), textCircleCentroids[:, 1].astype(int)]
    textAmp2 = forceAmpSign[textCircleCentroids2[:, 0].astype(int), textCircleCentroids2[:, 1].astype(int)]

    # Add spring forces
    # forceAmp,forceDir = addForces(circleAmp1,circleDir1,circleAmp2,circleDir2)
    forceAmp,forceDir = addForces(forceAmp,forceDir,ellipseAmp,ellipseDir)
    forceAmp,forceDir = addForces(forceAmp,forceDir,ellipseAmp2,ellipseDir2)
    forceSpring = forceAmp.copy()

    # Set friction objects here
    circleFriction1 = circleFriction(300,300,100,100,0.6)
    circleFriction2 = circleFriction(400,400,100,200,0.4,np.deg2rad(-20))
    circleFriction3 = circleFriction(200,200,100,200,0.1,np.deg2rad(20))

    # Add friction coeffs
    frictionForces = addFriction(circleFriction1,circleFriction2)
    frictionForces = addFriction(frictionForces,circleFriction3).reshape(height,width)

    circleCircum1, frictionCircum1 = findCircleCircumferenceFrictionOrdered(frictionForces,300,300,100,100)
    circleCircum2, frictionCircum2 = findCircleCircumferenceFrictionOrdered(frictionForces,400,400,100,200,np.deg2rad(-20))
    ellipseCircum, frictionCircum3 = findCircleCircumferenceFrictionOrdered(frictionForces,200,200,100,200,np.deg2rad(20))

    circum = np.concatenate((circleCircum1, circleCircum2, ellipseCircum), axis=0)
    frictionCircum = np.concatenate((frictionCircum1, frictionCircum2, frictionCircum3), axis=0)
    amplitudeCircum = forceAmpSign[circum[:, 0], circum[:, 1]]
    # frictionColor = frictionCircum / np.max(frictionCircum)
    frictionColor = cMap(ind,frictionForces)

    # Condensed quiver plot
    quiverDirX = skimage.measure.block_reduce(forceDir[:,1].reshape(height,width), block_size=(int((height)/(quiverH)),int((width)/(quiverW))), func = np.mean).reshape(quiverW*quiverH,1)
    quiverDirY = skimage.measure.block_reduce(forceDir[:,0].reshape(height,width), block_size=(int((height)/(quiverH)),int((width)/(quiverW))), func = np.mean).reshape(quiverW*quiverH,1)

    quiverDir = (np.concatenate((quiverDirX,quiverDirY),axis=1))

    indQ1 = skimage.measure.block_reduce(ind[1].reshape(height,width), block_size=(int((height)/(quiverH)),int((width)/(quiverW))), func = np.median)#.reshape(quiverW*quiverH,1)
    indQ2 = skimage.measure.block_reduce(ind[0].reshape(height,width), block_size=(int((height)/(quiverH)),int((width)/(quiverW))), func = np.median)#.reshape(quiverW*quiverH,1)

    fig2, ax2 = plt.subplots()
    ax2.set_title("2D Representation")
    Q = ax2.quiver(indQ1,indQ2,quiverDir[:,1],-quiverDir[:,0],pivot = 'mid', width = 0.01,units ='inches', headwidth = 5)
    C = ax2.contour(ind[1],ind[0],frictionForces,cmap=cm.cool)
    V = ax2.imshow(forceAmp)
    S = ax2.scatter(textCircleCentroids[:,1],textCircleCentroids[:,0],alpha=textCircleAlpha,marker='X',color='orange')
    S2 = ax2.scatter(textCircleCentroids2[:,1],textCircleCentroids2[:,0],alpha=textCircleAlpha2,marker='X',color='orange')
    fig2.colorbar(V, ax=ax2)
    fig2.colorbar(C, ax=ax2)

    ## 3D Surface plots ##

    fig = mlab.figure(10)
    mlab.mesh(ind[1], ind[0], forceAmpSign[ind[1],ind[0]]*50, scalars=frictionForces[ind[1],ind[0]], colormap="cool")
    # mlab.surf(forceAmpSign, warp_scale='auto', colormap='black-white')
    mlab.points3d(textCircleCentroids2[:, 0], textCircleCentroids2[:, 1], textAmp2 * 50 + 10, colormap='Oranges',opacity=textCircleAlpha2)
    mlab.points3d(textCircleCentroids[:, 0], textCircleCentroids[:, 1], textAmp * 50 + 2, colormap='Oranges',opacity=textCircleAlpha)

    plt.show()
    mlab.show()


    if haptics:
        haplyBoard = Board("test", port, 0)
        device = Device(5, haplyBoard)
        pantograph = Pantograph(hardware_version)
        device.set_mechanism(pantograph)

        #device.add_encoder(1, CCW, 240, 10978, 2)
        #device.add_encoder(2, CCW, -60, 10978, 1)
        if hardware_version == 3:  #for version 3.1 NOT 3.0
            device.add_actuator(1, CCW, 2)
            device.add_actuator(2, CCW, 1)
            device.add_encoder(1, CCW, 168, 4880, 2)
            device.add_encoder(2, CCW, 12, 4880, 1)
        else:
            device.add_actuator(1, CCW, 2)
            device.add_actuator(2, CW, 1)
            device.add_encoder(1, CCW, 241, 10752, 2)
            device.add_encoder(2, CW, -61, 10752, 1)

        device.device_set_parameters()
        i = 0

        while (True):
            # if (haplyBoard.data_available()):
            #     device.device_read_data()
            #     motorAngle = device.get_device_angles()
            #     device_position = device.get_device_position(motorAngle)
            #     # print("Device position: " + str(device_position))
            #     i += 1
            # forces = [0, 10]
            # device.set_device_torques(forces)
            # device.device_write_torques()
            applyForces(haplyBoard,device, forceAmp,forceDir)
            time.sleep(0.001)


if __name__ == "__main__":
    main()
