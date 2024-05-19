# Visualizing Haptics, a Debugging Tool for Learners

The python code provided can be used to simulate and visualise static haptic effects and interactions in 2D and 3D.
Still a WIP and can only be applied on round objects for now. Additional functionality will be added. 

## 2D Visualization: 
![alt text](https://github.com/emmanuelwilson/VisualizingHaptics/blob/main/2Drepresentation.PNG)

## 3D Visualization: 
![alt text](https://github.com/emmanuelwilson/VisualizingHaptics/blob/main/3Drepresentation.gif)

## Forces: 

### Springs: 
``` circleSpring(Xc, Yc, x radius, y radius, k, tilt angle) ``` 

Allows you to define a spring force at location (Xc, Yc) with a radius of your choice (can define 2 radii to form an ellipse), a force value k, and tilt or rotation angle of the shape.
Which will output a N by M by 2 matrix with force amplitude and direction of force. 

### Friction: 
``` circleFriction(Xc,Yc,rx,ry, friction ,tilt) ```

Allows you to create an area of friction within the haptic space. The center of the circle will be located at (Xc,Yc) with radius/radii rx and ry with a friction value between 0-1 and a tilt or rotation angle of the shape. 
This will output an N by M matrix which contains the friction values of this shape. 

### Texture: 
``` circleTextureScatter(Xc,Yc,rx,ry, textMax, textureDensity,tilt) ```

Allows you to define a circle of texture located at (Xc,Yc) with radius/radii rx/ry, with a maximum force of textMax and a event density of textureDensity and a rotation of the shape of tilt. 
This function is a placeholder until a more indepth view of the texture interactions can be made. 

## Interactions: 

There are specific functions used for adding forces together in order to properly view them. 

For 2D visualization of spings:
``` addForces ``` 

For 3D visualization of springs: 
``` addForcesSign ```

For Firction: 
``` addFriction ```


