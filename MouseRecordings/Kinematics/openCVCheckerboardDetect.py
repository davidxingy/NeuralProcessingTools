import numpy as np
import cv2 as cv
import glob
import matplotlib as mpl

# termination criteria
criteria = (cv.TERM_CRITERIA_EPS + cv.TERM_CRITERIA_MAX_ITER, 30, 0.001)

# prepare object points, like (0,0,0), (1,0,0), (2,0,0) ....,(6,5,0)
objp = np.zeros((8*6,3), np.float32)
objp[:,:2] = np.mgrid[0:6,0:8].T.reshape(-1,2)

# Arrays to store object points and image points from all the images.
objpoints = [] # 3d point in real world space
imgpoints = [] # 2d points in image plane.
usedFrames = [] # which images had checkerboard detection

#images = glob.glob('X:\David\ArenaRecordings\D050-120925-ArenaRecording\Video\TmpCalibrationFrameImages\*.tiff')
images = range(2360)


for imNum in images:
    
    fName = "X:\David\ArenaRecordings\D054-012126-ArenaRecording\Video\TmpCalibrationFrameImages\cam2_usedframe_" + str(imNum+1) + ".tiff"
    img = cv.imread(fName)
    gray = cv.cvtColor(img, cv.COLOR_BGR2GRAY)

    # Find the chess board corners
    ret, corners = cv.findChessboardCorners(gray, (6,8), None)

    # If found, add object points, image points (after refining them)
    if ret == True:
        objpoints.append(objp)

        corners2 = cv.cornerSubPix(gray,corners, (11,11), (-1,-1), criteria)
        imgpoints.append(corners2)
        usedFrames.append(imNum+1)
        
        # Draw and display the corners
        cv.drawChessboardCorners(img, (6, 8), corners2, ret)
        cv.imshow("img", img)
        cv.waitKey(200)
        cv.destroyAllWindows()
        
        #fig, ax = mpl.pyplot.subplots()
        #ax.imshow(img)
        #mpl.pyplot.plot(corners2[:,0,0],corners2[:,0,1],'.')
        
output=dict();
output['Points'] = imgpoints
output['Frames'] = usedFrames
  
np.save("Points.npy",imgpoints)
np.save("Frames.npy",usedFrames)
