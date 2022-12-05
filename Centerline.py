import math
import nrrd
import numpy  as np
from matplotlib import pyplot as plt

def create_tables(RAYS):
    sintable = []
    costable = []
    for i in range(RAYS + 1):
        sintable.append(math.sin(i/(180/math.pi)))
        costable.append(math.cos(i/(180/math.pi)))
    return sintable, costable

def boundary (I, px, py, X, Y):
    old_pos = I[int(px), int(py)] 
    new_pos = I[int(X), int(Y)]

    diff = (old_pos - new_pos) if (old_pos > new_pos) else (new_pos - old_pos) 
            

    if(diff > 200):
        return True

    return False
        
def limitExceed(I, X, Y):
    if( X < 0 or 
        Y < 0 or 
        int(X) > I.shape[0]-1 or 
        int(Y) > I.shape[1]-1):
        return True
    return False 


def find_path(I, px, py, dir, dist):
    pts = []

    for angle in range(-89 , 90):
        ax = sintable[(dir + angle) % 360]
        ay = costable[(dir + angle) % 360]

        X = px + dist * ax
        Y = py + dist * ay

        if(boundary(I, px, py, X, Y) or limitExceed(I, X, Y)):
            continue

        pts.append([int(X), int(Y)])

    return pts
        

def find_center(pts):
    sumx = 0
    sumy = 0

    for [x,y] in pts:
        sumx = sumx + x
        sumy = sumy + y

    return int(sumx / len(pts)), int(sumy / len(pts))



def raycast_at_angle(I, px, py, dir, angle):
        
        ax = sintable[(dir + angle) % 360]
        ay = costable[(dir + angle) % 360]

        X = px 
        Y = py 
        
        while(1):
            
            X = X + ax
            Y = Y + ay

            if(int(X) > I.shape[0]-1  or int(Y) > I.shape[1]-1 or int(X) < 0 or int(Y) < 0):
                break
            
            old_pos = I[int(px), int(py)] 
            new_pos = I[int(X-ax), int(Y-ay)]

            #[on_true] if [expression] else [on_false] 

            diff = (old_pos - new_pos) if (old_pos > new_pos) else (new_pos - old_pos) 
            

            if(diff > 200):
                break
        
        return int(X-ax), int(Y-ay)



# Location : D:\HD1_Backup\14806362\KiTS\KiTS\K1

readdata, header = nrrd.read('D://HD1_Backup//14806362//KiTS//KiTS//K1//K1.nrrd')
print(readdata.shape)

#(512,512,1059)
#print(header)
#print(type(readdata))
#print(np.min(readdata))
#print(np.max(readdata))

plt.imshow(readdata[256,:,:], cmap='gray')
#Test Candidate
#plt.imshow(readdata[:,256,:], cmap='gray')
#plt.plot(400,250,'ro', markersize=1 )
plt.plot(650,310,'ro', markersize=1 )

slice = readdata[256,:,:]

direction = 180
sintable, costable = create_tables(360)

px = 310
py = 650

while(1):

    b1_x, b1_y = raycast_at_angle(slice, px, py, direction, -90)
    b2_x, b2_y = raycast_at_angle(slice, px, py, direction, +90)

    new_ptx = int((b1_x + b2_x) / 2)
    new_pty = int((b1_y + b2_y) / 2)

    print('Boundary 1:', b1_x, b1_y)
    print('Boundary 2:', b2_x, b2_y)
    print('Centre Candidate:',new_ptx, new_pty)

    dist = int(( math.fabs(b1_x - new_ptx) + math.fabs(b1_x - new_ptx) ) / 2)

    #plt.plot(b1_y,b1_x,'b*', markersize=1 )
    #plt.plot(b2_y,b2_x,'b*', markersize=1 )
    plt.plot(new_pty,new_ptx,'y*', markersize=1 )

    pts = find_path(slice, new_ptx, new_pty, 180, dist)

    if(len(pts) == 0):
        break

    '''
    for i in range(len(pts)):
        plt.plot(pts[i][1], pts[i][0], 'g*', markersize=1)
    '''
    px, py = find_center(pts)
    #plt.plot(py, px, 'ro', markersize=2)
    

plt.show()
