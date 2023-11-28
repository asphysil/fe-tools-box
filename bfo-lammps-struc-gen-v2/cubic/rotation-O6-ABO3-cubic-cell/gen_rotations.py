import numpy as np 


def rotation_about_z(alpha_z):

	alpha_z=alpha_z*(np.pi/180.0) # deg to radian
	rmatrix = np.zeros((3,3))

	rmatrix[0][0] = np.cos(alpha_z)
	rmatrix[0][1] = np.sin(alpha_z)
	rmatrix[1][0] = -rmatrix[0][1] 
	rmatrix[1][1] = rmatrix[0][0]
	rmatrix[2][2] = 1.0
	return rmatrix

def rotation_about_y(alpha_y):
	
	alpha_y=alpha_y*(np.pi/180.0) # deg to radian
	rmatrix = np.zeros((3,3))

	rmatrix[0][0] = np.cos(alpha_y)
	rmatrix[0][2] = np.sin(alpha_y)
	rmatrix[2][0] = -rmatrix[0][2] 
	rmatrix[2][2] = rmatrix[0][0]
	rmatrix[1][1] = 1.0
	return rmatrix
def rotation_about_x(alpha_x):

	alpha_x=alpha_x*(np.pi/180.0) # deg to radian
	rmatrix = np.zeros((3,3))

	rmatrix[1][1] = np.cos(alpha_x)
	rmatrix[1][2] = -np.sin(alpha_x)
	rmatrix[2][1] = -rmatrix[1][2] 
	rmatrix[2][2] = rmatrix[1][1]
	rmatrix[0][0] = 1.0
	return rmatrix	

def rotation_about_u(alpha, uvec):
    alpha_rad=alpha*(np.pi/180.0) # deg to radian
    rmatrix = np.zeros((3,3))
    a = np.cos(alpha_rad)
    b = np.sin(alpha_rad)
    ux = uvec[0]
    uy = uvec[1]
    uz = uvec[2]

    rmatrix[0][0] = a + (ux**2*(1.0-a))
    rmatrix[0][1] = (ux*uy*(1-a))- (uz*b) 
    rmatrix[0][2] = (ux*uz*(1-a)) + uy*b

    rmatrix[1][0] = (uy*ux*(1.0-a)) + uz*b
    rmatrix[1][1] = a + (uy**2*(1.0-a))
    rmatrix[1][2] = (uy*uz*(1-a))-(ux*b)

    rmatrix[2][0] = (uz*ux*(1.0-a))-(uy*b)
    rmatrix[2][1] = (uz*uy*(1.0-a)) + (ux*b)
    rmatrix[2][2] = a + (uz**2*(1.0-a))

    return rmatrix

def rotation_vec(rmatrix, vec):
    rvec = np.matmul(rmatrix,vec)
    return rvec 

# # angle in degree
# alpha_x = 60
# alpha_y = 60
# alpha_z = 60


# matx=rotation_about_x(alpha_x)
# maty=rotation_about_y(alpha_y)
# matz=rotation_about_z(alpha_z)

# A = np.matmul(maty,matz)
# B = np.matmul(matx, A)
# v = np.zeros(3)

# v[0]=1
# v[1]=1
# v[2]=1

# #print(v, maty)

# l = np.linalg.norm(v)
# uvec = v/l 

# C = rotation_about_u(alpha_x, uvec)

# print(np.matmul(v, C))
