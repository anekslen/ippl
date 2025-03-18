import numpy as np

def ReadData(filename):

    with open(filename, 'r') as file:
        lines = file.readlines()
        data = [line.strip().split(',') for line in lines]
    
    Rmin = float(data[1][1])
    Rmax = float(data[1][2])
    npointsR = int(data[2][1])

    R = np.linspace(Rmin, Rmax, npointsR)

    Thetamin = float(data[3][1]) / 360 * 2 * np.pi
    Thetamax = float(data[3][2]) / 360 * 2 * np.pi
    npointsTheta = int(data[4][1])

    Theta = np.linspace(Thetamin, Thetamax, npointsTheta)

    Z = [float(value) for value in data[5][1:]]

    Ekin = [float(value) for value in data[6][1:]]

    Vsplit = []
    for line in data[7:]:
        if line[0] == "Vsplit":
            Vsplit.append([int(num) for num in line[1:4]])

    return R, Theta, Z, Ekin, Vsplit


def CreateOutput(outputfilename, RPoints, ThetaPoints, ZPoints, Ekin, Vsplit):
    with open(outputfilename, 'w') as file:
        file.write("R,Theta,Z,E_kin, Vsplit_x, Vsplit_y, Vsplit_z\n")
        
        for i in RPoints:
            for j in ThetaPoints:
                for k in ZPoints:
                    for e in Ekin:
                        for v in Vsplit:

                            file.write(f"{i},{j},{k},{e},{v[0]},{v[1]},{v[2]}\n")

r, theta, z, Ekin, vsplit = ReadData("pointdata.txt")

CreateOutput("output.csv", r, theta, z, Ekin, vsplit)