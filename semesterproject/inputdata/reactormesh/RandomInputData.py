import numpy as np

def ReadData(filename):

    with open(filename, 'r') as file:
        lines = file.readlines()
        data = [line.strip().split(',') for line in lines]
    
    Rmin = float(data[1][1])
    Rmax = float(data[1][2])
    npointsR = int(data[2][1])

    Thetamin = float(data[3][1])
    Thetamax = float(data[3][2])
    npointsTheta = int(data[4][1])

    Z = [float(value) for value in data[5][1:]]

    Ekin = data[6][1]

    Vsplit = []
    for line in data[7:]:
        if line[0] == "Vsplit":
            Vsplit.append([int(num) for num in line[1:4]])

    nPoints = npointsR * npointsTheta

    return Rmin, Rmax, Thetamin, Thetamax, Z, Ekin, Vsplit, nPoints


def CreateOutput(RMin, RMax, ThetaMin, ThetaMax, Z, Ekin, VSplit, NPoints):
    ThetaMin_rad = ThetaMin / 360 * 2 * np.pi
    ThetaMax_rad = ThetaMax / 360 * 2 * np.pi

    x = [RMin * np.cos(ThetaMin_rad), RMax * np.cos(ThetaMax_rad), RMin * np.cos(ThetaMax_rad), RMax * np.cos(ThetaMin_rad)]
    y = [RMin * np.sin(ThetaMin_rad), RMax * np.sin(ThetaMax_rad), RMin * np.sin(ThetaMax_rad), RMax * np.sin(ThetaMin_rad)]

    xmin = min(x)
    ymin = min(y)

    xmax = max(x)
    ymax = max(y)

    for z in Z:
        for Vsplit in VSplit:
            outputfilename = f"eulerRuns/Rand_R{RMin}_R{RMax}_T{ThetaMin}_T{ThetaMax}_Z{z}_Vspile{Vsplit[0]}{Vsplit[1]}{Vsplit[2]}_Np{NPoints}.csv"
            
            with open(outputfilename, 'w') as file:
                file.write("R,Theta,Z,E_kin, Vsplit_x, Vsplit_y, Vsplit_z\n")

                npoints = 0

                while npoints < NPoints:
                    x = np.random.uniform(xmin, xmax)
                    y = np.random.uniform(ymin, ymax)

                    r = np.sqrt(x**2 + y**2)
                    theta = np.arctan2(y, x) / 2 / np.pi * 360

                    if r >= RMin and r <= RMax and theta >= ThetaMin and theta <= ThetaMax:
                        npoints += 1
                        file.write(f"{r},{theta / 360 * 2 * np.pi},{z},{Ekin},{Vsplit[0]},{Vsplit[1]},{Vsplit[2]}\n")

rmin, rmax, thetamin, thetamax, z, Ekin, vsplit, npoints = ReadData("pointdataInner.txt")
CreateOutput(rmin, rmax, thetamin, thetamax, z, Ekin, vsplit, npoints)

rmin, rmax, thetamin, thetamax, z, Ekin, vsplit, npoints = ReadData("pointdataOuter.txt")
CreateOutput(rmin, rmax, thetamin, thetamax, z, Ekin, vsplit, npoints)