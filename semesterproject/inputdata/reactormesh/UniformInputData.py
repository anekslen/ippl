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

    return Rmin, Rmax, npointsR, Thetamin, Thetamax, npointsTheta, Z, Ekin, Vsplit


def CreateOutput(RMin, RMax, npointsR, ThetaMin, ThetaMax, npointsTheta, Z, Ekin, VSplit):

    RPoints = np.linspace(RMin, RMax, npointsR)
    ThetaPoints = np.linspace(ThetaMin / 360 * 2 * np.pi ,ThetaMax / 360 * 2 * np.pi, npointsTheta)

    NPoints = npointsR * npointsTheta

    for z in Z:
        for Vsplit in VSplit:
            outputfilename = f"eulerRuns/Uniform_R{RMin}_R{RMax}_T{ThetaMin}_T{ThetaMax}_Z{z}_Vspile{Vsplit[0]}{Vsplit[1]}{Vsplit[2]}_Np{NPoints}.csv"
            
            with open(outputfilename, 'w') as file:
                file.write("R,Theta,Z,E_kin, Vsplit_x, Vsplit_y, Vsplit_z\n")
                
                for r in RPoints:
                    for theta in ThetaPoints:
                        file.write(f"{r},{theta},{z},{Ekin},{Vsplit[0]},{Vsplit[1]},{Vsplit[2]}\n")

rmin, rmax, npointsr, thetamin, thetamax, npointstheta, z, Ekin, vsplit = ReadData("pointdataInner.txt")
CreateOutput(rmin, rmax, npointsr, thetamin, thetamax, npointstheta, z, Ekin, vsplit)

rmin, rmax, npointsr, thetamin, thetamax, npointstheta, z, Ekin, vsplit = ReadData("pointdataOuter.txt")
CreateOutput(rmin, rmax, npointsr, thetamin, thetamax, npointstheta, z, Ekin, vsplit)