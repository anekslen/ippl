import numpy as np
import os

def ReadData(filename):

    with open(filename, 'r') as file:
        lines = file.readlines()
        data = [line.strip().split(',') for line in lines]
    
    Rmin = float(data[1][1])
    Rmax = float(data[1][2])

    Thetamin = float(data[3][1])
    Thetamax = float(data[3][2])

    Z = [float(value) for value in data[5][1:]]

    Ekin = data[6][1]

    Vsplit = []
    for line in data[7:]:
        if line[0] == "Vsplit":
            Vsplit.append([int(num) for num in line[1:4]])

    return Rmin, Rmax, Thetamin, Thetamax, Z, Ekin, Vsplit


def CreateOutput(RMin, RMax, ThetaMin, ThetaMax, Z, Ekin, VSplit):
    npoints = np.zeros((len(Z), len(VSplit)))
    step = 0.0007

    ThetaMin_rad = ThetaMin / 360 * 2 * np.pi
    ThetaMax_rad = ThetaMax / 360 * 2 * np.pi

    x_range = [RMin * np.cos(ThetaMin_rad), RMax * np.cos(ThetaMax_rad), RMin * np.cos(ThetaMax_rad), RMax * np.cos(ThetaMin_rad)]
    y_range = [RMin * np.sin(ThetaMin_rad), RMax * np.sin(ThetaMax_rad), RMin * np.sin(ThetaMax_rad), RMax * np.sin(ThetaMin_rad)]

    xmin = min(x_range)
    ymin = min(y_range)

    xmax = max(x_range)
    ymax = max(y_range)

    for i, z in enumerate(Z):
        for j, Vsplit in enumerate(VSplit):
            outputfilename = f"eulerRuns/Grid_R{RMin}_R{RMax}_T{ThetaMin}_T{ThetaMax}_Z{z}_Vspile{Vsplit[0]}{Vsplit[1]}{Vsplit[2]}.csv"
            
            with open(outputfilename, 'w') as file:
                file.write("X,Y,Z,E_kin, Vsplit_x, Vsplit_y, Vsplit_z\n")

                for x in np.arange(xmin, xmax + step, step):
                    for y in np.arange(ymin, ymax + step, step):
                        r = np.sqrt(x**2 + y**2)
                        theta = np.arctan2(y, x) / 2 / np.pi * 360

                        if r >= RMin and r <= RMax and theta >= ThetaMin and theta <= ThetaMax:
                            npoints[i][j] += 1
                            file.write(f"{x},{y},{z},{Ekin},{Vsplit[0]},{Vsplit[1]},{Vsplit[2]}\n") 

            new_outputfilename = outputfilename.replace(".csv", f"_Np{int(npoints[i][j])}.csv")
            os.rename(outputfilename, new_outputfilename)

            print(f"Created {npoints[i][j]} points for {new_outputfilename}")

rmin, rmax, thetamin, thetamax, z, Ekin, vsplit = ReadData("pointdataInner.txt")
CreateOutput(rmin, rmax, thetamin, thetamax, z, Ekin, vsplit)

rmin, rmax, thetamin, thetamax, z, Ekin, vsplit = ReadData("pointdataOuter.txt")
CreateOutput(rmin, rmax, thetamin, thetamax, z, Ekin, vsplit)