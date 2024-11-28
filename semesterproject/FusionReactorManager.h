#ifndef IPPL_FUSION_REACTOR_MANAGER_H
#define IPPL_FUSION_REACTOR_MANAGER_H

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkCellArray.h>
#include <vtkIdList.h>
#include <vtkGradientFilter.h>
#include <vtkArrayCalculator.h>

#include <vector>

#include "UnstructuredFieldContainer.hpp"
#include "UnstructuredGridManager.h"
#include "UParticleContainer.h"

template <typename T, unsigned Dim>
class FusionReactorManager
    : public UnstructuredGridManager<T, Dim, ParticleContainer<T, Dim>>, UnstructuredFieldContainer<T, Dim>> {
public:
    using UParticleContainer_t = UParticleContainer<T, Dim>;
    using UnstructuredFieldContainer_t = UnstructuredFieldContainer<T, Dim>;

    FusionReactorManager(size_type totalP_, int nt_, std::string& stepMethod_)
        : UnstructuredGridManager<T, Dim, ParticleContainer<T, Dim>, UnstructuredFieldContainer<T, Dim>>(size_type totalP_, int nt_, std::string& stepMethod_)

    ~FusionReactorManager(){}

    void pre_run(const char* filename) override {
        Inform m("Pre Run");

        this->dt_m     = 0.05;  // TODO: how to choose timestep?
        this->it_m     = 0;
        this->time_m   = 0.0;

        m << "Discretization:" << endl
          << "nt " << this->nt_m << " Np= " << this->totalP_m << endl;

        // Create an unstructured field container
        this->setUFieldContainer( std::make_shared<UnstructuredFieldContainer_t>(filename) );

        // Initialize variables for dummy mesh
        this->nr_m = 100;                                                       // TODO: this is chosen randomly
        for (unsigned i = 0; i < Dim; i++) {
            this->domain_m[i] = ippl::Index(this->nr_m[i]);
        }
        this->UFieldContainer_t->getGridBounds(&this->rmin_m, &this->rmax_m);   // set mesh bounds to the bounds of the unstructured grid TODO: check this
        this->hr_m = this->rmax_m / this->nr_m;                                 // calculate mesh spacing
        this->origin_m = this->rmin_m;

        m << "Origin = " << this->origin_m << " hr = " << this->hr_m << " rmin = " << this->rmin_m << " rmax = " << this->rmax_m << endl;

        // Create mesh
        this->mesh_m(this->domain_m, this->hr_m, this->origin_m);

        // Initialize variables for dummy field layout
        this->decomp_m = {false, false, false};                                 // No parallel decomposition needed
        this->isAllPeriodic_m = false;

        // Create dummy field layout
        this->fl_m(MPI_COMM_WORLD, this->domain_m, this->decomp_m, this->isAllPeriodic_m);

        this->setParticleContainer( std::make_shared<ParticleContainer_t>( this->mesh_m, this->fl_m) );

        initializeParticles();

        // Get B value at the initial particle positions    // TODO: not needed probably
        this->pcontainer_m->B = ufc->getB(this->pcontainer_m->RReal);

        m << "Initialization done";
    }

    void initializeParticles(){
        Inform m("Initialize Particles");

        std::shared_ptr<ParticleContainer_t> pc = this->pcontainer_m;
        std::shared_ptr<UnstructuredFieldContainer_t> ufc    = this->ufcontainer_m;

        static IpplTimings::TimerRef particleCreation = IpplTimings::getTimer("particlesCreation");
        IpplTimings::startTimer(particleCreation);

        pc->create(this->totalP_m);

        assert(this->totalP_m == 1);    // TODO: change this to read initial particle positions from a file

        // Set initial particle positions
        pc->RCylReal(0) = Vector_t(0.318, M_PI / 2.0, 0.43);                // Initial real position in cylindrical coordinates
        pc->RReal(0) = Vector_t(                                            // Initial real position in cartesian coordinates
            pc->RCylReal(0)[0] * Kokkos::cos(pc->RCylReal(0)[1]),
            pc->RCylReal(0)[0] * Kokkos::sin(pc->RCylReal(0)[1]),
            pc->RCylReal(0)[2]
        );

        if(static_cast<int>(pc->RCylReal(0)[1] / (M_PI / 4.0)) % 2 == 0){   // Initial reference position in cylindrical coordinates
            pc->RCylRef(0) = Vector_t(pc->RCylReal(0)[0], Kokkos::fmod(pc->RCylReal(0)[1], M_PI / 2.0), Kokkos::abs(pc->RCylReal(0)[2]));
        }
        else{
            pc->RCylRef(0) = Vector_t(pc->RCylReal(0)[0], M_PI / 4.0 - Kokkos::fmod(pc->RCylReal(0)[1], M_PI / 4.0), Kokkos::abs(pc->RCylReal(0)[2]));
        }

        pc->R(0) = pc->RReal(0) = Vector_t(                 // Initial reference position in cartesian coordinates
            pc->RCylRef(0)[0] * Kokkos::cos(pc->RCylRef(0)[1]),
            pc->RCylRef(0)[0] * Kokkos::sin(pc->RCylRef(0)[1]),
            pc->RCylRef(0)[2]
        );

        // Check if the initial particle position is inside the domain
        int cellId = ufc->GetGridCell(pc->R(0));
        assert(cellId != -1);

        // Set initial particle charge in C and mass in kg
        pc->Q(0) = 1.602e-19;  
        pc->mass(0) = 1.673e-27;

        // Set initial particle energy in J
        pc->Ek(0) = 100 * 1.602e-19;

        // Set initial particle velocity in m/s
        double absv = Kokkos::sqrt(2.0 * pc->Ek(0) / pc->mass(0));
        pc->V(0) = Vector_t(absv/2, absv/2, 0.0);

        // Calculate initial change of velocity
        ufc->interpolateField("B_field", pc->R(0), pc->B(0), cellId);
        pc->dV(0) = pc->Q(0) * pc->V(0).cross(pc->B(0));

        IpplTimings::stopTimer(particleCreation);
        m << "particles created and initial conditions assigned " << endl;
    }

    void advance() override {
        if (this->stepMethod_m == "LeapFrog") {
            LeapFrogStep();
        }
	else{
            throw IpplException(TestName, "Step method is not set/recognized!");
        }
    }

    void LeapFrogStep(){
        // LeapFrog time stepping https://en.wikipedia.org/wiki/Leapfrog_integration
        static IpplTimings::TimerRef PTimer           = IpplTimings::getTimer("pushVelocity");
        static IpplTimings::TimerRef RTimer           = IpplTimings::getTimer("pushPosition");

        double dt                               = this->dt_m;
        std::shared_ptr<ParticleContainer_t> pc = this->pcontainer_m;
        std::shared_ptr<UnstructuredFieldContainer_t> ufc    = this->ufcontainer_m;

        IpplTimings::startTimer(PTimer);
        pc->V = pc->V + 0.5 * dt * pc->dV;
        IpplTimings::stopTimer(PTimer);

        // drift
        IpplTimings::startTimer(RTimer);
        pc->RReal = pc->RReal + dt * pc->P;
        IpplTimings::stopTimer(RTimer);

        // Since the particles have moved calculate new reference position, check if particle is still in the grid
        IpplTimings::startTimer(updateTimer);
        // Position in cylindrical coordinates
        for(int i = 0; i < pc->getTotalNum(); ++i){
            pc->RCylReal(i) = Vector_t(
                Kokkos::sqrt(Kokkos::pow(pc->RReal(i)[0], 2) + Kokkos::pow(pc->RReal(i)[1], 2)),
                Kokkos::atan2(pc->RReal(i)[1], pc->RReal(i)[0]),
                pc->RReal(i)[2]
            );
        }

        // Reference position in cylindrical coordinates
        if(static_cast<int>(pc->RCylReal(0)[1] / (M_PI / 4.0)) % 2 == 0){
            for (int i = 0; i < pc->getTotalNum(); ++i) {
                pc->RCylRef(i) = Vector_t(pc->RCylReal(i)[0], Kokkos::fmod(pc->RCylReal(i)[1], M_PI / 2.0), Kokkos::abs(pc->RCylReal(i)[2]));
            }
        }
        else{
            for(int i = 0; i < pc->getTotalNum(); ++i){
                pc->RCylRef(i) = Vector_t(pc->RCylReal(i)[0], M_PI / 4.0 - Kokkos::fmod(pc->RCylReal(i)[1], M_PI / 4.0), Kokkos::abs(pc->RCylReal(i)[2]));
            }
        }

        // Reference position in cartesian coordinates
        for(int i = 0; i < pc->getTotalNum(); ++i){
            pc->R(i) = Vector_t(
                pc->RCylRef(i)[0] * Kokkos::cos(pc->RCylRef(i)[1]),
                pc->RCylRef(i)[0] * Kokkos::sin(pc->RCylRef(i)[1]),
                pc->RCylRef(i)[2]
            );
        }

        // Check if the particle position is inside the domain
        Vector_t<int, pc->getTotalNum()> cellIds;
        for(int i = 0; i < pc->getTotalNum(); ++i){
            cellIds[i] = ufc->GetGridCell(pc->R(i));
            assert(cellIds[i] != -1);
        }

        IpplTimings::stopTimer(updateTimer);

        // Calculate new magnetic field at the new particle position
        for(int i = 0; i < pc->getTotalNum(); ++i){
            ufc->interpolateField("B_field", pc->R(i), pc->B(i), cellIds[i]);
        }

        // Calculate new change of velocity
        for(int i = 0; i < pc->getTotalNum(); ++i){
            pc->dV(i) = pc->Q(i) * pc->V(i).cross(pc->B(i));
        }

        // kick
        IpplTimings::startTimer(PTimer);
        pc->P = pc->P + 0.5 * dt * pc->dV;
        IpplTimings::stopTimer(PTimer);
    }

    void dump() override {
        static IpplTimings::TimerRef dumpDataTimer = IpplTimings::getTimer("dumpData");
        IpplTimings::startTimer(dumpDataTimer);
        dumpPositions();
        IpplTimings::stopTimer(dumpDataTimer);
    }

    template <typename View>
    void dumpPositions() {
        Inform m("dumpPositions");

        for (int i = 0; i < this->pcontainer_m->size(); ++i) {
            m << "Time " << this->time_m << "Particle " << i << ": Position = " << this->pcontainer_m->R(i)
              << ", Velocity = " << this->pcontainer_m->V(i) << std::endl;
        }
};
#endif

/*
// Define functions to read the grid from a file
vtkSmartPointer<vtkUnstructuredGrid> readGrid(const char* filename)
{
    // Create a reader .vtk file
    vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    reader->SetFileName(filename);
    reader->Update();

    // Get the unstructured grid from the reader
    vtkSmartPointer<vtkUnstructuredGrid> grid = reader->GetOutput();
    assert(grid);
    if (!grid) {
        std::cerr << "Error: Failed to read the mesh from the file!" << std::endl;
        return NULL;
    }

    // Get the points (nodes) of the mesh
    vtkPoints* points = grid->GetPoints();
    if (!points) {
        std::cerr << "Error: Failed to get the points from the mesh!" << std::endl;
        return NULL;
    }

    // Get the B vector field array from the mesh
    vtkDataArray* bFieldArray = grid->GetPointData()->GetArray("B_Field");
    if (!bFieldArray) {
        std::cerr << "Error: Failed to get the B_Field array from the mesh!" << std::endl;
        return NULL;
    }

    // Get the cells (elements) of the mesh
    vtkCellArray* cells = grid->GetCells();
    if (!cells) {
        std::cerr << "Error: Failed to get the cells from the mesh!" << std::endl;
        return NULL;
    }

    return grid;
}

// Define functions to write a grid to a file
void writeGrid(vtkSmartPointer<vtkUnstructuredGrid> grid, const char* filename) {
    // Write the grid to a file
    vtkSmartPointer<vtkUnstructuredGridWriter> writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
    writer->SetFileName(filename);
    writer->SetInputData(grid);
    writer->Write();
}

// Define functions to write a vector field to a file
void writeField(vtkSmartPointer<vtkUnstructuredGrid> grid, const char* field_name, const char* filename, bool Magnitude = false, const char* magnitude_field_name = "Magnitude") {
    // Check that the field exists in the grid
    if (!grid->GetPointData()->HasArray(field_name)) {
        std::cerr << "Error: Field " << field_name << " does not exist in the grid!" << std::endl;
        return;
    }

    // Check that the field is a vector field
    if (grid->GetPointData()->GetArray(field_name)->GetNumberOfComponents() != 3) {
        std::cerr << "Error: Field " << field_name << " is not a vector field!" << std::endl;
        return;
    }

    // Check that the magnitude field exists in the grid
    if (Magnitude && !grid->GetPointData()->HasArray(magnitude_field_name)) {
        std::cerr << "Error: Magnitude field " << magnitude_field_name << " does not exist in the grid!" << std::endl;
        return;
    }

    // Check that the magnitude field is a scalar field
    if (Magnitude && grid->GetPointData()->GetArray(magnitude_field_name)->GetNumberOfComponents() != 1) {
        std::cerr << "Error: Magnitude field " << magnitude_field_name << " is not a scalar field!" << std::endl;
        return;
    }

    // Write the curl values to csv file
    std::ofstream file(filename);

    if (Magnitude) {
        file << "Point_Index," << field_name << "_X," << field_name << "_Y," << field_name << "_Z," << magnitude_field_name << std::endl;
    }
    else {
        file << "Point_Index," << field_name << "_X," << field_name << "_Y," << field_name << "_Z" << std::endl;
    }

    // Write the curl values to csv file
    for (int i = 0; i < grid->GetPoints()->GetNumberOfPoints(); i++) {
        double curl[3];
        double magnitude;
        grid->GetPointData()->GetArray(field_name)->GetTuple(i, curl);
        if (Magnitude) {
            magnitude = grid->GetPointData()->GetArray(magnitude_field_name)->GetTuple1(i);
        }

        // Set the format for floating-point values
        file << std::scientific << std::setprecision(4) << std::setw(10);

        if (Magnitude) {
            file << i << ","
                << std::setw(10) << curl[0] << ","
                << std::setw(10) << curl[1] << ","
                << std::setw(10) << curl[2] << ","
                << std::setw(10) << magnitude << std::endl;
        }
        else {
            file << i << ","
                << std::setw(10) << curl[0] << ","
                << std::setw(10) << curl[1] << ","
                << std::setw(10) << curl[2] << std::endl;
        }
    }
    file.close();
}

// Define a function to fint the gridcell in which a point is located
int CoordinateToGridCell(vtkSmartPointer<vtkUnstructuredGrid> grid, double x, double y, double z) {
    // Get the points (nodes) of the mesh
    vtkPoints* points = grid->GetPoints();
    if (!points) {
        std::cerr << "Error: Failed to get the points from the mesh!" << std::endl;
        return -1;
    }

    // Get the cells (elements) of the mesh
    vtkCellArray* cells = grid->GetCells();
    if (!cells) {
        std::cerr << "Error: Failed to get the cells from the mesh!" << std::endl;
        return -1;
    }

    // Find the cell that contains the point (x, y, z)
    vtkIdType cellId = -1;
    double p[3] = { x, y, z };
    double pcoords[3];
    double weights[8];
    int subId;
    double dist2 = -1.0;
    // Optimization if previous cell is known (set cellId and cell to previous cell)
    // TODO: Not thread save, check this !!!!!!
    cellId = grid->FindCell(p, NULL, 0, 0, dist2, subId, pcoords, weights);
    return cellId;
}

// Define functions to calculate the curl of the a vector field
vtkSmartPointer<vtkUnstructuredGrid> calculateCurl(vtkSmartPointer<vtkUnstructuredGrid> grid, const char* field_name, const char* output_field_name = "Vorticity") {
    // Check that the field exists in the grid
    if (!grid->GetPointData()->HasArray(field_name)) {
        std::cerr << "Error: Field " << field_name << " does not exist in the grid!" << std::endl;
        return NULL;
    }

    // Check that the field is a vector field
    if (grid->GetPointData()->GetArray(field_name)->GetNumberOfComponents() != 3) {
        std::cerr << "Error: Field " << field_name << " is not a vector field!" << std::endl;
        return NULL;
    }

    // Get the points (nodes) of the mesh
    vtkPoints* points = grid->GetPoints();
    if (!points) {
        std::cerr << "Error: Failed to get the points from the mesh!" << std::endl;
        return NULL;
    }

    // Get the B vector field array from the mesh
    vtkDataArray* bFieldArray = grid->GetPointData()->GetArray(field_name);
    if (!bFieldArray) {
        std::cerr << "Error: Failed to get the " << field_name << " array from the mesh!" << std::endl;
        return NULL;
    }

    // Calculate the curl of the vector field at each point
    vtkSmartPointer<vtkGradientFilter> gradientFilter = vtkSmartPointer<vtkGradientFilter>::New();
    gradientFilter->SetInputData(grid);
    gradientFilter->SetInputScalars(vtkDataObject::FIELD_ASSOCIATION_POINTS, field_name);
    gradientFilter->SetComputeGradient(false);   // Disable the gradient computation
    gradientFilter->SetComputeVorticity(true);   // Enable the curl computation
    gradientFilter->Update();

    // Get the output unstructured grid from the gradient filter
    vtkSmartPointer<vtkUnstructuredGrid> curlGrid = gradientFilter->GetUnstructuredGridOutput();

    // Get the array from the curl grid
    vtkDataArray* curlArray = curlGrid->GetPointData()->GetArray("Vorticity");
    if (!curlArray) {
        std::cerr << "Error: Failed to get the curl array from the mesh!" << std::endl;
        return NULL;
    }
    else {
        curlArray->SetName(output_field_name);
    }

    return curlGrid;
}
// Define functions to calculate the magnitude of a vector field
void calculateMagnitude(vtkSmartPointer<vtkUnstructuredGrid> grid, const char* field_name, const char* output_field_name = "Magnitude") {
    // Check that the field exists in the grid
    if (!grid->GetPointData()->HasArray(field_name)) {
        std::cerr << "Error: Field " << field_name << " does not exist in the grid!" << std::endl;
        return;
    }

    // Check that the field is a vector field
    if (grid->GetPointData()->GetArray(field_name)->GetNumberOfComponents() != 3) {
        std::cerr << "Error: Field " << field_name << " is not a vector field!" << std::endl;
        return;
    }

    // Calculate the magnitude of the curl at each point
    std::string function = "mag(" + std::string(field_name) + ")";

    vtkSmartPointer<vtkArrayCalculator> magnitudeCalculator = vtkSmartPointer<vtkArrayCalculator>::New();
    magnitudeCalculator->SetInputData(grid);
    magnitudeCalculator->SetAttributeTypeToPointData();
    magnitudeCalculator->AddVectorArrayName(field_name);
    magnitudeCalculator->SetFunction(function.c_str());    // Calculate the magnitude of the curl
    magnitudeCalculator->SetResultArrayName(output_field_name);
    magnitudeCalculator->Update();

    // Add the magnitude of the curl to the curl grid
    grid->GetPointData()->AddArray(magnitudeCalculator->GetUnstructuredGridOutput()->GetPointData()->GetArray(output_field_name));
}
*/






/*
// Interpolate a field at a given point
void interpolateField(vtkSmartPointer<vtkUnstructuredGrid> grid, const char* field_name, double x, double y, double z, std::vector<int> &interpolatedField, int cellID = -1, double weights[8] = NULL)
{
    // Check that the field exists in the grid
    if (!grid->GetPointData()->HasArray(field_name)) {
        std::cerr << "Error: Field " << field_name << " does not exist in the grid!" << std::endl;
        return;
    }

    // Check that the field is a vector field
    if (grid->GetPointData()->GetArray(field_name)->GetNumberOfComponents() != 3) {
        std::cerr << "Error: Field " << field_name << " is not a vector field!" << std::endl;
        return;
    }

    if (cellID == -1 || weights == NULL) {
        // Find the cell that contains the point (x, y, z)
        double p[3] = { x, y, z };
        double pcoords[3];
        double weights[8];
        int subId;
        double dist2 = -1.0;
        // Optimization if previous cell is known (set cellId and cell to previous cell)
        // TODO: Not thread save, check this !!!!!!
        cellID = grid->FindCell(p, NULL, 0, 0, dist2, subId, pcoords, weights);

        // Access the field
        vtkDataArray* fieldArray = grid->GetPointData()->GetArray(field_name);
        if (!fieldArray) {
            std::cerr << "Error: Failed to get the " << field_name << " array from the mesh!" << std::endl;
            return;
        }

        // Get the point indices of the cell in which the point (x, y, z) is located
        vtkIdList* pointIds = grid->GetCell(cellID)->GetPointIds();
        
        // Interpolate the field at the point (x, y, z)
        vtkIdType interpolationValues[3];

        fieldArray->InterpolateTuple(interpolationValues, pointIds, fieldArray, weights);

        for(int i = 0; i < 3; i++) {
            interpolatedField[i] = interpolationValues[i];
        }
        
        return;
    }
}
*/

#endif