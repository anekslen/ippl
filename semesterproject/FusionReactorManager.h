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
    using DummyFieldContainer_t = FieldContainer<T, Dim>;

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


        // Create dummy mesh and field layout for the particle container

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

        m << "Initialization done";
    }

    void initializeParticles(){
        Inform m("Initialize Particles");
        static IpplTimings::TimerRef particleCreation = IpplTimings::getTimer("particlesCreation");
        IpplTimings::startTimer(particleCreation);

        this->pcontainer_m->create(this->totalP_m);

        assert(this->totalP_m == 1);    // TODO: change this to read initial particle positions from a file

        // Set initial particle positions
        this->pcontainer_m->RCylReal(0) = Vector_t(0.318, M_PI / 2.0, 0.43);                // Initial real position in cylindrical coordinates
        this->pcontainer_m->RReal(0) = Vector_t(                                            // Initial real position in cartesian coordinates
            this->pcontainer_m->RCylReal(0)[0] * Kokkos::cos(this->pcontainer_m->RCylReal(0)[1]),
            this->pcontainer_m->RCylReal(0)[0] * Kokkos::sin(this->pcontainer_m->RCylReal(0)[1]),
            this->pcontainer_m->RCylReal(0)[2]
        );

        if(static_cast<int>(this->pcontainer_m->RCylReal(0, 1) / (M_PI / 4.0)) % 2 == 0){   // Initial reference position in cylindrical coordinates
            this->pcontainer_m->RCylRef(0) = Vector_t(this->pcontainer_m->RCylReal(0)[0], Kokkos::fmod(this->pcontainer_m->RCylReal(0)[1], M_PI / 2.0), Kokkos::abs(this->pcontainer_m->RCylReal(0)[2]));
        }
        else{
            this->pcontainer_m->RCylRef(0) = Vector_t(this->pcontainer_m->RCylReal(0)[0], M_PI / 4.0 - Kokkos::fmod(this->pcontainer_m->RCylReal(0)[1], M_PI / 4.0), Kokkos::abs(this->pcontainer_m->RCylReal(0)[2]));
        }

        this->pcontainer_m->R(0) = this->pcontainer_m->RReal(0) = Vector_t(                 // Initial reference position in cartesian coordinates
            this->pcontainer_m->RCylRef(0)[0] * Kokkos::cos(this->pcontainer_m->RCylRef(0)[1]),
            this->pcontainer_m->RCylRef(0)[0] * Kokkos::sin(this->pcontainer_m->RCylRef(0)[1]),
            this->pcontainer_m->RCylRef(0)[2]
        );

        // Set initial particle charge in C and mass in kg
        this->pcontainer_m->Q(0) = 1.602e-19;  
        this->pcontainer_m->mass(0) = 1.673e-27;

        // Set initial particle energy in J
        this->pcontainer_m->Ek(0) = 100 * 1.602e-19;

        // Set initial particle velocity in m/s
        double absv = Kokkos::sqrt(2.0 * this->pcontainer_m->Ek(0) / this->pcontainer_m->mass(0));
        this->pcontainer_m->V(0) = Vector_t(absv/2, absv/2, 0.0);

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
        std::shared_ptr<UnstructuredFieldContainer_t> fc    = this->fcontainer_m;

        IpplTimings::startTimer(PTimer);
        pc->P = pc->P - 0.5 * dt * pc->E;
        IpplTimings::stopTimer(PTimer);

        // drift
        IpplTimings::startTimer(RTimer);
        pc->R = pc->R + dt * pc->P;
        IpplTimings::stopTimer(RTimer);

        // Since the particles have moved spatially update them to correct processors
        IpplTimings::startTimer(updateTimer);
        pc->update();
        IpplTimings::stopTimer(updateTimer);

        size_type totalP        = this->totalP_m;
        int it                  = this->it_m;
        bool isFirstRepartition = false;
        if (this->loadbalancer_m->balance(totalP, it + 1)) {
                IpplTimings::startTimer(domainDecomposition);
                auto* mesh = &fc->getRho().get_mesh();
                auto* FL = &fc->getFL();
                this->loadbalancer_m->repartition(FL, mesh, isFirstRepartition);
                IpplTimings::stopTimer(domainDecomposition);
        }

        // scatter the charge onto the underlying grid
        this->par2grid();

        // Field solve
        IpplTimings::startTimer(SolveTimer);
        this->fsolver_m->runSolver();
        IpplTimings::stopTimer(SolveTimer);

        // gather E field
        this->grid2par();

        // kick
        IpplTimings::startTimer(PTimer);
        pc->P = pc->P - 0.5 * dt * pc->E;
        IpplTimings::stopTimer(PTimer);
    }

    void dump() override {
        static IpplTimings::TimerRef dumpDataTimer = IpplTimings::getTimer("dumpData");
        IpplTimings::startTimer(dumpDataTimer);
        dumpLandau(this->fcontainer_m->getE().getView());
        IpplTimings::stopTimer(dumpDataTimer);
    }

    template <typename View>
    void dumpLandau(const View& Eview) {
        const int nghostE = this->fcontainer_m->getE().getNghost();

        using index_array_type = typename ippl::RangePolicy<Dim>::index_array_type;
        double localEx2 = 0, localExNorm = 0;
        ippl::parallel_reduce(
            "Ex stats", ippl::getRangePolicy(Eview, nghostE),
            KOKKOS_LAMBDA(const index_array_type& args, double& E2, double& ENorm) {
                // ippl::apply<unsigned> accesses the view at the given indices and obtains a
                // reference; see src/Expression/IpplOperations.h
                double val = ippl::apply(Eview, args)[0];
                double e2  = Kokkos::pow(val, 2);
                E2 += e2;

                double norm = Kokkos::fabs(ippl::apply(Eview, args)[0]);
                if (norm > ENorm) {
                    ENorm = norm;
                }
            },
            Kokkos::Sum<double>(localEx2), Kokkos::Max<double>(localExNorm));

        double globaltemp = 0.0;
        ippl::Comm->reduce(localEx2, globaltemp, 1, std::plus<double>());

        double fieldEnergy =
            std::reduce(this->fcontainer_m->getHr().begin(), this->fcontainer_m->getHr().end(), globaltemp, std::multiplies<double>());

        double ExAmp = 0.0;
        ippl::Comm->reduce(localExNorm, ExAmp, 1, std::greater<double>());

        if (ippl::Comm->rank() == 0) {
            std::stringstream fname;
            fname << "data/FieldLandau_";
            fname << ippl::Comm->size();
            fname << "_manager";
            fname << ".csv";
            Inform csvout(NULL, fname.str().c_str(), Inform::APPEND);
            csvout.precision(16);
            csvout.setf(std::ios::scientific, std::ios::floatfield);
            if ( std::fabs(this->time_m) < 1e-14 ) {
                csvout << "time, Ex_field_energy, Ex_max_norm" << endl;
            }
            csvout << this->time_m << " " << fieldEnergy << " " << ExAmp << endl;
        }
        ippl::Comm->barrier();
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