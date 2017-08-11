#include "graph.h"
#include "program.h"
#include "visualizer/visualizer.h"
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

using namespace simit;


void make_mesh(int nx, int ny, double lx, double ly, Set& points, Set& hinges, Set& faces) {
// create a uniform mesh of nx x ny gridpoints over a lx x ly rectangle

    vector<vector<ElementRef> > ref_node(nx, vector<ElementRef>(ny));

    FieldRef<double,3> x = points.getField<double,3>("x");
    FieldRef<bool> fixed = points.getField<bool>("fixed");
    FieldRef<int>    idx = points.getField<int>("idx");

    int cnt = 0;
    for(int j=0; j<ny; j++) {
        for(int i=0; i<nx; i++) {
            ref_node[i][j]    = points.add();
            x(ref_node[i][j]) = { lx*i/(nx-1), ly*j/(ny-1), 0.0 };
            idx(ref_node[i][j]) = cnt;
            fixed(ref_node[i][j]) = false;
            cnt++;
        }
    }

    // connect gridpoints by making triangles of them
    // 
    //  (i,j+1)  (i+1, j+1)
    //       +-------+-------+
    //       |      /|      /|
    //       |    ,' |    ,' |
    //       |   /   |   /   |
    //       | ,'    | ,'    |
    //       +-------+-------+ 
    //  (i,j)    (i+1, j)
    //
    for(int j=0; j<ny-1; j++) {
        for(int i=0; i<nx-1; i++) {
            ElementRef node1 = ref_node[ i ][ j ];
            ElementRef node2 = ref_node[i+1][ j ];
            ElementRef node3 = ref_node[i+1][j+1];
            ElementRef element1 = faces.add(node1, node2, node3);

            ElementRef node4 = ref_node[i+1][j+1];
            ElementRef node5 = ref_node[ i ][j+1];
            ElementRef node6 = ref_node[ i ][ j ];
            ElementRef element2 = faces.add(node4, node5, node6);

            if (j>0) {
              ElementRef node7 = ref_node[ i ][j-1];
              hinges.add(node1, node2, node3, node7);
            }
            if (i<nx-2) {
              ElementRef node8 = ref_node[i+2][j+1];
              hinges.add(node2, node3, node1, node8);
            }
            
            hinges.add(node3, node1, node2, node5);

            

            
        }
    }

    // fix boundary nodes
    for(int i=0; i<nx; i++) 
        fixed(ref_node[i][0]) = true;      // south edge
    for(int i=0; i<nx; i++) 
        fixed(ref_node[i][ny-1]) = true;  // north edge

}


/************************************************************************//**
 * \brief Dump results as a surface plot in an ASCII-vtk file
 * \param filename    [in] VTK-filename. Use paraview to view this later
 * \param node        [in] Simit set to hold the nodes
 * \param faces       [in] Simit set to hold the elements (connects 3 nodes)
 ***************************************************************************/
void plot_results(string filename, Set& points, Set& faces) {
    // full documentation of the VTK file format can be found at
    // http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
    
    // fetch fields
    FieldRef<double,3> x = points.getField<double,3>("x");
    FieldRef<double,3> u = points.getField<double,3>("u");
    FieldRef<int>      idx = points.getField<int>("idx");

    // write header
    ofstream out(filename);
    out << "# vtk DataFile Version 3.0" << endl;
    out << "simit results " << endl;
    out << "ASCII " << endl;
    out << "DATASET UNSTRUCTURED_GRID " << endl << endl;

    // write node coordinates 
    out << "POINTS " << points.getSize() << " float" << endl;
    for(auto n : points)
        out << x(n)(0) << " " << x(n)(1) << " " << x(n)(2) << endl;
    out << endl;

    // write element connections
    out << "CELLS " << faces.getSize() << " " << faces.getSize()*4 << endl;
    for(auto f : faces) {
        out << "3 ";
        for(int j=0; j<3; j++) {
            ElementRef pt = faces.getEndpoint(f, j);
            out << idx(pt) << " ";
        }
        out << endl;
    }
    out << endl;

    // tag all elements as triangles
    out << "CELL_TYPES " << faces.getSize() << endl;
    for(auto f : faces)
        out << "5" << endl; // VTK_TRIANGLE element type
    out << endl;

    // write FEM solution
    out << "POINT_DATA " << points.getSize() << endl;
    out << endl;
    out << "SCALARS solution float" << endl;
    out << "LOOKUP_TABLE default" << endl;
    for(auto n : points)
        out << sqrt(u(n)(0)*u(n)(0) + u(n)(1)*u(n)(1) + u(n)(2)*u(n)(2)) << endl;
    out << endl;

}

void compute_loc_dir(Set& faces) {

    // fetch fields
    FieldRef<double,3,2> loc_dir  = faces.getField<double,3,2>("loc_dir");

    for(auto f : faces) { 
        loc_dir(f) = {1.0, 0.0,
                      0.0, 1.0,
                      0.0, 0.0};
    }
    

}

void set_props(Set& faces, Set& hinges) {

    // fetch fields
    FieldRef<double> Nu = faces.getField<double>("Nu");
    FieldRef<double> E = faces.getField<double>("E");
    FieldRef<double> Rho = faces.getField<double>("Rho");

    FieldRef<double> B = hinges.getField<double>("B");

    for(auto f : faces) { 
        Nu(f) = 0.3;
        E(f) = 1.0e2;
        Rho(f) = 1.0;
    }    

    for(auto h : hinges) { 
        B(h) = 1.0e-4;
    }    

}

int main(int argc, char **argv)
{

    init("cpu", sizeof(double));

    // Create a graph
    Set points;
    Set hinges(points, points, points, points);
    Set faces(points, points, points);

    // The fields of the points set 
    FieldRef<double,3> x     = points.addField<double,3>("x");
    FieldRef<double,3> u     = points.addField<double,3>("u");
    FieldRef<double,3> v     = points.addField<double,3>("v");
    FieldRef<double>   m     = points.addField<double>("m");
    FieldRef<bool>     fixed = points.addField<bool>("fixed");
    FieldRef<int>      idx   = points.addField<int>("idx");

    // The fields of the hinge set 
    FieldRef<double>   theta0   = hinges.addField<double>("theta0");
    FieldRef<double>   he0         = hinges.addField<double>("he0");
    FieldRef<double>   B         = hinges.addField<double>("B");

    // The fields of the face set 
    FieldRef<double,3,2> loc_dir  = faces.addField<double,3,2>("loc_dir");
    FieldRef<double,2,3> precomp_w  = faces.addField<double,2,3>("precomp_w");
    FieldRef<double> area = faces.addField<double>("area");
    FieldRef<double> Nu = faces.addField<double>("Nu");
    FieldRef<double> E = faces.addField<double>("E");
    FieldRef<double> Rho = faces.addField<double>("Rho");

    double timeout = 0.01;

    // create a mesh
    make_mesh(100, 100, 1.0, 1.0, points, hinges, faces);
    compute_loc_dir(faces);
    set_props(faces, hinges);
    
    plot_results("e0.vtk", points, faces);

  // Compile program and bind arguments
  Program program;
  program.loadFile("../shell.sim");

  Function precompute = program.compile("initializeClothPhysics");
  precompute.bind("points",  &points);
  precompute.bind("hinges", &hinges);
  precompute.bind("faces", &faces);
  precompute.bind("timeout", &timeout);
  
  Function timestep   = program.compile("main");
  timestep.bind("points",  &points);
  timestep.bind("hinges", &hinges);
  timestep.bind("faces", &faces);  
  timestep.bind("timeout", &timeout);

  precompute.runSafe();

  char filename[32];

  plot_results("output/shell_000.vtk", points, faces);

  for(int i=1; i<100; i++) {
    timestep.runSafe();       // Run the timestep function
    sprintf(filename, "output/shell_%03d.vtk", i);

    plot_results(filename, points, faces);
  }
    

}