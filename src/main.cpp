#include <iostream>
#include <fstream>
#include "core/universe.hpp"
#include "core/evolution.hpp"
#include "core/csolver.hpp"
#include "core/logger.hpp"

int main()
{
    USING_STANDARD_NAMESPACES;
    using namespace physics;

    Universe u(3);

    double dt = 0.001;

    Particle p1;
    p1.r = Vector3d(1,0,0);
    p1.v = Vector3d(0,0,0);
    p1.m = 1;

    Particle p2;
    p2.r = Vector3d(2,0,0);
    p2.v = Vector3d(0,-1,0);
    p2.m = 1;

    Particle p3;
    p3.r = Vector3d(3,0,0);
    p3.v = Vector3d(0,0,0);
    p3.m = 1;
 
    //SpringSeries spring_coef(5,0.75);
    //InteractionSeries spring_force(&spring_coef);
    

    ForceGravity fg(Vector3d(0,-9.81,0));

    u.add_force(&fg);

    u.add_particle(p1);
    u.add_particle(p2);
    u.add_particle(p3);

    ConstraintPivot1P cnst(p1,2,0);
    ConstraintPivot2P cnst2(p1,p2,0,1);
    ConstraintPivot2P cnst3(p2,p3,1,2);

    u.add_constraint(&cnst);
    u.add_constraint(&cnst2);
    u.add_constraint(&cnst3);

    UEvolver solver;

    //solver.evolve_step(u,dt);

    double logtstep = 1.0/60.0;
    DefaultLogger logger(u,logtstep);
    solver.loggers.push_back(&logger);

    //u.add_interaction(&spring_force,0,1);


    for (int i = 0; i < 50000; i++){
        solver.evolve_step(u,dt);
    }

    string logout = logger.get_csv_output();

    ofstream outfile;
    outfile.open("out/pendulum3.csv");
    outfile << logout;
    outfile.close();





    std :: cout << "Test" << std :: endl;
    


    return 0;
}