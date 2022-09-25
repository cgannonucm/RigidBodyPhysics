#include <iostream>
#include <fstream>
#include "universe.h"
#include "evolution.h"
#include "constraint_solver.h"


int main()
{
    using namespace std;
    using namespace physics;
    using namespace Eigen; 

    Universe u(1);

    double dt = 0.001;

    Particle p1;
    p1.r = Vector3d(1,0,0);
    p1.v = Vector3d(0,0,0);
    p1.m = 1;

    //Particle p2;
    //p2.r = Vector3d(2,2,2);
    //p2.v = Vector3d(0,0,0);
    //p2.m = 1;
 
    SpringSeries spring_coef(5,0.75);
    InteractionSeries spring_force(&spring_coef);
    

    ForceGravity fg(Vector3d(0,-9.81,0));

    u.add_force(&fg);

    u.add_particle(p1);
    //u.add_particle(p2);

    ConstraintPivot1P cnst(p1,2,0);

    u.add_constraint(&cnst);

    UniverseSolver solver;

    solver.evolve_step(u,dt);



    //u.add_interaction(&spring_force,0,1);
    ofstream outfile;
    outfile.open("out.txt");

    for (int i = 0; i < 100000; i++){
        solver.evolve_step(u,dt);
        auto x0 = u.get_p(0).r[0];
        auto y0 = u.get_p(0).r[1];

        if (i % 200 == 0)
            outfile << "(" << x0 << "," << y0 << ")" << std :: endl;
    }

    outfile.close();





    std :: cout << "Test" << std :: endl;
    


    return 0;
}
 

/*int main() {

    ConstraintSolver solver;

    Particle p1;
    p1.m = 1;
    p1.r = *new Vector3d(1,1,1);
    p1.v = *new Vector3d(0,0,0);

    Universe u;

    vector<Vector3d> dummy_force = {* new Vector3d(0,0,0)}; 

    u.add_particle(&p1);

    ConstraintPivot1P pc1(&p1,2);
    
    u.add_constraint(&pc1);

    solver.solve(dummy_force,u);


    return 0;
};*/

/*
int main() {
    
    Universe u;
    ForceGravity fg;

    SpringSeries spring_series(5,1);
    InteractionSeries spring1(&spring_series);
    InteractionSeries spring2(&spring_series);
    

    Particle p1;
    p1.m = 1;
    p1.r = *new Vector3d(0,0,0);
    p1.v = *new Vector3d(0,0,0);

    Particle p2;
    p2.m = 1;
    p2.r = *new Vector3d(0,0,0.5);
    p2.v = *new Vector3d(0,0,0);

    Particle p3;
    p3.m = 1;
    p3.r = *new Vector3d(0,0,1.5);
    p3.v = *new Vector3d(0,0,0);
    

    auto id1 = u.add_particle(&p1);
    auto id2 = u.add_particle(&p2);
    auto id3 = u.add_particle(&p3);

    u.add_interaction(id1,id2,&spring1);
    
    //u.add_force(&fg);

    //auto f = u.get_forces()[0];

    for (int i = 0; i < 500; i++){
        evolve_step(u,0.01);
        //std::cout << "(" << u.clock << "," << p1.r[2] << ")" << std::endl;
        std::cout << "(" << u.clock << "," << p2.r[2] - p1.r[2] << ")" << std::endl;
    }
    
    return 0;
};


*/