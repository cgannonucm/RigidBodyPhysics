 #include <iostream>
#include "evolution.h"


using namespace physics_evolution;

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
    u.add_interaction(id2,id3,&spring2);
    
    //u.add_force(&fg);

    //auto f = u.get_forces()[0];

    for (int i = 0; i < 2000; i++){
        evolve_step(u,0.01);
        std::cout << "(" << u.clock << "," << p2.r[2] - p1.r[2] << ")" << std::endl;
    }
    
    return 0;
};