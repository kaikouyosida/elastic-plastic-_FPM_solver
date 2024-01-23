#include<stdio.h>
#include<stdlib.h>
#include"type.h"
#include"fpm.h"
#include"field.h"
#include"internal_force.h"

extern Global global;
extern Option option;

void analize_by_NewtonRapdon(){
    double time;
    double time_increment;
    init_field();
    
    time_increment = option.Delta_time;
    time = time_increment;
    for(int time_step = 0; time_step < option.N_timestep; time_step++){
        for(int iteration_step = 0; iteration_step < 1000; iteration_step++){
            zero_fill_displacement_increments();
            
        }
    }
    
    break_field();
}