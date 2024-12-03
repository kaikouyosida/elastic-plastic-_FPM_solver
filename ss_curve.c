#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include"ss_curve.h"
#include"type.h"

extern Global global;
extern SS_CURVE ss_curve;

double get_hardening_stress(const double equivalent_plastic_strain){
    
    if(ss_curve.num == 0)
        return 1.0e+100;
    
    if(ss_curve.num == 1)
        return ss_curve.stress[0];
    
    if (equivalent_plastic_strain < ss_curve.plastic_strain[0])
    {   
        double sy_i0 = ss_curve.stress[0];
        double sy_i1 = ss_curve.stress[1];

        double dep = ss_curve.plastic_strain[1] - ss_curve.plastic_strain[0];
        double dep_i0 = equivalent_plastic_strain - ss_curve.plastic_strain[0];
        double dep_i1 = ss_curve.plastic_strain[1] - equivalent_plastic_strain;

        return (dep_i1 * sy_i0 + dep_i0 * sy_i1) / dep;
    }
    
    if (equivalent_plastic_strain > ss_curve.plastic_strain[ss_curve.num - 1])
    {   
        double sy_i0 = ss_curve.stress[ss_curve.num - 2];
        double sy_i1 = ss_curve.stress[ss_curve.num - 1];

        double dep = ss_curve.plastic_strain[ss_curve.num - 1] - ss_curve.plastic_strain[ss_curve.num - 2];
        double dep_i0 = equivalent_plastic_strain - ss_curve.plastic_strain[ss_curve.num - 2];
        double dep_i1 = ss_curve.plastic_strain[ss_curve.num - 1] - equivalent_plastic_strain;

        return (dep_i1 * sy_i0 + dep_i0 * sy_i1) / dep;
    }
    
    for (int i = 0; i < ss_curve.num - 1; i++)
        if (equivalent_plastic_strain >= ss_curve.plastic_strain[i]
            && equivalent_plastic_strain <= ss_curve.plastic_strain[i + 1])
        {   
            double sy_i0 = ss_curve.stress[i];
            double sy_i1 = ss_curve.stress[i + 1];
            
            double dep = ss_curve.plastic_strain[i + 1] - ss_curve.plastic_strain[i];
            double dep_i0 = equivalent_plastic_strain - ss_curve.plastic_strain[i];
            double dep_i1 = ss_curve.plastic_strain[i + 1] - equivalent_plastic_strain;
            //printf("%5d  %+15.14e\n", i, (dep_i1 * sy_i0 + dep_i0 * sy_i1) / dep);
            return (dep_i1 * sy_i0 + dep_i0 * sy_i1) / dep;
        }
    #if 0
        printf("Warning: Failed to get stress in ss curve");
    #endif

    return nan("");
}

double get_hardening_modulus(const double equivalent_plastic_strain){

    if (ss_curve.num == 0)
        return 0.0;

    if (ss_curve.num == 1)
        return 0.0;

    if (equivalent_plastic_strain < ss_curve.plastic_strain[0])
    {
        double dsy = ss_curve.stress[1] - ss_curve.stress[0];
        double dep = ss_curve.plastic_strain[1] - ss_curve.plastic_strain[0];

        return dsy / dep;
    }

    if (equivalent_plastic_strain > ss_curve.plastic_strain[ss_curve.num - 1])
    {
        double dsy = ss_curve.stress[ss_curve.num - 1] - ss_curve.stress[ss_curve.num - 2];
        double dep = ss_curve.plastic_strain[ss_curve.num - 1] - ss_curve.plastic_strain[ss_curve.num - 2];

        return dsy / dep;
    }

    for (int i = 0; i < ss_curve.num - 1; i++)
        if (equivalent_plastic_strain >= ss_curve.plastic_strain[i]
            && equivalent_plastic_strain <= ss_curve.plastic_strain[i + 1])
        {
            double dsy = ss_curve.stress[i + 1] - ss_curve.stress[i];
            double dep = ss_curve.plastic_strain[i + 1] - ss_curve.plastic_strain[i];

            return dsy / dep;
        }

#if 0
    printError("Warning: failed to get hardening modulus in SS curve\n");
#endif

    return nan("");
}

