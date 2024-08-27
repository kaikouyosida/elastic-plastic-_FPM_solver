//硬化則に基づいた相当応力を計算
double get_hardening_stress(double equivalent_plastic_strain);

//硬化曲線に基づいた硬化係数を計算
double get_hardening_modulus(const double equivalent_plastic_strain);