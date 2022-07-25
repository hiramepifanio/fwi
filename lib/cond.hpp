float cfl_condition(float dz, float dx, float vel_max);
bool verify_cfl_condition(float dz, float dx, float dt, float vel_max);
float spatial_condition(float dz, float dx, float vel_min);
bool verify_spatial_condition(float dz, float dx, float vel_min, float fp);