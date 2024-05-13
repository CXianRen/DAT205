set -x
mkdir -p build

cp -r scenes build/
cp  project/*.frag build/project/
cp  project/*.vert build/project/

cp  particle_system/shaders/*.frag build/particle_system/
cp  particle_system/shaders/*.vert build/particle_system/

# cd build/project && make && ./project
cd build/particle_system && make && ./particle_system