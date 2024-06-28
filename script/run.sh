set -x
mkdir -p build

cp -r scenes build/
cp  project/*.frag build/project/
cp  project/*.vert build/project/

cp  src/shaders/*.frag build/src/
cp  src/shaders/*.vert build/src/

# cd build/project && make && ./project
cd build/src && make && ./particle_system