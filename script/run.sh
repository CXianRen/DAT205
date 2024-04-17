set -x
mkdir -p build

cp -r scenes build/
cp  project/*.frag build/project/
cp  project/*.vert build/project/

cd build/project && make && ./project

