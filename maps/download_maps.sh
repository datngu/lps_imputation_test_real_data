cd maps
for i in {1..22}
do
    wget https://raw.githubusercontent.com/odelaneau/GLIMPSE/master/maps/genetic_maps.b38/chr${i}.b38.gmap.gz
done
cd ..
