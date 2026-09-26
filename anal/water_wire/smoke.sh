ff19_hewl=/home/wayyne/cphmd/amber24/ff19/ti/harris_2022/prot/hewl/val-take2/rex-cont

python3 calc_water_wire.py \
    -p ${ff19_hewl}/../../pre/*.prmtop \
    -t ${ff19_hewl}/arex.ph6_0.nc \
    --site-a ':36@OE1,OE2' \
    --site-b ':53@OD1,OD2' \
    --solvent-mask ':WAT' \
    --water-oxygen-mask ':WAT@O' \
    --endpoint-cutoff 3.4 \
    --water-cutoff 3.5 \
    --max-waters 3 \
    --closest-waters 64 \
    --image-anchor '^1' \
    --first-frame-step 510000 \
    --frame-step-interval 10000 \
    --trajin-start 1 \
    --trajin-stop 500 \
    --trajin-stride 1 \
    --jobs 1 \
    --cpptraj-threads 1 \
    --keep-intermediates \
    --workdir ff19_wire_test \
    -o ff19_ph6_wire_test.tsv
